#!/usr/bin/env python3
"""
Geometry optimisation + NVT (Langevin) MD driven by a MACE model,
optionally with Grimme-D3 dispersion added through ASE's SumCalculator.

Examples
--------
# minimise + 100 ps NVT at 300 K, trajectories -> struc_min.traj / struc_MD.traj
python run_mace.py -s struc.xyz -m MACE-matpes-pbe-omat-ft.model

# with D3(BJ)/PBE dispersion, hotter run, custom trajectory base name
python run_mace.py -s slab.xyz -m my.model --dispersion --xc pbe -T 500 --trj-name slab500

# relax the cell too, then MD
python run_mace.py -s bulk.cif --format cif -m my.model --relax-cell
"""

from __future__ import annotations

import argparse
import inspect
import logging
import sys
import time
from pathlib import Path

import torch
from ase import units
from ase.calculators.mixing import SumCalculator
from ase.filters import FrechetCellFilter
from ase.io import read, write
from ase.io.trajectory import Trajectory
from ase.md import Langevin
from ase.md.md import MDLogger
from ase.md.velocitydistribution import (
    MaxwellBoltzmannDistribution,
    Stationary,
    ZeroRotation,
)
from ase.optimize import QuasiNewton

LOG = logging.getLogger("mace-run")


# --------------------------------------------------------------------------- #
#  logging
# --------------------------------------------------------------------------- #
def setup_logging(logfile: str | None = None, level: str = "INFO") -> None:
    """Console + (optional) file logging with timestamps."""
    fmt = logging.Formatter(
        fmt="%(asctime)s | %(levelname)-7s | %(message)s", datefmt="%H:%M:%S"
    )
    LOG.setLevel(getattr(logging, level.upper()))
    LOG.handlers.clear()

    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    LOG.addHandler(sh)

    if logfile:
        fh = logging.FileHandler(logfile, mode="a")
        fh.setFormatter(fmt)
        LOG.addHandler(fh)
    LOG.propagate = False


class OptLogger:
    """Pretty per-step logger for an ASE optimiser."""

    def __init__(self, opt, target, interval: int = 1, tag: str = "OPT"):
        self.opt = opt          # the optimiser
        self.target = target    # atoms *or* filter the optimiser acts on
        self.interval = interval
        self.tag = tag
        self.t0 = time.time()
        self.header_done = False

    def __call__(self):
        n = self.opt.nsteps
        if n % self.interval:
            return
        f = self.target.get_forces()
        fmax = (f ** 2).sum(axis=1).max() ** 0.5
        e = self.target.get_potential_energy()
        if not self.header_done:
            LOG.info("[%s] %6s %12s %12s %9s", self.tag, "step", "E [eV]",
                     "fmax [eV/A]", "t [s]")
            self.header_done = True
        LOG.info("[%s] %6d %12.6f %12.5f %9.1f", self.tag, n, e, fmax,
                 time.time() - self.t0)


class MDStatusLogger:
    """Per-interval logger for an MD run (energies, instantaneous T, speed)."""

    def __init__(self, dyn, atoms, n_steps: int, interval: int = 100,
                 tag: str = "MD"):
        self.dyn = dyn
        self.atoms = atoms
        self.n_steps = n_steps
        self.interval = interval
        self.tag = tag
        self.t0 = time.time()
        self.header_done = False

    def __call__(self):
        n = self.dyn.nsteps
        a = self.atoms
        epot = a.get_potential_energy()
        ekin = a.get_kinetic_energy()
        ndof = 3 * len(a) - 3  # crude, ignores constraints
        temp = 2.0 * ekin / (ndof * units.kB)
        elapsed = time.time() - self.t0
        speed = n / elapsed if elapsed > 0 else 0.0
        eta = (self.n_steps - n) / speed if speed > 0 else float("nan")

        if not self.header_done:
            LOG.info("[%s] %8s %9s %13s %11s %13s %8s %9s %9s", self.tag,
                     "step", "t [ps]", "Epot [eV]", "Ekin [eV]", "Etot [eV]",
                     "T [K]", "step/s", "ETA [s]")
            self.header_done = True
        LOG.info("[%s] %8d %9.3f %13.5f %11.5f %13.5f %8.1f %9.2f %9.0f",
                 self.tag, n, self.dyn.get_time() / (1000 * units.fs),
                 epot, ekin, epot + ekin, temp, speed, eta)


# --------------------------------------------------------------------------- #
#  calculator
# --------------------------------------------------------------------------- #
def build_calculator(model: str,
                     device: str = "cuda",
                     default_dtype: str = "float64",
                     enable_cueq: bool = False,
                     dispersion: bool = False,
                     xc: str = "pbe",
                     damping: str = "bj",
                     dispersion_cutoff: float = 40.0 * units.Bohr,
                     head: str | None = None):
    from mace.calculators import MACECalculator

    kwargs = dict(model_paths=model, device=device, default_dtype=default_dtype)
    # `enable_cueq` only exists in newer mace-torch versions -> pass it only if
    # the installed signature accepts it.
    sig = inspect.signature(MACECalculator.__init__).parameters
    if head is not None:
        if "head" in sig:
            kwargs["head"] = head
        else:
            LOG.warning("Installed MACECalculator has no 'head' argument "
                        "-- ignoring it.")
    if enable_cueq:
        if "enable_cueq" in sig:
            kwargs["enable_cueq"] = True
        else:
            LOG.warning("Installed MACECalculator has no 'enable_cueq' "
                        "argument -- ignoring it.")
    if "model_paths" not in sig and "model_path" in sig:  # very old mace
        kwargs["model_path"] = kwargs.pop("model_paths")

    LOG.info("Building MACECalculator(model=%s, device=%s, dtype=%s, cueq=%s, "
             "head=%s)", model, device, default_dtype, enable_cueq, head)
    mace_calc = MACECalculator(**kwargs)

    if not dispersion:
        return mace_calc

    try:
        from torch_dftd.torch_dftd3_calculator import TorchDFTD3Calculator
    except ImportError as exc:
        raise ImportError(
            "D3 dispersion requested but 'torch-dftd' is not installed. "
            "Install it with: pip install torch-dftd"
        ) from exc

    LOG.info("Adding TorchDFTD3Calculator(xc=%s, damping=%s, cutoff=%.2f A)",
             xc, damping, dispersion_cutoff)
    d3_calc = TorchDFTD3Calculator(
        device=device,
        damping=damping,
        xc=xc,
        cutoff=dispersion_cutoff,
        dtype={"float64": torch.float64, "float32": torch.float32}.get(default_dtype),
    ) if "dtype" in inspect.signature(TorchDFTD3Calculator.__init__).parameters \
        else TorchDFTD3Calculator(device=device, damping=damping, xc=xc,
                                  cutoff=dispersion_cutoff)

    return SumCalculator([mace_calc, d3_calc])


# --------------------------------------------------------------------------- #
#  drivers
# --------------------------------------------------------------------------- #
def run_minimization(mol, stationary=True, trj_name="min", fmax=0.05,
                     steps=10_000, relax_cell=False, log_interval=1):
    if stationary:
        Stationary(mol)

    target = FrechetCellFilter(mol) if relax_cell else mol
    LOG.info("Starting %s minimisation: fmax=%.3f eV/A, max %d steps -> %s.traj",
             "cell+geometry" if relax_cell else "geometry", fmax, steps, trj_name)

    dyn = QuasiNewton(target, logfile=f"{trj_name}.log")
    traj = Trajectory(f"{trj_name}.traj", "w", mol)
    dyn.attach(traj.write, interval=1)
    dyn.attach(OptLogger(dyn, target, interval=log_interval, tag="MIN"),
               interval=log_interval)

    t0 = time.time()
    converged = dyn.run(fmax=fmax, steps=steps)
    traj.close()

    f = target.get_forces()
    LOG.info("Minimisation finished after %d steps in %.1f s "
             "(converged=%s, E=%.6f eV, fmax=%.5f eV/A)",
             dyn.nsteps, time.time() - t0, bool(converged),
             mol.get_potential_energy(),
             (f ** 2).sum(axis=1).max() ** 0.5)
    return mol


def run_NVT(mol, T=300, time_step=1.0, n_steps=100_000, stationary=True,
            trj_name="MD", constraints=None, init_T=None, friction=0.001,
            traj_interval=50, log_interval=100, zero_rotation=False):
    if constraints is not None:
        mol.set_constraint(constraints)

    if init_T is not None:
        LOG.info("Initialising velocities from Maxwell-Boltzmann at %.1f K", init_T)
        MaxwellBoltzmannDistribution(mol, temperature_K=init_T)
    if stationary:
        Stationary(mol)
    if zero_rotation:
        ZeroRotation(mol)

    if "H" in mol.get_chemical_symbols() and time_step >= 1.0:
        LOG.warning("Time step %.2f fs with hydrogens present -- consider "
                    "<= 0.5 fs or constraining X-H bonds.", time_step)

    LOG.info("Starting NVT Langevin: T=%.1f K, dt=%.2f fs, friction=%.4f, "
             "%d steps (%.1f ps) -> %s.traj",
             T, time_step, friction, n_steps, n_steps * time_step / 1000.0,
             trj_name)

    dyn = Langevin(mol, time_step * units.fs, temperature_K=T, friction=friction)

    traj = Trajectory(f"{trj_name}.traj", "a", mol)
    dyn.attach(traj.write, interval=traj_interval)
    dyn.attach(MDLogger(dyn, mol, f"{trj_name}.log", header=True, stress=False,
                        peratom=False, mode="a"), interval=log_interval)
    dyn.attach(MDStatusLogger(dyn, mol, n_steps, interval=log_interval, tag="MD"),
               interval=log_interval)

    t0 = time.time()
    dyn.run(n_steps)
    traj.close()
    LOG.info("MD finished: %d steps in %.1f s (%.2f steps/s)",
             n_steps, time.time() - t0, n_steps / max(time.time() - t0, 1e-9))
    return mol


# --------------------------------------------------------------------------- #
#  CLI
# --------------------------------------------------------------------------- #
def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="MACE (+optional D3) minimisation and NVT MD.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # structure / calculator
    p.add_argument("-s", "--structure", default="struc.xyz",
                   help="Structure file to read into `mol`.")
    p.add_argument("--format", default=None,
                   help="ASE format of the structure file (default: guess; "
                        "extxyz is assumed for .xyz).")
    p.add_argument("-m", "--model", default="MACE-matpes-pbe-omat-ft.model",
                   help="Path to the MACE model file.")
    p.add_argument("--head", default=None,
                   help="Name of the MACE head to use (multi-head models); "
                        "default: the model's default head.")
    p.add_argument("--device", default="cuda", choices=["cuda", "cpu", "mps"])
    p.add_argument("--dtype", default="float64", choices=["float64", "float32"],
                   help="MACE default_dtype.")
    p.add_argument("--enable-cueq", action="store_true",
                   help="Use cuEquivariance kernels (if supported).")

    # dispersion
    p.add_argument("--dispersion", action="store_true",
                   help="Add Grimme D3 on top of MACE via SumCalculator "
                        "(needs torch-dftd).")
    p.add_argument("--xc", default="pbe", help="Functional for the D3 parameters.")
    p.add_argument("--damping", default="bj", choices=["bj", "zero", "bjm", "zerom"],
                   help="D3 damping function.")
    p.add_argument("--dispersion-cutoff", type=float, default=40.0 * units.Bohr,
                   help="D3 cutoff in Angstrom.")

    # naming
    p.add_argument("--trj-name", default=None,
                   help="Base name for outputs. Default: the structure file "
                        "stem, giving <name>min.traj and <name>MD.traj.")

    # minimisation
    p.add_argument("--fmax", type=float, default=0.01, help="Force convergence.")
    p.add_argument("--min-steps", type=int, default=10000,
                   help="Max optimiser steps.")
    p.add_argument("--relax-cell", action="store_true",
                   help="Relax the cell too (FrechetCellFilter).")
    p.add_argument("--skip-min", action="store_true", help="Skip minimisation.")
    p.add_argument("--freeze-in-min", action="store_true",
                   help="Keep the input file's constraints (frozen atoms) "
                        "active during minimisation too, instead of the "
                        "default of lifting them for minimisation and only "
                        "reapplying them for MD.")

    # MD
    p.add_argument("--n-steps", type=int, default=100000, help="MD steps.")
    p.add_argument("-T", "--temperature", type=float, default=300.0,
                   help="Thermostat temperature in K.")
    p.add_argument("--init-T", type=float, default=None,
                   help="Initial Maxwell-Boltzmann temperature in K "
                        "(default: 0.5 * T).")
    p.add_argument("--timestep", type=float, default=1.0, help="MD time step in fs.")
    p.add_argument("--friction", type=float, default=0.001,
                   help="Langevin friction (ASE units).")
    p.add_argument("--traj-interval", type=int, default=50,
                   help="Write a frame every N MD steps.")
    p.add_argument("--skip-md", action="store_true", help="Skip the MD stage.")

    # misc
    p.add_argument("--log-interval", type=int, default=100,
                   help="Logging interval (MD steps; minimisation uses 1).")
    p.add_argument("--log-level", default="INFO",
                   choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    stem = Path(args.structure).stem
    base = args.trj_name if args.trj_name is not None else stem
    min_name, md_name = f"{base}min", f"{base}MD"

    setup_logging(logfile=f"{base}.run.log", level=args.log_level)
    LOG.info("=" * 78)
    LOG.info("Command: %s", " ".join(sys.argv))

    fmt = args.format
    if fmt is None and Path(args.structure).suffix == ".xyz":
        fmt = "extxyz"
    mol = read(args.structure, format=fmt)
    LOG.info("Read %s: %d atoms, formula %s, pbc=%s",
             args.structure, len(mol), mol.get_chemical_formula(), mol.pbc)

    # Preserve any constraints (e.g. frozen/fixed atoms) defined in the input
    # file. They are lifted for minimisation and reapplied for MD.
    orig_constraints = mol.constraints
    if orig_constraints:
        LOG.info("Found %d constraint(s) in %s (will be lifted for "
                 "minimisation and reapplied for MD).",
                 len(orig_constraints), args.structure)

    mol.calc = build_calculator(
        model=args.model,
        device=args.device,
        default_dtype=args.dtype,
        enable_cueq=args.enable_cueq,
        dispersion=args.dispersion,
        xc=args.xc,
        damping=args.damping,
        dispersion_cutoff=args.dispersion_cutoff,
        head=args.head,
    )
    LOG.info("Initial energy: %.6f eV", mol.get_potential_energy())

    if not args.skip_min:
        if orig_constraints and not args.freeze_in_min:
            LOG.info("Lifting constraints for minimisation.")
            mol.set_constraint()
        run_minimization(mol, fmax=args.fmax, steps=args.min_steps,
                         trj_name=min_name, relax_cell=args.relax_cell)
        write(f"{min_name}.xyz", mol, format="extxyz")
        LOG.info("Wrote relaxed structure to %s.xyz", min_name)

    if not args.skip_md:
        if orig_constraints:
            LOG.info("Reapplying %d constraint(s) for MD.", len(orig_constraints))
        init_T = args.init_T if args.init_T is not None else 0.5 * args.temperature
        run_NVT(mol,
                T=args.temperature,
                init_T=init_T,
                time_step=args.timestep,
                friction=args.friction,
                n_steps=args.n_steps,
                trj_name=md_name,
                constraints=orig_constraints or None,
                traj_interval=args.traj_interval,
                log_interval=args.log_interval)
        write(f"{md_name}_final.xyz", mol, format="extxyz")
        LOG.info("Wrote final MD structure to %s_final.xyz", md_name)

    LOG.info("Done.")


if __name__ == "__main__":
    main()