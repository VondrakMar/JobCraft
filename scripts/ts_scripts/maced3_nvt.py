import os

os.environ["TORCH_USE_CUDA_DSA"] = "1"
os.environ["CUDA_LAUNCH_BLOCKING"] = "1"

from pathlib import Path

import ase
import ase.io
import torch

import torch_sim as ts
from torch_sim.io import atoms_to_state
from torch_sim.telemetry import configure_logging, get_logger
from torch_sim.trajectory import TrajectoryReporter
from torch_SSPD.calculators.modified_maced3 import MaceD3Model
from torch_SSPD.utils.d3utils import extract_dftd3_parameters
import numpy as np

configure_logging(log_file="maced3_nvt.log")
log = get_logger(name="maced3_nvt")
XC_D3_PARAMS = {
    "r2scan": {"a1": 0.49484001, "s8": 0.78981345, "a2": 5.73083694, "s6": 1.0},
    "rpbe": {"a1": 0.1820, "s8": 0.8318, "a2": 4.0094, "s6": 1.0},
    "pbe": {"a1": 0.4289, "s8": 0.7875, "a2": 4.4407, "s6": 1.0},
}

MACE_HEAD = None
SEED = int(os.environ.get("SEED", "42"))
TEMPERATURE = float(os.environ.get("TEMPERATURE", "700"))  # K, same in all directions
XC = "rpbe"

MODEL_PATH = "trainer_it10_0.model"
INPUT_STRUCTURE = "struc.xyz"

N_MEMBERS = 50
N_STEPS_NVT = 5000
N_STEPS_MIN = 1000
TIME_STEP = 0.001  # ps (standard Langevin, not overdamped)
STATE_FREQUENCY = 5
D3_CUTOFF = 20.0
MIN_FMAX = 0.05

DTYPE = torch.float64

############## vibed freezing layer
def get_bottom_layer_indices(atoms, n_layers=2, symbol="Cu", tol=0.5):
    positions = atoms.get_positions()
    symbols = np.array(atoms.get_chemical_symbols())
    elem_indices = np.where(symbols == symbol)[0]

    if len(elem_indices) == 0:
        raise ValueError(f"No atoms with symbol '{symbol}' found in structure")

    z = positions[elem_indices, 2]
    order = np.argsort(z)
    sorted_indices = elem_indices[order]
    sorted_z = z[order]
    layer_id = np.zeros(len(sorted_z), dtype=int)
    for i in range(1, len(sorted_z)):
        if sorted_z[i] - sorted_z[i - 1] > tol:
            layer_id[i:] += 1

    layer_id = np.zeros(len(sorted_z), dtype=int)
    current = 0
    for i in range(1, len(sorted_z)):
        if sorted_z[i] - sorted_z[i - 1] > tol:
            current += 1
        layer_id[i] = current

    n_found = layer_id.max() + 1
    if n_layers > n_found:
        raise ValueError(
            f"Requested {n_layers} layers but only found {n_found} distinct "
            f"'{symbol}' layers with tol={tol} Å. Check tol or your structure."
        )

    bottom_mask = layer_id < n_layers
    frozen = sorted(sorted_indices[bottom_mask].tolist())

    log.info(
        "Identified %d '%s' layers; freezing bottom %d layer(s) -> %d atoms "
        "(local indices %s)",
        n_found, symbol, n_layers, len(frozen), frozen,
    )
    return frozen

def freeze_bottom_layers(state, atoms, n_layers=2, symbol="Cu", tol=0.5):
    from torch_sim.constraints import FixAtoms

    local_frozen = get_bottom_layer_indices(
        atoms, n_layers=n_layers, symbol=symbol, tol=tol
    )
    n_atoms_per_system = len(atoms)

    local = torch.as_tensor(local_frozen, dtype=torch.long, device=state.device)
    offsets = (
        torch.arange(state.n_systems, dtype=torch.long, device=state.device)
        * n_atoms_per_system
    )
    global_indices = (offsets.unsqueeze(1) + local.unsqueeze(0)).flatten()

    state.constraints.append(FixAtoms(atom_idx=global_indices))

    log.info(
        "Froze bottom %d '%s' layer(s) across all %d systems: %d atoms/system, "
        "%d atoms total",
        n_layers, symbol, state.n_systems, len(local_frozen), len(global_indices),
    )
    return state
##############
def set_seed(seed, device):
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    generator = torch.Generator(device=device)
    generator.manual_seed(seed)
    log.info("Random seed set to %d on device %s", seed, device)
    return generator


def load_state(path, n_members, device, dtype):
    input_struct = ase.io.read(path, format="extxyz")
    input_struct.info["total_charge"] = 0
    input_struct.info["spin"] = 0
    log.info(
        "Loaded structure '%s' with %d atoms; replicating into %d members",
        path,
        len(input_struct),
        n_members,
    )

    atoms_list = [input_struct.copy() for _ in range(n_members)]
    atomic_numbers = torch.Tensor(
        [atom.get_atomic_numbers() for atom in atoms_list]
    ).ravel()

    state = atoms_to_state(
        atoms_list,
        device,
        dtype,
        system_extras_map={"spin": "spin", "charge": "total_charge"},
        atom_extras_map={},
    )
    return input_struct, state, atomic_numbers


def get_d3_xc_params(xc):
    try:
        return XC_D3_PARAMS[xc]
    except KeyError:
        raise ValueError(f"Unknown xc functional: {xc!r}") from None


def build_d3_params(param_file):
    log.info("Extracting DFT-D3 reference parameters")
    params = extract_dftd3_parameters()
    torch.save(params, param_file)
    log.info("Saved DFT-D3 parameters to '%s'", param_file)
    return params


def build_maced3_model(model_path, head, state, atomic_numbers, device, dtype, xc,
                       d3_params, d3_cutoff):
    xc_params = get_d3_xc_params(xc)
    log.info(
        "Loading MACE-D3 model from '%s' (head=%s, xc=%s, d3_cutoff=%g, params=%s)",
        model_path,
        head,
        xc,
        d3_cutoff,
        xc_params,
    )
    return MaceD3Model(
        model=model_path,
        device=device,
        head=head,
        dtype=dtype,
        #system_idx=state.system_idx,
        #atomic_numbers=atomic_numbers,
        d3_cutoff=d3_cutoff,
        d3_params=d3_params,
        **xc_params,
    )


def make_reporter(prefix, n_members, frequency):
    files = [f"traj/{prefix}_{i}.h5md" for i in range(n_members)]
    return TrajectoryReporter(files, state_frequency=frequency)



def run_nvt(state, model, temperature, n_steps, time_step, reporter):
    log.info(
        "Starting NVT Langevin dynamics: %d steps, dt=%g ps, T=%g K",
        n_steps, time_step, temperature,
    )
    md_state = ts.integrate(
        system=state,
        model=model,
        integrator=ts.Integrator.nvt_langevin,
        n_steps=n_steps,
        temperature=temperature,
        timestep=time_step,
        trajectory_reporter=reporter,
        autobatcher=False,
    )
    log.info("NVT dynamics complete")
    return md_state


def run_minimization(state, model, n_steps, reporter, fmax):
    log.info(
        "Starting FIRE minimization: max_steps=%d, fmax=%g eV/A",
        n_steps,
        fmax,
    )
    final_state = ts.optimize(
        system=state,
        model=model,
        optimizer=ts.Optimizer.fire,
        convergence_fn=ts.generate_force_convergence_fn(fmax),
        max_steps=n_steps,
        trajectory_reporter=reporter,
        autobatcher=False,
    )
    log.info("Minimization complete")
    return final_state



def main():
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    log.info("=" * 70)
    log.info("NVT Langevin Simulation")
    log.info("head=%s seed=%s T=%s xc=%s", MACE_HEAD, SEED, TEMPERATURE, XC)
    log.info("=" * 70)
    log.info(
        "Run config: N_members=%d, nvt_steps=%d, dt=%g, T=%g K, dtype=%s",
        N_MEMBERS, N_STEPS_NVT, TIME_STEP, TEMPERATURE, DTYPE,
    )

    set_seed(SEED, device)

    input_struct, state, atomic_numbers = load_state(
        INPUT_STRUCTURE, N_MEMBERS, device, DTYPE
    )

    d3_params = build_d3_params(Path.cwd() / "dftd3_parameters.pt")
    model = build_maced3_model(
        MODEL_PATH, MACE_HEAD, state, atomic_numbers, device, DTYPE, XC,
        d3_params, D3_CUTOFF,
    )

    nvt_reporter = make_reporter("nvt", N_MEMBERS, STATE_FREQUENCY)
    premin_reporter = make_reporter("premin", N_MEMBERS, STATE_FREQUENCY)

    state = freeze_bottom_layers(state, input_struct, n_layers=2, symbol="Cu")

    premin_state = run_minimization(
        state, model, N_STEPS_MIN, premin_reporter, MIN_FMAX
    )
    run_nvt(premin_state, model, TEMPERATURE, N_STEPS_NVT, TIME_STEP, nvt_reporter)

    log.info("All stages finished successfully")


if __name__ == "__main__":
    main()
