from mace.calculators import mace_mp
import numpy as np
from ase.calculators.calculator import Calculator, all_changes
from ase import build, units
from ase.md import Langevin
from ase.io.trajectory import Trajectory
from ase.md.velocitydistribution import (
    MaxwellBoltzmannDistribution,
    Stationary,
    ZeroRotation)
from ase.optimize import QuasiNewton, MDMin
from ase.filters import FrechetCellFilter
from ase.io import read,write
from ase.constraints import FixAtoms, Hookean
from ase.calculators.mixing import SumCalculator



def run_minimization(mol,
                    # calculator, # calculator should be attached before
                    stationary=True,
                    trj_name = "bfgs_ls",
                    fmax=0.05,steps=100000000000):
    if stationary:
        Stationary(mol)
    # mol.calc = calculator
    dyn = QuasiNewton(atoms=mol, trajectory=f'{trj_name}.traj')#, restart=f'{trj_name}.pckl')
    dyn.run(fmax=fmax,steps=steps)


def run_minimization_cell(mol,
                    # calculator, # calculator should be attached before
                    stationary=True,
                    trj_name = "bfgs_ls",
                    fmax=0.05,
                    steps=100000000000):
    if stationary:
        Stationary(mol)
    # mol.calc = calculator
    ecf = FrechetCellFilter(mol)
    qn = QuasiNewton(ecf)
    traj = Trajectory(f"{trj_name}.traj","w",mol)
    qn.attach(traj)
    qn.run(fmax=fmax,steps=steps) 
    # dyn = QuasiNewton(atoms=mol, trajectory=f'{trj_name}.traj')#, restart=f'{trj_name}.pckl') qn.run(fmax=fmax,steps=steps)


def run_NVT(mol,
            T=300,
            time_step = 1.0,
            n_steps=100000,
            stationary=True,
            trj_name = "mdNVT",
            constraints=None,
            init_T=None):
    if constraints is not None:
        mol.set_constraint(constraints)
    if init_T is not None:
        MaxwellBoltzmannDistribution(mol, init_T * units.kB)
    if stationary:g
        Stationary(mol)
    if "H" in mol.symbols and time_step >= 1.0:
        print("Be aware you are running time step larger than what I would do for structure with hydrogens")
    dyn = Langevin(mol, time_step * units.fs, T * units.kB, 0.001)
    traj = Trajectory(trj_name + '.traj', 'a', mol)
    dyn.attach(traj.write, interval=50)
    dyn.run(n_steps)

from mace.calculators import MACECalculator    
if __name__== "__main__":
    atoms = read("struc.xyz",format="extxyz")
    mace_calc = mace_mp(model="MACE-matpes-pbe-omat-ft.model",device="cuda",dispersion=False,enable_cueq=True,default_dtype="float64")   
    #mace_calc = MACECalculator(model_path='CuOrun3_naive.model', device='cuda')
    mol.calc = mace_calc
    run_minimization(mol,fmax=0.01,steps=10000,trj_name=f"gasMin")
    run_NVT(mol,
            n_steps=100000,
            T=300,
            init_T= 0.5*300,
            trj_name=f"gasMD")
