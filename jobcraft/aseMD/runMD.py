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
    qn.run(fmax=fmax,steps) 
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
    if stationary:
        Stationary(mol)
    if "H" in mol.symbols and time_step >= 1.0:
        print("Be aware you are running time step larger than what I would do for structure with hydrogens")
    dyn = Langevin(mol, time_step * units.fs, T * units.kB, 0.001)
    traj = Trajectory(trj_name + '.traj', 'a', mol)
    dyn.attach(traj.write, interval=50)
    dyn.run(n_steps)

    
if __name__== "__main__":
    # atoms = read("def1Vs_part3.xyz@:",format="extxyz")
    atoms = read("def1Vs_part1.xyzMD2.traj@::50")
    macemp = mace_mp(model="MACE-matpes-pbe-omat-ft.model",device="cuda",dispersion=False,enable_cueq=True,default_dtype="float64")   
    for id_mol,mol in enumerate(atoms): 
        mol.calc = macemp
        run_minimization_cell(mol,fmax=0.001,steps=1000,trj_name=f"struc_min{id_mol}")
        run_NVT(mol,
                n_steps=50000,
                T=400,
                init_T= 2*400,
                trj_name=f"strucMD{id_mol}")
                     
