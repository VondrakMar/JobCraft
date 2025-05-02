from mace.calculators import mace_mp
from ase import build, units
from ase.md import Langevin
from ase.io.trajectory import Trajectory
from ase.md.velocitydistribution import (
    MaxwellBoltzmannDistribution,
    Stationary,
    ZeroRotation)
from ase.optimize import QuasiNewton, MDMin
from ase.io import read,write
from ase.constraints import FixAtoms, Hookean

def run_minimization(mol,calculator,stationary=True,trj_name = "bfgs_ls",fmax=0.05):
    if stationary:
        Stationary(mol)
    mol.calc = calculator
    dyn = QuasiNewton(atoms=mol, trajectory=f'{trj_name}.traj', restart=f'{trj_name}.pckl')
    dyn.run(fmax=fmax)

def run_NVT(mol,
            calculator,
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
    mol.calc = calculator
    dyn = Langevin(mol, time_step * units.fs, T * units.kB, 0.001)
    traj = Trajectory(trj_name + '.traj', 'a', mol)
    dyn.attach(traj.write, interval=50)
    dyn.run(n_steps)





if __name__== "__main__":
    c = FixAtoms(indices=[117,127])
    atoms = read("test.xyz",format="extxyz")
    macemp = mace_mp(model="https://github.com/ACEsuit/mace-foundations/releases/download/mace_matpes_0/MACE-matpes-pbe-omat-ft.model",device="cuda",dispersion=False)
    run_minimization(atoms,macemp)
    run_NVT(atoms,macemp,init_T = 600)


