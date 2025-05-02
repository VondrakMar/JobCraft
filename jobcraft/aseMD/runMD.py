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
from ase.calculators.mixing import SumCalculator


class HarmonicRestraint(Calculator):
    implemented_properties = ['energy', 'forces']

    def __init__(self, i, j, k, r0):
        super().__init__()
        self.i = i # atoms for constrain 
        self.j = j # atoms for constrain
        self.k = k # in eV/A^2
        self.r0 = r0 # distance between them

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        pos = atoms.get_positions()
        rij = pos[self.j] - pos[self.i]
        dist = np.linalg.norm(rij)
        direction = rij / dist if dist != 0 else np.zeros(3)
        delta = dist - self.r0
        force = -self.k * delta * direction

        forces = np.zeros_like(pos)
        forces[self.i] -= force
        forces[self.j] += force

        energy = 0.5 * self.k * delta**2

        self.results = {
            'energy': energy,
            'forces': forces,
        }


def run_minimization(mol,
                    # calculator, # calculator should be attached before
                    stationary=True,
                    trj_name = "bfgs_ls",
                    fmax=0.05):
    if stationary:
        Stationary(mol)
    mol.calc = calculator
    dyn = QuasiNewton(atoms=mol, trajectory=f'{trj_name}.traj', restart=f'{trj_name}.pckl')
    dyn.run(fmax=fmax)

def run_NVT(mol,
            # calculator, # calculator should be attached before 
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
    # mol.calc = calculator
    dyn = Langevin(mol, time_step * units.fs, T * units.kB, 0.001)
    traj = Trajectory(trj_name + '.traj', 'a', mol)
    dyn.attach(traj.write, interval=50)
    dyn.run(n_steps)


if __name__== "__main__":
    c = FixAtoms(indices=[117,127])
    atoms = read("test.xyz",format="extxyz")
    macemp = mace_mp(model="https://github.com/ACEsuit/mace-foundations/releases/download/mace_matpes_0/MACE-matpes-pbe-omat-ft.model",device="cuda",dispersion=False)
    calc_constrain = HarmonicRestraint(i=117, j=127, k=5.0, r0=3.0)
    atoms.calc = sumcalculator([macemp, restraint])
    run_minimization(atoms)#,macemp)
    run_NVT(atoms,macemp,init_T = 600)


