from ase.io import read,write
from ase.visualize import view
import numpy as np
from ase.build import molecule
from ase import units

def make_water(density, super_cell=[3, 3, 3]):
    """
    Generates a supercell of water molecules with a desired density.
    Density in g/cm^3!!!
    This is copied from GAP tutorial: https://libatoms.github.io/GAP/gap_fitting_tutorial.html
    """
    h2o = molecule('H2O')
    a = np.cbrt((sum(h2o.get_masses()) * units.m ** 3 * 1E-6 ) / (density * units.mol))
    h2o.set_cell((a, a, a))
    h2o.set_pbc((True, True, True))
    return h2o.repeat(super_cell)


from ase import Atoms
from ase.data import atomic_masses
from ase.units import kJ, mol, Angstrom
import numpy as np


NA = 6.02214076e23  

def density_g_cm3(atoms: Atoms) -> float:
    masses = atoms.get_masses()  # amu per atom
    mass_total_amu = masses.sum()
    mass_g = mass_total_amu / NA  # since 1 mol of 1 amu = 1 g
    vol_A3 = atoms.get_volume()       # Å^3
    vol_cm3 = vol_A3 * 1e-24
    return mass_g / vol_cm3

if __name__ == "__main__":
    water = make_water(1.0, [3, 3, 3])
    rho = density_g_cm3(water)
    print(f"ρ(water) = {rho:.3f} g/cm³")
    
    a = 5.43
    si = Atoms("Si8", cell=[a, a, a], pbc=True)
    rho = density_g_cm3(si)
    print(f"ρ(Si) = {rho:.3f} g/cm³")

    h2o = Atoms("H2O", positions=[[0,0,0],[0.96,0,0],[0,0.96,0]],cell=[10,10,10], pbc=True)
    rho = density_g_cm3(h2o)
    print(f"ρ(H2O box) = {rho:.3f} g/cm³")
