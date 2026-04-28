import ase.io
from pathlib import Path
from ase.visualize import view

from jobcraft.structures.interfaces import solvate_slab
from jobcraft.aseMD.runMD import run_minimization

runMD = True
data_path = "../../data/"
water = ase.io.read(f"{data_path}water.pdb")
slab = ase.io.read(f"{data_path}Cu65O6.cif")
slab.cell[2] = slab.cell[2]*4
# dmf = ase.io.read(f"{data_path}dmf.xyz") # uncomment to test with different solvent 

slab_solv =solvate_slab(slab,water,slab_element="Cu",N_mols=30,sphere_factor = 1.4)
if runMD:
    from mace.calculators import mace_mp    
    macemp = mace_mp(model="mh-1",head="omat_r2scan")
    slab_solv.calc = macemp
    run_minimization(slab_solv,fmax=0.05,steps=100,trj_name=f"struc_min")
#view(slab_solv)
