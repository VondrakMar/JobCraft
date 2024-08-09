aims_command="raven"
aims_species="raven/species"

import sys
from ase import Atoms
from ase.calculators.aims import Aims
from ase.io import read
mols = read(f"{sys.argv[1]}", format="extxyz")

aims_kwargs = {
    'xc': 'pbe0',
    'relativistic': ("atomic_zora", "scalar"),
    'compute_forces': True,
#    'charge': 1,
#    'fixed_spin_moment': 1,
    'override_illconditioning': True,
}

calc = Aims(output=["hirshfeld","hartree_multipoles"],
        command=f"srun -N {sys.argv[2]} -n {sys.argv[3]} {aims_command}",
            species_dir=f"{aims_species}",#"/u/mvondrak/software/fhi-aims.231212_1/species_defaults/defaults_2020/tight/",
            **aims_kwargs)
mols.calc = calc
mols.get_potential_energy()

