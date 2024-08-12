import numpy as np
import ase.io
import os


def prep_aims_ase_file(base_file=None,
                       final_name="runAims.py",
                       aims_command="raven",
                       aims_species="raven/species/"):
    if base_file == None:
        this_file_path = os.path.dirname(os.path.realpath(__file__))
        base_file = os.path.join(this_file_path,"../templates/aims_ase_template")
        with open(base_file, 'r') as base_file:
            base_lines = base_file.read()
    else:
        with open(base_file, 'r') as base_file:
            base_lines = base_file.read()
    with open(f'{final_name}', 'w') as final_file:
        final_file.write(f'aims_command="{aims_command}"\n')
        final_file.write(f'aims_species="{aims_species}"\n')
        final_file.write("\n"+base_lines)

    
def prep_aims_file(mol = None,
                   aims_species="raven/species/"):
    from ase.calculators.aims import Aims
    # for writing basic input file for aims I am using ASE.
    aims_kwargs = {
        'xc': 'pbe0',
        'relativistic': ("atomic_zora", "scalar"),
        'compute_forces': True,
        'override_illconditioning': True,
    }
    aims_outputs = ["hirshfeld","hartree_multipoles"]

    calc = Aims(output=aims_outputs,
            command=None,
            species_dir=f"{aims_species}",
            **aims_kwargs)

    calc.write_input(mol,"control.in")

