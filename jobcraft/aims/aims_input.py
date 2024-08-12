import numpy as np
import ase.io
import os

aims_for_raven = '''
export OMP_NUM_THREADS=1
export MKL_DYNAMIC=FALSE
export MKL_NUM_THREADS=1
module load parallel/201807
module load anaconda/3/2020.02
module load intel/21.4.0 mkl/2021.4 impi/2021.4
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:${INTEL_HOME}/compiler/2021.4.0/linux/compiler/lib/intel64_lin/"
'''
aims_path_raven = '/u/mvondrak/software/fhi-aims.231212_1/'
aims_exec_raven="build/aims.231212_1.scalapack.mpi.x"
aims_species_raven='species_defaults/defaults_2020/'


aims_for_viper = '''
export OMP_NUM_THREADS=1
ulimit -s unlimited
module load parallel/201807
module load anaconda/3/2023.03
module load gcc/13 openmpi/4.1 mkl/2024.0
'''
aims_exec_viper = 'build.gnu_mkl_ompi/aims.231212_1.scalapack.mpi.x'
aims_path_viper = '/u/mvondrak/software/fhi-aims.231212_1/'
aims_species_viper = 'species_defaults/defaults_2020/'


aims_path_linux = '/home/mvondrak/software/fhi-aims.231212_1/species_defaults/defaults_2020' 



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

