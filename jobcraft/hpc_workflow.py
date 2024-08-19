import os
import my_presets
import file_creation
import aims.aims_input
'''
the script prep_aims.sh is used as bash prep_aims.sh aimsRun.py temp at_the_same_time in_one_file
herein I am calculating how many of at_the_same_time and in_one_file should be used based on sbatch setting
this is usable for everything
this script assumes you are using GNU parallel
WARNING: This code is tested only with my version of ASE. It probably will still work with 3.22.1 version of ASE 
from the main repository, but the newest one will definetly screem in problems
example of usage 
python /path/to/jobcraft/preparation_script.py raven --usedN 4 -N 1 -n 4 --prep_submit --strucs waters.xyz --method="aims" --aims_basis="tight" --wall_time="1:00:00"
'''

class HPC_job():
    def __init__(self,
                 usedN = 1,
                 N = 1,
                 n = 1,
                 method = "aims",
                 hpc_setting = "raven",
                 path_to_species=None):
        # Raven settings
        self.method = method
        self.hpc_setting = hpc_setting
        if self.method == "aims":
            if self.hpc_setting == "raven":
                self.PRESETS_FOR_HEADER = aims.aims_input.aims_for_raven
                self.AIMS_PATH = aims.aims_input.aims_path_raven
                self.AIMS_EXEC = f"{self.AIMS_PATH}{aims.aims_input.aims_exec_raven}"
                if path_to_species == None:
                    self.AIMS_SPECIEC = f"{self.AIMS_PATH}{aims.aims_input.aims_species_raven}"
                else:
                    self.AIMS_SPECIEC = path_to_species
            elif self.hpc_setting == "viper":
                self.PRESETS_FOR_HEADER = aims.aims_input.aims_for_viper
                self.AIMS_PATH = aims.aims_input.aims_path_viper
                self.AIMS_EXEC = f"{self.AIMS_PATH}{aims.aims_input.aims_exec_viper}"
                if path_to_species == None:
                    self.AIMS_SPECIEC = f"{self.AIMS_PATH}{aims.aims_input.aims_species_viper}"
                else:
                    self.AIMS_SPECIEC = path_to_species


        if hpc_setting == "raven":
            self.NTASKS_PER_NODE = 72 # this has to be change
            self.CPUS_PER_NODE = 1
            self.CPUS_PER_NODE_HW = 72 # I fuck up and this has to be here together with the setting before until I will unfuck it
            self.MEMORY_PER_NODE_GB = 240 
            self.MEMORY_PER_NODE_MB = 240000 

        # Viper setting
        elif hpc_setting =="viper":
            self.NTASKS_PER_NODE = 128 # this has to be change, it works only for AIMS settings
            self.CPUS_PER_NODE = 1
            self.CPUS_PER_NODE_HW = 128 # I fuck up and this has to be here together with the setting before until I will unfuck it
            self.MEMORY_PER_NODE_GB = 480 # I took this from the viper website from the table
            self.MEMORY_PER_NODE_MB = 480000 
            self.PRESETS_FOR_HEADER = my_presets.aims_for_viper

        self.submitted_nodes = usedN
        self.node_per_job = N 
        self.cpu_per_job = n 
        if self.node_per_job == 1:
            n_cpus = self.submitted_nodes*self.CPUS_PER_NODE_HW
            self.at_the_same_time = n_cpus/self.cpu_per_job 
            self.cpu_per_job_to_srun = self.cpu_per_job
            if (not self.at_the_same_time.is_integer()):
                print(f"The resulting number of jobs at the same time is not integer, you will waste resources. at_the_same_time variable si {at_the_same_time}") 
        elif self.node_per_job >1:
            # here I am assuming, the full nodes are used
            self.at_the_same_time = self.submitted_nodes/self.node_per_job
            self.cpu_per_job_to_srun = self.cpu_per_job*self.node_per_job
            if (not self.at_the_same_time.is_integer()):
                print(f"The resulting number of jobs at the same time is not integer, you will waste resources. at_the_same_time variable si {at_the_same_time}") 

    def prep_submit_header(self,
                           wall_time="1:00:00"):
        '''
        This function is universal for all methods, only thing dependable on used cacl method is PRESETS_FOR_HEADER. 
        '''
        header_file = open("header_file.temp","w")
        header_file.write("#!/bin/bash -l\n")
        header_file.write("#SBATCH -o ./tjob.out.%j\n") # stardart output file
        header_file.write("#SBATCH -e ./tjob.err.%j\n") # error output file
        header_file.write("#SBATCH -D ./\n") # working directory 
        # header_file.write(f"#SBATCH -J {job_name}\n") # name of the file
        header_file.write("\n")
        header_file.write(f"#SBATCH --nodes={self.submitted_nodes}\n")
        header_file.write(f"#SBATCH --ntasks-per-node={self.NTASKS_PER_NODE}\n")
        header_file.write(f"#SBATCH --cpus-per-task={self.CPUS_PER_NODE}\n")
        header_file.write(f"#SBATCH --time={wall_time}\n")
        header_file.write("\n")
        header_file.write(f"{self.PRESETS_FOR_HEADER}")
        #header_file.write("#SBATCH --cpus-per-task=1\n")
        header_file.close()

    def prep_ase_file(self,
                      aims_basis,
                      aims_run_file):
        if self.method == "aims":
            aims.aims_input.prep_aims_ase_file(final_name=aims_run_file,
                                     aims_command=f"{self.AIMS_EXEC}",
                                     aims_species=f"{self.AIMS_SPECIEC}{aims_basis}")

    def prep_aims_ase_folders(self,
                     strucs,
                     aims_run_file,
                     strucs_format,
                     strucs_ext,
                     per_file):
        with open("header_file.temp","r") as head_file:
            head_data= head_file.read()
        import ase.io
        import shutil
        mols = ase.io.read(f"{strucs}@:",format=f"{strucs_format}")
        counting_digits = len(str(len(mols)))+1
        prev_sub_mol = 0
        count = -1
        for id_mol,mol in enumerate(mols):
            if id_mol%per_file == 0:
                if id_mol != 0:
                    slurm_file.close()
                    paral_file.close()
                count+=1
                slurm_file = open(f"submit_file{count}.sl","w")
                slurm_file.write(head_data)
                slurm_file.write(f"parallel --delay 0.2 --joblog task.log --progress -j {self.at_the_same_time} < paral_file{count}")
                paral_file = open(f"paral_file{count}","w")
            dir_name = f"struc{id_mol:0{counting_digits}}/"
            struc_file_name = f"struc{id_mol:0{counting_digits}}{strucs_ext}" 
            ase.io.write(struc_file_name,mol,format=f"{strucs_format}")
            os.mkdir(dir_name)
            shutil.move(struc_file_name,dir_name)
            shutil.copy(aims_run_file,dir_name)
            paral_file.write(f"cd {dir_name}; python {aims_run_file} {struc_file_name} {self.node_per_job} {self.cpu_per_job_to_srun}\n")


    def prep_aims_folders(self,
                          strucs,
                          strucs_format,
                          strucs_ext,
                          geometry_lines=[],
                          aims_basis="light",
                          per_file=64,
                          all_control_same = True,
                          aims_kwargs_dict=None):
        with open("header_file.temp","r") as head_file:
            head_data= head_file.read()
        import ase.io
        import shutil
        mols = ase.io.read(f"{strucs}@:",format=f"{strucs_format}")
        aims_command=f"{self.AIMS_EXEC}"
        aims_species=f"{self.AIMS_SPECIEC}{aims_basis}"
        counting_digits = len(str(len(mols)))+1
        prev_sub_mol = 0
        count = -1
        if all_control_same:
            aims.aims_input.prep_aims_file(mols[0],aims_species,aims_kwargs_dict)
        for id_mol,mol in enumerate(mols):
            if id_mol%per_file == 0:
                if id_mol != 0:
                    slurm_file.close()
                    paral_file.close()
                count+=1
                slurm_file = open(f"submit_file{count}.sl","w")
                slurm_file.write(head_data)
                slurm_file.write(f"parallel --delay 0.2 --joblog task.log --progress -j {int(self.at_the_same_time)} < paral_file{count}")
                paral_file = open(f"paral_file{count}","w")
            dir_name = f"struc{id_mol:0{counting_digits}}/"
            struc_file_name = f"struc{id_mol:0{counting_digits}}{strucs_ext}" 
            ase.io.write(struc_file_name,mol,format=f"{strucs_format}")
            ase.io.write("temp.in",mol,format=f"aims")
            #################
            with open('temp.in', 'r') as file:
                temp_control = file.readlines()
            to_which_line = 5 # proablby be aware if ASE will change number of lines it putting in the geometry.in file 
            for geometry_line in geometry_lines:
                temp_control.insert(to_which_line, f'{geometry_line}\n')
                to_which_line += 1
            with open('geometry.in', 'w') as file:
                file.writelines(temp_control)

            ################
            os.mkdir(dir_name)
            shutil.move(struc_file_name,dir_name)
            shutil.move("geometry.in",dir_name)
            if not all_control_same:
                aims.aims_input.prep_aims_file(mol,aims_species)
            shutil.copy("control.in",dir_name)
            paral_file.write(f"cd {dir_name}; srun -N {self.node_per_job} -n {self.cpu_per_job_to_srun} {aims_command} >> aims.out; python -c \"import sys; from jobcraft.aims.aims_output import read_aims_output; import ase.io; from jobcraft.file_creation import save_results_to_xyz; res = read_aims_output(mol_file_name=f'{{sys.argv[1]}}.xyz', properties=['energy', 'forces', 'hirshfeld']); mol = ase.io.read(f'{{sys.argv[1]}}.xyz', format='extxyz'); save_results_to_xyz(mol, res)\" {dir_name[:-1]}\n")
            # paral_file.write(f"cd {dir_name}; srun -N {self.node_per_job} -n {self.cpu_per_job_to_srun} {aims_command}\n; python3 -c 'import sys; from jobcraft.aims.aims_output import read_aims_output; import ase.io; from jobcraft.file_creation import save_results_to_xyz; res = read_aims_output(mol_file_name="struc00100.xyz", properties=["energy", "forces", "hirshfeld"]); mol = ase.io.read(f"{sys.argv[1]}.xyz", format="extxyz"); save_results_to_xyz(mol, res)' {dir_name[:-1]}.xyz")
