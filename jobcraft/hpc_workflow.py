import os
import my_presets
import file_creation
'''
the script prep_aims.sh is used as bash prep_aims.sh aimsRun.py temp at_the_same_time in_one_file
herein I am calculating how many of at_the_same_time and in_one_file should be used based on sbatch setting
this is usable for everything
this script assumes you are using GNU parallel
WARNING: This code is tested only with my version of ASE. It probably will still work with 3.22.1 version of ASE 
from the main repository, but the newest one will definetly screem in problems
example of usage 
python hpc_workflow.py raven --usedN 4 -N 1 -n 4 --prep_submit --strucs waters.xyz --method="aims" --aims_basis="tight" --wall_time="1:00:00"
'''

class HPC_job():
    def __init__(self,
                 usedN = 1,
                 N = 1,
                 n = 1,
                 method = "aims",
                 hpc_setting = "raven"):
        # Raven settings
        if hpc_setting == "raven":
            self.NTASKS_PER_NODE = 72 # this has to be change
            self.CPUS_PER_NODE = 72
            self.MEMORY_PER_NODE_GB = 240 
            self.MEMORY_PER_NODE_MB = 240000 
            self.PRESETS_FOR_HEADER = my_presets.aims_for_raven
            self.AIMS_EXEC = my_presets.aims_exec_raven
            self.AIMS_PATH = my_presets.aims_path_raven
            self.AIMS_SPECIEC = my_presets.aims_species_raven
        # Viper setting
        elif hpc_setting =="viper":
            self.NTASKS_PER_NODE = 128 # this has to be change, it works only for AIMS settings
            self.CPUS_PER_NODE = 128
            self.MEMORY_PER_NODE_GB = 480 # I took this from the viper website from the table
            self.MEMORY_PER_NODE_MB = 480000 
            self.PRESETS_FOR_HEADER = my_presets.aims_for_viper

        self.submitted_nodes = usedN
        self.node_per_job = N 
        self.cpu_per_job = n 
        self.method = method
        if self.node_per_job == 1:
            n_cpus = self.submitted_nodes*self.CPUS_PER_NODE
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
        header_file = open("header_file.temp","w")
        header_file.write("#!/bin/bash -l\n")
        header_file.write("#SBATCH -o ./tjob.out.%j\n") # stardart output file
        header_file.write("#SBATCH -e ./tjob.err.%j\n") # error output file
        header_file.write("#SBATCH -D ./\n") # working directory 
        # header_file.write(f"#SBATCH -J {job_name}\n") # name of the file
        header_file.write("\n")
        header_file.write(f"#SBATCH --nodes={self.submitted_nodes}\n")
        header_file.write(f"#SBATCH --ntasks-per-node={self.NTASKS_PER_NODE}\n")
        header_file.write(f"#SBATCH --cpus-per-task=1\n")
        header_file.write(f"#SBATCH --time={wall_time}\n")
        header_file.write("\n")
        header_file.write(f"{self.PRESETS_FOR_HEADER}")
        #header_file.write("#SBATCH --cpus-per-task=1\n")
        header_file.close()

    def prep_ase_file(self,
                      aims_basis,
                      aims_run_file):
        if self.method == "aims":
            file_creation.prep_aims_ase_file(final_name=aims_run_file,
                                     aims_command=f"{self.AIMS_PATH}{self.AIMS_EXEC}",
                                     aims_species=f"{self.AIMS_PATH}{self.AIMS_SPECIEC}{aims_basis}")
    
    def prep_folders(self,
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
