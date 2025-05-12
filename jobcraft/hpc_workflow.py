import os
import my_presets
import file_creation
import aims.aims_input
import numpy as np
'''
the script prep_aims.sh is used as bash prep_aims.sh aimsRun.py temp at_the_same_time in_one_file
herein I am calculating how many of at_the_same_time and in_one_file should be used based on sbatch setting
this is usable for everything
this script assumes you are using GNU parallel
WARNING: This code is tested only with my version of ASE. It probably will still work with 3.22.1 version of ASE 
from the main repository, but the newest one will definetly screem in problems
'''

class HPC_job():
    def __init__(self,
                 usedN = 1,
                 N = 1,
                 n = 1,
                 method = "aims",
                 hpc_setting = "raven",
                 path_to_species=None,
                 head_temp_name = "header_file",
                 submit_file_name = "submit_file",
                 paral_file_name = "paral_file",
                 diff_Ncpu=False):
        # Raven settings
        self.method = method
        self.head_temp_name = head_temp_name
        self.hpc_setting = hpc_setting
        self.diff_Ncpu = diff_Ncpu
        self.submit_file_name = submit_file_name
        self.paral_file_name = paral_file_name
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
        if self.diff_Ncpu and self.node_per_job > 1:
            self.cpu_per_job = n 
        else:
            self.cpu_per_job = self.CPUS_PER_NODE_HW
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
        This function is universal for all methods, only thing dependable on used calc method is PRESETS_FOR_HEADER. 
        '''
        header_file = open(f"{self.head_temp_name}.temp","w")
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
        with open(f"{self.head_temp_name}.temp","r") as head_file:
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
                slurm_file = open(f"{self.submit_file_name}{count}.sl","w")
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
                          atoms_lines=[],
                          atoms_indeces=[],
                          aims_basis="light",
                          per_file=64,
                          all_control_same = True,
                          aims_kwargs_dict=None):
        if isinstance(atoms_lines, str):
            atoms_lines = [atoms_lines]
        assert len(atoms_indeces) == len(atoms_lines) or len(atoms_lines) == 1, "atoms_indeces and atoms_lines are different lenght"
        if len(atoms_indeces) > 0:
            atoms_indeces = [int(tmp) for tmp in atoms_indeces]
            # Indeces has to be sorted, otherwise the lines will move and it won't be placed in the correct spot
            sorted_indices = np.argsort(-np.array(np.array(atoms_indeces)))  # minus for descending sort            
            atoms_indeces = [atoms_indeces[sorted_indice] for sorted_indice in sorted_indices]
            if len(atoms_lines) > 1:
                atoms_lines = [atoms_lines[sorted_indice] for sorted_indice in sorted_indices]
        with open(f"{self.head_temp_name}.temp","r") as head_file:
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
                slurm_file = open(f"{self.submit_file_name}{count}.sl","w")
                slurm_file.write(head_data)
                slurm_file.write(f"parallel --delay 0.2 --joblog task.log --progress -j {int(self.at_the_same_time)} < {self.paral_file_name}{count}")
                paral_file = open(f"paral_file{count}","w")
            dir_name = f"struc{id_mol:0{counting_digits}}/"
            struc_file_name = f"struc{id_mol:0{counting_digits}}{strucs_ext}" 
            ase.io.write(struc_file_name,mol,format=f"{strucs_format}")
            ase.io.write("temp.in",mol,format=f"aims")
            #################
            with open('temp.in', 'r') as file:
                temp_geometry = file.readlines()
            appended_header = 5 # probably be aware if ASE will change number of lines it putting in the geometry.in file 
            for geometry_line in geometry_lines:
                temp_geometry.insert(appended_header, f'{geometry_line}\n')
                appended_header += 1
            if len(atoms_lines) == 1:
                cur_line = atoms_lines[0]
                for atoms_indx in atoms_indeces:
                    temp_geometry.insert(atoms_indx+appended_header+1, f'{cur_line}\n') # + 1 because the keywords inside of the geometry.in are applied on the previous line
            with open('geometry.in', 'w') as file:
                file.writelines(temp_geometry)

            ################
            os.mkdir(dir_name)
            shutil.move(struc_file_name,dir_name)
            shutil.move("geometry.in",dir_name)
            if not all_control_same:
                aims.aims_input.prep_aims_file(mol,aims_species)
            shutil.copy("control.in",dir_name)
            # paral_file.write(f"cd {dir_name}; srun -N {self.node_per_job} -n {self.cpu_per_job_to_srun} {aims_command} >> aims.out; python -c \"import sys; from jobcraft.aims.aims_output import read_aims_output; import ase.io; from jobcraft.file_creation import save_results_to_xyz; res = read_aims_output(mol_file_name=f'{{sys.argv[1]}}.xyz', properties=['energy', 'forces', 'hirshfeld']); mol = ase.io.read(f'{{sys.argv[1]}}.xyz', format='extxyz'); save_results_to_xyz(mol, res)\" {dir_name[:-1]}\n")
            # paral_file.write(f"cd {dir_name}; srun -N {self.node_per_job} -n {self.cpu_per_job_to_srun} {aims_command}\n; python3 -c 'import sys; from jobcraft.aims.aims_output import read_aims_output; import ase.io; from jobcraft.file_creation import save_results_to_xyz; res = read_aims_output(mol_file_name="struc00100.xyz", properties=["energy", "forces", "hirshfeld"]); mol = ase.io.read(f"{sys.argv[1]}.xyz", format="extxyz"); save_results_to_xyz(mol, res)' {dir_name[:-1]}.xyz")
            paral_file.write(f"cd {dir_name}; srun -N {self.node_per_job} -n {self.cpu_per_job_to_srun} {aims_command} >> aims.out\n")

    def prep_only_aims_submit_restart(self,
                                      prep_name_folders="struc",
                                        per_file=64):
        '''
        This is the code that will prep submission for unfinished jobs
        At this moment I am distiqusing if the job should be restart only by 
        "Have a nice day" not in lines[line] 
        or 
        "scf_solver: SCF cycle not converged"
        TODO: allow to resubmit SCF not converged structures
        
        TODO: This should allows to change srun setting, but will used previous control.in and geometry.in 
        '''

        #if isinstance(atoms_lines, str):
        #    atoms_lines = [atoms_lines]
        len_prep = len(prep_name_folders)
        aims_command=f"{self.AIMS_EXEC}"
        with open(f"{self.head_temp_name}.temp","r") as head_file:
            head_data= head_file.read()
        import subprocess
        if os.path.exists("end_of_aimsout_files"):
            os.remove("end_of_aimsout_files")
        command = "for a in */; do echo $a >> end_of_aimsout_files; tail -2 ${a}aims.out | head -1 >> end_of_aimsout_files ; done"
        subprocess.run(command, shell=True, capture_output=True, text=True) # this is way faster than open each file in python
        lines = open("end_of_aimsout_files","r").readlines()
        prev_line = lines[0]
        to_restart = []
        for line in range(1,len(lines)):
            if prep_name_folders in lines[line] and prep_name_folders in prev_line:
                to_restart.append(prev_line) 
            if prep_name_folders not in lines[line] and ("Have a nice day" not in lines[line] and "scf_solver: SCF cycle not converged" not in lines[line]) and prep_name_folders in prev_line:
                to_restart.append(prev_line)
            prev_line = lines[line]
        count = -1
        for id_mol,mol_name in enumerate(to_restart):
            if (mol_name[-1]=="\n"):
                mol_name = mol_name[:-1]
                # id_mol = int(mol_name[len_prep:-2])
            # else:
            # id_folder = int(mol_name[len_prep:])
            if id_mol%per_file == 0:
                if id_mol != 0:
                    slurm_file.close()
                    paral_file.close()
                count+=1
                slurm_file = open(f"{self.submit_file_name}{count}.sl","w")
                slurm_file.write(head_data)
                slurm_file.write(f"parallel --delay 0.2 --joblog task.log --progress -j {int(self.at_the_same_time)} < {self.paral_file_name}{count}")
                paral_file = open(f"{self.paral_file_name}{count}","w")
            paral_file.write(f"cd {mol_name}; mv aims.out aims_tmp.out ; srun -N {self.node_per_job} -n {self.cpu_per_job_to_srun} {aims_command} >> aims.out\n")
