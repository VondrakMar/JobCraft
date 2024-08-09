import os
import my_presets
import file_creation
'''
the script prep_aims.sh is used as bash prep_aims.sh aimsRun.py temp at_the_same_time in_one_file
herein I am calculating how many of at_the_same_time and in_one_file should be used based on sbatch setting
this is usable for everything
this script assumes you are using GNU parallel
'''
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("hpc",help="Setup which HPC you want to use, current settings are viper and raven",type=str,choices=["viper","raven"])
parser.add_argument("--usedN",help="How many nodes are used in one submission script",type=int)
parser.add_argument("-N",help="How many nodes per job, if this is 1, -n has to be setup",type=int)
parser.add_argument("-n",help="How many cpus is used per job, in srun the value after -n, when using multiple nodes this is ignored for now",type=int)
parser.add_argument("--prep_submit",help="If true, script will create an header for the basic FHI-Aims job based on provided values",action="store_true")
parser.add_argument("--job_name","-J",help="Name that will be put into the header of submit file",type=str,default="aims_job")
parser.add_argument("--strucs",help="File from which DFT folders should be prepared",type=str)
parser.add_argument("--strucs_ext",help="Extension of the strucs file",type=str,default="extxyz")
parser.add_argument("--method",help="Name of code you want to use",type=str,choices=["aims"])
parser.add_argument("--aims_basis",choices=["light","intermediate","tight"])
parser.add_argument("--wall_time",help="Wall time in a slurm script, has to be provided in the string form hh:mm:ss",type=str)

args = parser.parse_args()

###########
per_file = 64
extension = ".xyz"
aims_run_file = "aims_run.py"
########3
hpc_setting = args.hpc 
# Raven setting
if hpc_setting == "raven":
    CPUS_PER_NODE = 72
    MEMORY_PER_NODE_GB = 240 
    MEMORY_PER_NODE_MB = 240000 
    PRESETS_FOR_HEADER = my_presets.aims_for_raven
    AIMS_EXEC = my_presets.aims_exec_raven
    AIMS_PATH = my_presets.aims_path_raven
    AIMS_SPECIEC = my_presets.aims_species_raven
# Viper setting
elif hpc_setting =="viper":
    CPUS_PER_NODE = 128
    MEMORY_PER_NODE_GB = 480 # I took this from the viper website from the table
    MEMORY_PER_NODE_MB = 480000 
    PRESETS_FOR_HEADER = my_presets.aims_for_viper


submitted_nodes = args.usedN
node_per_job = args.N 
cpu_per_job = args.n 
if node_per_job == 1:
    n_cpus = submitted_nodes*CPUS_PER_NODE
    at_the_same_time = n_cpus/cpu_per_job 
    cpu_per_job_to_srun = cpu_per_job
    if (not at_the_same_time.is_integer()):
        print(f"The resulting number of jobs at the same time is not integer, you will waste resources. at_the_same_time variable si {at_the_same_time}") 
elif node_per_job >1:
    # here I am assuming, the full nodes are used
    at_the_same_time = submitted_nodes/node_per_job
    cpu_per_job_to_srun = cpu_per_job*node_per_job
    if (not at_the_same_time.is_integer()):
        print(f"The resulting number of jobs at the same time is not integer, you will waste resources. at_the_same_time variable si {at_the_same_time}") 


# for aims
at_the_same_time = int(at_the_same_time)
print(f"bash prep_aims.sh aimsRun.py temp {at_the_same_time} {per_file}")
if args.prep_submit:
    job_name=args.job_name
    file_creation.prep_submit_header(submitted_nodes,
                                     CPUS_PER_NODE,
                                     PRESETS_FOR_HEADER,
                                     wall_time = args.wall_time)

if args.method == "aims":
    file_creation.prep_aims_ase_file(final_name=aims_run_file,
                                     aims_command=f"{AIMS_PATH}{AIMS_EXEC}",
                                     aims_species=f"{AIMS_PATH}{AIMS_SPECIEC}{args.aims_basis}")



prep_folders = True
if prep_folders:
    with open("header_file.temp","r") as head_file:
        head_data= head_file.read()
    import ase.io
    import shutil
    mols = ase.io.read(f"{args.strucs}@:",format=f"{args.strucs_ext}")
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
            slurm_file.write(f"parallel --delay 0.2 --joblog task.log --progress -j {at_the_same_time} < paral_file{count}")
            paral_file = open(f"paral_file{count}","w")
        dir_name = f"struc{id_mol:0{counting_digits}}/"
        struc_file_name = f"struc{id_mol:0{counting_digits}}{extension}" 
        ase.io.write(struc_file_name,mol,format=f"{args.strucs_ext}")
        os.mkdir(dir_name)
        shutil.move(struc_file_name,dir_name)
        shutil.copy(aims_run_file,dir_name)
        paral_file.write(f"cd {dir_name}; python {aims_run_file} {struc_file_name} {node_per_job} {cpu_per_job_to_srun}\n")
