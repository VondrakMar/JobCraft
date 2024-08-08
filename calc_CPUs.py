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
args = parser.parse_args()


hpc_setting = args.hpc 
# Raven setting
if hpc_setting == "raven":
    CPUS_PER_NODE = 72

# Viper setting
elif hpc_setting =="viper":
    CPUS_PER_NODE = 128


submitted_nodes = args.usedN
node_per_submit = args.N 
cpu_per_job = args.n 
if node_per_submit == 1:
    n_cpus = submitted_nodes*CPUS_PER_NODE
    at_the_same_time = n_cpus/cpu_per_job 
    if (not at_the_same_time.is_integer()):
        print(f"The resulting number of jobs at the same time is not integer, you will waste resources. at_the_same_time variable si {at_the_same_time}") 
elif node_per_submit >1:
    # here I am assuming, the full nodes are used
    at_the_same_time = submitted_nodes/node_per_submit
    if (not at_the_same_time.is_integer()):
        print(f"The resulting number of jobs at the same time is not integer, you will waste resources. at_the_same_time variable si {at_the_same_time}") 


# for aims
at_the_same_time = int(at_the_same_time)
print(f"bash prep_aims.sh aimsRun.py temp {at_the_same_time} {64}")
