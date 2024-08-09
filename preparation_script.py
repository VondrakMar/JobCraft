import hpc_workflow
import formats
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("hpc",help="Setup which HPC you want to use, current settings are viper and raven",type=str,choices=["viper","raven"])
parser.add_argument("--usedN",help="How many nodes are used in one submission script",type=int)
parser.add_argument("-N",help="How many nodes per job, if this is 1, -n has to be setup",type=int)
parser.add_argument("-n",help="How many cpus is used per job, in srun the value after -n, when using multiple nodes this is ignored for now",type=int)
parser.add_argument("--prep_submit",help="If true, script will create an header for the basic FHI-Aims job based on provided values",action="store_true")
parser.add_argument("--job_name","-J",help="Name that will be put into the header of submit file",type=str,default="aims_job")
parser.add_argument("--strucs",help="File from which DFT folders should be prepared",type=str)
parser.add_argument("--strucs_format",help="Extension of the strucs file",type=str,default="extxyz")
parser.add_argument("--method",help="Name of code you want to use",type=str,choices=["aims"])
parser.add_argument("--wall_time",help="Wall time in a slurm script, has to be provided in the string form hh:mm:ss",default="1:00:00",type=str)
parser.add_argument("--per_submit",help="How many jobs will run in 1 submit", type=int,default=64)

################# AIMS INPUT
parser.add_argument("--aims_basis",choices=["light","intermediate","tight"])


args = parser.parse_args()
strucs_ext = formats.ext_to_name(args.strucs_format)
strucs_format = args.strucs_format

###########
aims_run_file = "aims_run.py"
########
per_file = args.per_submit
hpc_setting = args.hpc

my_job = hpc_workflow.HPC_job(
    usedN = args.usedN,
    N = args.N,
    n = args.n,
    method = args.method
)
if args.prep_submit:
    my_job.prep_submit_header(wall_time = args.wall_time)
my_job.prep_ase_file(aims_basis=args.aims_basis,aims_run_file=aims_run_file)
my_job.prep_folders(strucs=args.strucs,
                    aims_run_file=aims_run_file,
                     strucs_format=strucs_format,
                     strucs_ext=strucs_ext,
                     per_file=per_file)