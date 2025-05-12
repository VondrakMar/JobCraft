import hpc_workflow
import formats
import argparse
import json

parser = argparse.ArgumentParser()
parser.add_argument("hpc",help="Setup which HPC you want to use, current settings are viper and raven",type=str,choices=["viper","raven"])
parser.add_argument("--usedN",help="How many nodes are used in one submission script",type=int)
parser.add_argument("-N",help="How many nodes per job, if this is 1, -n has to be setup",type=int)
parser.add_argument("-n",help="How many cpus is used per job, in srun the value after -n, when using multiple nodes this is ignored for now",type=int)
parser.add_argument("--diff_Ncpu",help="Use this if you use more than 1 node for 1 job, but you for some reason do not want use ALL CPU cores on these nodes for the submission, but want to have a custon srun -n",action="store_true")
parser.add_argument("--prep_submit",help="If true, script will create an header for the basic FHI-Aims job based on provided values",action="store_true")
#parser.add_argument("--job_name","-J",help="Name that will be put into the header of submit file",type=str,default="aims_job")
#parser.add_argument("--strucs",help="File from which DFT folders should be prepared",type=str)
parser.add_argument("--strucs_format",help="Extension of the strucs file",type=str,default="extxyz")
parser.add_argument("--method",help="Name of code you want to use",type=str,choices=["aims","min_ase"],default="aims")
parser.add_argument("--wall_time",help="Wall time in a slurm script, has to be provided in the string form hh:mm:ss",default="1:00:00",type=str)
parser.add_argument("--per_submit",help="How many jobs will run in 1 submit", type=int,default=64)
#parser.add_argument("--use_ase",help="if the ase should be used for dft calculation. Be aware, this is not implemented for all dft code, and ASE has some limitation in certain cases",action="store_true")
################# AIMS INPUT
# parser.add_argument("--aims_basis",choices=["light","intermediate","tight"],type=str)
# parser.add_argument("--aims_species_path",default=None,type=str)
# parser.add_argument("--aims_geometry_lines",nargs='+',default=[])
# parser.add_argument("--aims_atoms_lines",nargs='+',default=[])
# parser.add_argument("--aims_atoms_indeces",nargs='+',default=[])
# parser.add_argument("--aims_kwargs",type=str,default=None)

args = parser.parse_args()
# strucs_ext = formats.ext_to_name(args.strucs_format)
# strucs_format = args.strucs_format
# aims_kwargs = json.loads(args.aims_kwargs)

###########


########
per_file = args.per_submit
hpc_setting = args.hpc
my_job = hpc_workflow.HPC_job(
    usedN = args.usedN,
    N = args.N,
    n = args.n,
    method = args.method,
    hpc_setting = args.hpc,
    path_to_species=None,
    head_temp_name = "head_resubmit",
    submit_file_name = "submit_resubmit",
    paral_file_name = "paral_resubmit",
    diff_Ncpu=args.diff_Ncpu
)
# if args.prep_submit:
my_job.prep_submit_header(wall_time = args.wall_time)
'''
I wrote a restart file only for aims run without ASE (TODO: and I will get rid of all ASE bits for aims here in the future)
'''
my_job.prep_only_aims_submit_restart(
                        prep_name_folders="struc",
                        per_file=per_file)