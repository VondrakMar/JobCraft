#!/bin/bash
#ASC --vanilla
#SBATCH -p zen4_0768
#SBATCH --qos zen4_0768
#SBATCH --ntasks-per-node=40
#SBATCH --threads-per-core=1
#SBATCH --time=12:00:00
#SBATCH --job-name=vasp_batch
#SBATCH --output=vasp_batch_%j.out
#SBATCH --error=vasp_batch_%j.err
#SBATCH --nodelist=n3020-011
#SBATCH --mem-per-cpu=4G
 
module --force purge
module load ASC/2023.06
module load VASP
export OMP_NUM_THREADS=1
 
VASP_EXEC="vasp_std"
NPROCS=40
 
mpirun --bind-to none -n "$NPROCS" "$VASP_EXEC"
