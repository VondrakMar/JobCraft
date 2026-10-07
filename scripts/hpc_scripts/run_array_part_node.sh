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
module load  ASC/2023.06
module load VASP
export OMP_NUM_THREADS=1

FOLDERS_PER_JOB=10
PARENT_DIR="$SLURM_SUBMIT_DIR"
PREFIX="struc"
VASP_EXEC="vasp_std"
NPROCS=40

cd "$PARENT_DIR" || exit 1

mapfile -t ALL_FOLDERS < <(ls -d ${PREFIX}* 2>/dev/null | sort)
TOTAL=${#ALL_FOLDERS[@]}

if [ "$TOTAL" -eq 0 ]; then
    echo "ERROR: No folders matching ${PREFIX}* found in $PARENT_DIR"
    exit 1
fi

START=$(( SLURM_ARRAY_TASK_ID * FOLDERS_PER_JOB ))
END=$(( START + FOLDERS_PER_JOB - 1 ))
if [ "$END" -ge "$TOTAL" ]; then
    END=$(( TOTAL - 1 ))
fi
for i in $(seq "$START" "$END"); do
    FOLDER="${ALL_FOLDERS[$i]}"
    if [ -z "$FOLDER" ]; then continue; fi

    cd "$PARENT_DIR/$FOLDER" || { echo "Cannot cd into $FOLDER"; continue; }

    for f in INCAR POSCAR KPOINTS POTCAR; do
        if [ ! -f "$f" ]; then
            echo "    WARNING: $f missing in $FOLDER — skipping"
            cd "$PARENT_DIR" || exit 1
            continue 2
        fi
    done

    mpirun --bind-to none -n "$NPROCS" "$VASP_EXEC"

    EXIT_CODE=$?

    if [ $EXIT_CODE -eq 0 ]; then
        echo "    Finished $FOLDER successfully at $(date)"
    else
        echo "    VASP failed in $FOLDER (exit $EXIT_CODE)"
    fi

    cd "$PARENT_DIR" || exit 1
done
