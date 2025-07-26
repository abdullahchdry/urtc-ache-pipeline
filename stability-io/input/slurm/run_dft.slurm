#!/bin/bash
#SBATCH --job-name=medoid_dft
#SBATCH --array=1,4,5,8,10,11      # only medoids 1,4,5,8,10,11
#SBATCH --nodes=1                  # all threads on one node
#SBATCH --ntasks=1                 # one ORCA process per task
#SBATCH --cpus-per-task=16        # matches %pal nprocs in .inp files
#SBATCH --gres=gpu:1               # to access the GPU partition
#SBATCH --partition=gpuA100x4      # Delta GPU queue (change as needed) 
#SBATCH --account=bcpv-delta-gpu   # change based on your account
#SBATCH --mem=40G                  # ~2 GB per core
#SBATCH --time=24:00:00            # 24 hour limit
#SBATCH --output=compound_medoid_%A_%a.out
#SBATCH --error=compound_medoid_%A_%a.err

module load orca/5.0.1
module load openmpi               # change based on hpc cluster specifications

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMPI_MCA_rmaps_base_oversubscribe=1

echo "===== MEDOID ${SLURM_ARRAY_TASK_ID} START ====="
echo "Host: $(hostname)"
echo "Date: $(date)"

# point at the input in its folder, but write output here
INP="compound_medoid_${SLURM_ARRAY_TASK_ID}/compound_medoid_${SLURM_ARRAY_TASK_ID}.inp"
OUT="compound_medoid_${SLURM_ARRAY_TASK_ID}.out"

if [[ ! -f "$INP" ]]; then
  echo "ERROR: Input $INP missing, skipping."
  exit 1
fi

ORCABIN=$(which orca)
echo "Running: $ORCABIN $INP"
$ORCABIN $INP > "$OUT"

echo "===== MEDOID ${SLURM_ARRAY_TASK_ID} DONE at $(date) ====="
