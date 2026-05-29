#!/bin/env bash
#
#SBATCH --partition=gpu
#SBATCH --job-name=af3_pipeline
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
##SBATCH --gres=gpu:nvidia_h100_nvl:1
#SBATCH --output=./jobLogs/%x-%j-%a.log
#SBATCH --array=1-4

# # CHANGE THESE PATHS AS NEEDED
MASTER_FASTA=./masterFasta.fasta
JOB_TABLE=./AlphaFoldJobList.csv
OUTPUT_DIR=./output


# if not running with sge task array, set to 5
taskID=${SLURM_ARRAY_TASK_ID:-5}

# necessary for the GetContactsPAE.R
# necessary for the GetContactsPAE.R
module purge
conda deactivate 2>/dev/null || true
module load CBI
module load r
module load nvidia/cuda/12.9.1
module load apptainer/1.4.1
#
# Compute cap for A100 GPU is 8.0 (40 or 80 GB), for A40 GPU is 8.6 (48 GB).
#
t0=$(date --rfc-3339=seconds)

# optional: check GPU
echo "SLURM_JOB_ID: $SLURM_JOB_ID"
echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"
echo "CUDA_VISIBLE_DEVICES: $CUDA_VISIBLE_DEVICES"
nvidia-smi

echo ./run_alphafold3.py \
	--jobTable=$JOB_TABLE \
	--job_id=$taskID \
	--master_fasta=$MASTER_FASTA \
	--output_dir=$OUTPUT_DIR \
	--ipsae_pae_threshold=10 \
	--ipsae_dist_threshold=10 \
	--nSeeds=5

./run_alphafold3.py \
	--jobTable=$JOB_TABLE \
	--job_id=$taskID \
	--master_fasta=$MASTER_FASTA \
	--output_dir=$OUTPUT_DIR \
	--ipsae_pae_threshold=10 \
	--ipsae_dist_threshold=10 \
	--nSeeds=5

t1=$(date --rfc-3339=seconds)
echo "Duration: $t0 -- $t1"
