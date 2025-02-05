#!/bin/env bash
#
#$ -q gpu.q
#$ -N af3   # CHANGE THIS -- any name you want
#$ -cwd
###$ -l h_rt=24:00:00
#$ -l h_rt=2:00:00
#$ -l mem_free=64G
#$ -l scratch=50G
#$ -l compute_cap=80,gpu_mem=40G
#$ -j y
#$ -o ./jobLogs/$JOB_NAME-$JOB_ID-$TASK_ID.log

#$ -t 1-47     # CHANGE THIS - match numbers in jobTable  ## job array with xx tasks

# if not running with sge task array, set to 5
taskID=${SGE_TASK_ID:-5}

# necessary for the GetContactsPAE.R
module load CBI
module load r
#
# Compute cap for A100 GPU is 8.0 (40 or 80 GB), for A40 GPU is 8.6 (48 GB).
#
t0=$(date --rfc-3339=seconds)

echo "QUEUE: $QUEUE"
echo "SGE_GPU: $SGE_GPU"
export CUDA_VISIBLE_DEVICES=$SGE_GPU

echo ./run_alphafold3.py \
	--jobTable=AlphaFoldJobList.csv \
	--job_id=$taskID \
	--master_fasta=masterFasta.fasta \
	--output_dir=./output \
	--nSeeds=5

./run_alphafold3.py \
	--jobTable=AlphaFoldJobList.csv \
	--job_id=$taskID \
	--master_fasta=masterFasta.fasta \
	--output_dir=./output \
	--nSeeds=5

t1=$(date --rfc-3339=seconds)
echo "Duration: $t0 -- $t1"

