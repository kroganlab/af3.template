#!/bin/env bash
#
##$ -q member.q
#$ -N af3.msaOnly
#$ -cwd
#$ -l h_rt=24:00:00
#$ -l mem_free=124G
#$ -l scratch=80G
#$ -l compute_cap=80
#$ -j y
#$ -o ./jobLogs/$JOB_NAME-$JOB_ID-$TASK_ID.log

#$ -t 1-161           # CHANGE THIS - match numbers in jobTable  ## job array with xx tasks

# if not running with sge task array, set to 5
taskID=${SGE_TASK_ID:-5}

module load CBI
module load r

t0=$(date --rfc-3339=seconds)

echo "QUEUE: $QUEUE"

echo ./run_alphafold3.py \
	--jobTable=./pten.jobTable.txt \
	--job_id=$taskID \
	--master_fasta=./pten_preys.fa \
	--output_dir=./outDir \
	--nSeeds=5 \
	--run_inference=false \
	--model_preset=multimer

./run_alphafold3.py \
	--jobTable=./pten.jobTable.txt \
	--job_id=$taskID \
	--master_fasta=./pten_preys.fa \
	--output_dir=./outDir \
	--nSeeds=5 \
	--run_inference=false \
	--model_preset=multimer


t1=$(date --rfc-3339=seconds)
echo "Duration: $t0 -- $t1"