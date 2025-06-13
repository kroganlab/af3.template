#!/bin/env bash
#
#$ -q member.q
#$ -N af3.msaOnly
#$ -cwd
#$ -l h_rt=48:00:00 #extend as long as needed
#$ -l mem_free=10G
#s -pe smp 12
##$ -l scratch=80G
#$ -j y
#$ -o ./jobLogs/$JOB_NAME-$JOB_ID-$TASK_ID.log

#$ -t 1-16           # CHANGE THIS - match numbers in jobTable  ## job array with xx tasks

# if not running with sge task array, set to 5
taskID=${SGE_TASK_ID:-5}

# necessary for the GetContactsPAE.R
module load CBI
module load r

t0=$(date --rfc-3339=seconds)

echo "QUEUE: $QUEUE"

echo ./run_alphafold3.py \
	--jobTable=./pten.jobTable.txt \
	--job_id=$taskID \
	--master_fasta=./pten_preys.fa \
	--output_dir=./output \
	--nSeeds=5 \
	--run_inference=False
	
echo ./run_alphafold3.py \
	--jobTable=./pten.jobTable.txt \
	--job_id=$taskID \
	--master_fasta=./pten_preys.fa \
	--output_dir=./output \
	--nSeeds=5 \
	--run_inference=False

t1=$(date --rfc-3339=seconds)
echo "Duration: $t0 -- $t1"