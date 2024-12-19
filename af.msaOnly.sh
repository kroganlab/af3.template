#!/bin/env bash
#
#$ -q !gpu.q
#$ -N af3.msaOnly
#$ -cwd
#$ -l h_rt=24:00:00
#$ -l mem_free=124G
#$ -l scratch=80G
#$ -l compute_cap=80
#$ -j y
#$ -o /wynton/group/krogan/mgordon/projects/112624_MGordon_AF3_pipeline/output/jobLogs/$JOB_NAME-$JOB_ID-$TASK_ID.log

#$ -t 1-16           # CHANGE THIS - match numbers in jobTable  ## job array with xx tasks

# if not running with sge task array, set to 5
taskID=${SGE_TASK_ID:-5}

# necessary for the GetContactsPAE.R
module load CBI
module load r

t0=$(date --rfc-3339=seconds)

echo "QUEUE: $QUEUE"

echo ./run_alphafold3.py \
	--jobTable=/wynton/group/krogan/mgordon/projects/112624_MGordon_AF3_pipeline/data/2024_12_17_af.genMSA.jobtable.csv \
	--job_id=$taskID \
	--master_fasta=/wynton/group/krogan/mgordon/projects/112624_MGordon_AF3_pipeline/data/121624_HIV.HS.uniprot.seqs.fa \
	--output_dir=/wynton/group/krogan/mgordon/projects/112624_MGordon_AF3_pipeline/output/HIV_HS \
	--nSeeds=5 \
	--run_inference=False
	

./run_alphafold3.py \
	--jobTable=/wynton/group/krogan/mgordon/projects/112624_MGordon_AF3_pipeline/data/2024_12_17_af.genMSA.jobtable.csv \
	--job_id=$taskID \
	--master_fasta=/wynton/group/krogan/mgordon/projects/112624_MGordon_AF3_pipeline/data/121624_HIV.HS.uniprot.seqs.fa \
	--output_dir=/wynton/group/krogan/mgordon/projects/112624_MGordon_AF3_pipeline/output/HIV_HS \
	--nSeeds=5 \
	--run_inference=False 

t1=$(date --rfc-3339=seconds)
echo "Duration: $t0 -- $t1"