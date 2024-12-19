#!/bin/env bash
#
#$ -q gpu.q
#$ -N HIV_HS   # CHANGE THIS -- any name you want
#$ -cwd
###$ -l h_rt=24:00:00
#$ -l h_rt=2:00:00
#$ -l mem_free=64G
#$ -l scratch=50G
#$ -l compute_cap=80,gpu_mem=40G
#$ -j y
#$ -o /wynton/group/krogan/mgordon/projects/112624_MGordon_AF3_pipeline/output/jobLogs/$JOB_NAME-$JOB_ID-$TASK_ID.log

#$ -t 4          # CHANGE THIS - match numbers in jobTable  ## job array with xx tasks

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
	--jobTable=/wynton/group/krogan/mgordon/projects/112624_MGordon_AF3_pipeline/data/2024_12_16.hivHu.newJobtable.csv \
	--job_id=$taskID \
	--master_fasta=/wynton/group/krogan/mgordon/projects/112624_MGordon_AF3_pipeline/data/121624_HIV.HS.uniprot.seqs.fa \
	--output_dir=/wynton/group/krogan/mgordon/projects/112624_MGordon_AF3_pipeline/output/HIV_HS \
	--nSeeds=5

# ./run_alphafold3.py \
# 	--jobTable=/wynton/group/krogan/mgordon/projects/112624_MGordon_AF3_pipeline/data/2024_12_16.hivHu.newJobtable.csv \
# 	--job_id=$taskID \
# 	--master_fasta=/wynton/group/krogan/mgordon/projects/112624_MGordon_AF3_pipeline/data/121624_HIV.HS.uniprot.seqs.fa \
# 	--output_dir=/wynton/group/krogan/mgordon/projects/112624_MGordon_AF3_pipeline/output/HIV_HS \
# 	--nSeeds=5

# for testing 
./run_alphafold3.py \
	--jobTable=/wynton/group/krogan/mgordon/projects/112624_MGordon_AF3_pipeline/data/pten.jobTable.txt \
	--job_id=$taskID \
	--master_fasta=/wynton/group/krogan/mgordon/projects/112624_MGordon_AF3_pipeline/data/pten_preys.fa \
	--output_dir=/wynton/group/krogan/mgordon/projects/112624_MGordon_AF3_pipeline/output/HIV_HS \
	--nSeeds=5

t1=$(date --rfc-3339=seconds)
echo "Duration: $t0 -- $t1"