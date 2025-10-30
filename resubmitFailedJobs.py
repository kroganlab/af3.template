#!/usr/bin/python3
#
#$ -S /usr/bin/python3
#$ -q gpu.q
#$ -N alphafold 
#$ -cwd
#$ -l h_rt=2:00:00
#$ -l mem_free=64G
#$ -l scratch=50G
#$ -l compute_cap=80,gpu_mem=40G
###$ -pe smp 2

import os
import glob
import re
import argparse
from fileinput import FileInput

"""Script to detect and re AlphaFold jobs using the af3.template setup."""

parser = argparse.ArgumentParser(description='Resubmit failed AlphaFold3 jobs')

parser.add_argument(
    '--work_dir', type = str, required = True,
    help = 'working directory where af3.template is located')

parser.add_argument(
    '--jobs_script', required=True, default= None,
    help = "AF jobs submission script used in initial run")

parser.add_argument(
    '--new_jobTable', required=False, default= "resubmit.AlphaFoldJobList.csv",
    help = "Name of new job table to create. csv file with columns ID, seq1.name, seq2.name, seq3.name etc..")

parser.add_argument(
    '--new_runtime', required=True, default= "24:00:00",
    help = "New runtime in SGE format HH:MM:SS for resubmitted jobs")  

args = parser.parse_args()

failed_runs = glob.glob(os.path.join(args.work_dir, "*.msaOut*"))

print(f"Found {len(failed_runs)} incomplete runs in {args.work_dir}\nResubmitting jobs with runtime increased to {args.new_runtime}")

with open(os.path.join(args.work_dir,args.new_jobTable), mode="w", encoding="utf-8") as rerun_jobs:
    for i,msa in enumerate(failed_runs):
        msa = msa.strip()
        ppi = os.path.basename(msa).split(".msaOut")[0]
        rerun_jobs.write(f"{i+1},{ppi.replace('__',',')}\n")

print(f"Copying af3 job submission script...\ncp {os.path.join(args.work_dir, args.jobs_script)} {os.path.join(args.work_dir, 'resubmit_jobs.sh')}") 
#os.popen(f"cp {os.path.join(args.work_dir, args.jobs_script)}  {os.path.join(args.work_dir, 'resubmit_jobs.sh')}") 

#os.popen(f"sed -i 's/#$ -l h_rt=2:00:00/#$ -l h_rt={args.new_runtime}/g' {os.path.join(args.work_dir, 'resubmit_jobs.sh')}") 
#os.popen(f"sed -i -E 's/[#][$] -N.+/#$ -N resubmit/g' {os.path.join(args.work_dir, 'resubmit_jobs.sh')}") 
#os.popen(f"sed -i -E 's/--jobTable=.+/--jobTable={os.path.join(args.work_dir, args.new_jobTable)}/g' {os.path.join(args.work_dir, 'resubmit_jobs.sh')}") 

# regex matches with replacements
mod_rgx = [(r"^#\$ -l h_rt=.+", f"#$ -l h_rt={args.new_runtime}\n"),
           (r"^#\$ -t.+", f"#$ -t 1-{len(failed_runs)}"),
           (r"\t--jobTable=.+", f"\t--jobTable {os.path.join(args.work_dir, args.new_jobTable)}\n"),
           (r"^#\$ -N.+", "#$ -N alphafold_resubmit\n")]

with open(os.path.join(args.work_dir,  args.jobs_script), mode="r", encoding="utf-8") as job_in:
    af_job = job_in.readlines() # read in as a list

    for line in af_job:
        for rgx, repl in mod_rgx:
            if re.match(rgx, line):
                print(f"Modifying line: {line.strip()} --> {repl.strip()}")
                af_job[af_job.index(line)] = repl

print(f"Writing modified job script to {os.path.join(args.work_dir, 'resubmit_jobs.sh')}")
with open(os.path.join(args.work_dir, 'resubmit_jobs.sh'), mode="w", encoding="utf-8") as job_out:
    job_out.writelines(af_job)

print(f"Submitting repeat jobs:\nqsub {os.path.join(args.work_dir, 'resubmit_jobs.sh')}")
os.popen(f"qsub {os.path.join(args.work_dir, 'resubmit_jobs.sh')}")
