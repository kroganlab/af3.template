#! /usr/bin/env bash
#$ -S /bin/bash
#$ -cwd

module load CBI
module load r
Rscript installPackages.R
