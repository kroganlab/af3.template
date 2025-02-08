## Krogan Lab Alphafold 3 pipeline compatible with UCSF Wynton HPC 

#### Quickstart: Running Pipeline

1. Clone the github repo to a Wynton working directory
```
git clone https://github.com/kroganlab/af3.template.git 
```

2. Move into `af3.template`. This will be your project working directory
```
cd  af3.template
```

2b. *First time only you run this on Wynton*, make sure you have the R packages you need. This is only necessary for some of the post processing after AlphaFold completes.

```
bash installPackages.R
```

3. Make a new `AlphaFoldJobList.csv` file (match format to `pten.jobTable.txt`)

4. Make a new `masterFasta.fasta` file (match format to `pten_preys.fa`)

5. Edit submission script `af.jobs.sh`, most importantly the number of tasks, but also new file names or job names if desired

6. Submit job with `qsub af.jobs.sh`

7. View queue and job status with `qstat`

#### Postprocessing 

1. Once the job is finised, check your output folder. Each PPI should have its own subdirectory. Within the output directory, the key output files are: `{PPI_name}_summaryScores.csv`, `{PPI_name}.msa.png` and `{PPI_name}.pae.png`

2. Sometimes AF jobs fail to complete, which can be often due to timeout. To check for incomplete runs, you can use the following bash one-liner. Just `cd` into the output directory and run:
```
for dir in *;do if [[ ! -e ./$dir/${dir}_model.cif ]]; then echo $dir;fi;  done
```
This checks for the existance of the top scoring model, and prints the name of the directory if this isn't found (you can remove the `!` symbol to print names of completed runs)

3. You can capture these incomplete runs and convert to a new `AlphaFoldJobList_remainingRuns.csv` for a fresh submission, using the following one-liner:
 (assumption here is you are searching from the parent directory of `outputDir` i.e `./outputDir/protein1__protein2/outputFiles`)
```
find ./outputDir -maxdepth 1 -type d '!' -exec test -e "{}/ranking_scores.csv" ';' -print | cut -d '/' -f3 | awk '{FS="__";OFS=","}  {print(NR-1,toupper($1),toupper($2))}' > ./AlphaFoldJobList_remainingRuns.csv
```
Be sure to edit the `af.jobs.sh` script and extend the job runtime beyond two hours to avoid another timeout!
