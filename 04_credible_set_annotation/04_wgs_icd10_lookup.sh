#!/bin/sh
#SBATCH --partition=compute
#SBATCH --account=sc-users
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=2
#SBATCH --array=1-763
#SBATCH --mail-type=FAIL
#SBATCH --output="<path_to_file>"%x-%j.out

## change to the relevant director
cd "<path_to_file>"

## Use Array Index to select features
echo "Job ID: $SLURM_ARRAY_TASK_ID"
icd10="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $1}' ../input/icd10.gwas.stats)"

## do look-up: may well contain SNPs not needed
zgrep -f ../input/Lookup.SNPs.WGS.stats.20260421.txt "<path_to_file>" > lookup.${icd10}.tsv

