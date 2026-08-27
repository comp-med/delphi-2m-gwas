#!/bin/sh

## script to run drug colocalisation - one task per genomic region

#SBATCH --partition=compute
#SBATCH --job-name=coloc_drugs
#SBATCH --account=sc-users
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=10G
#SBATCH --array=421,435,56,152,421%50
#SBATCH --output="<path_to_file>"%x-%A-%2a.out

## change directory
cd "<path_to_file>"

## Use Array Index to select features
echo "Job ID: $SLURM_ARRAY_TASK_ID"

## Do some logging
printf -v date '%(%Y-%m-%d %H:%M:%S)T\n' -1
echo "Starting script at: $date"

## This is the container to be used
R_CONTAINER='"<path_to_container>"'

# This is the script that is executed
# Get with rstudioapi::getSourceEditorContext()$path
R_SCRIPT='03_run_drug_coloc.R'

# Enter all directories you need, simply in a comma-separated list
BIND_DIR=""<path_to_file>","<path_to_file>""

## The container - the array index is now a REGION index (1..460)
singularity exec \
--bind $BIND_DIR \
$R_CONTAINER Rscript $R_SCRIPT $SLURM_ARRAY_TASK_ID

## Some more logging
printf -v date '%(%Y-%m-%d %H:%M:%S)T\n' -1
echo "Finishing script at: $date"
echo "Done!"