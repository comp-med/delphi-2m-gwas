#!/bin/sh

## reconstruct exposures from the Delphi embeddings (single multi-core job)

#SBATCH --partition=compute
#SBATCH --job-name=recon_r2
#SBATCH --account=sc-users
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=8G
#SBATCH --output="<path_to_file>"%x-%j.out

## change directory
cd "<path_to_file>"
  
## some logging
echo "Node ID: $SLURM_NODELIST"
echo "CPUs per task: $SLURM_CPUS_PER_TASK"
printf -v date '%(%Y-%m-%d %H:%M:%S)T\n' -1
echo "Starting script at: $date"

## This is the container to be used
R_CONTAINER='"<path_to_container>"'

# This is the script that is executed
# Get with rstudioapi::getSourceEditorContext()$path
R_SCRIPT='09_compute_r2_embeddings.R'

# Enter all directories you need, simply in a comma-separated list
BIND_DIR=""<path_to_file>","<path_to_file>""

## the R script reads $SLURM_CPUS_PER_TASK (passed through by singularity) to size mclapply
singularity exec \
--bind $BIND_DIR \
$R_CONTAINER Rscript $R_SCRIPT

## some more logging
printf -v date '%(%Y-%m-%d %H:%M:%S)T\n' -1
echo "Finishing script at: $date"
echo "Done!"