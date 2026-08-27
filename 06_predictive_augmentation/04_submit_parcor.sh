#!/bin/sh

## run partial correlation network among lifetime ICD-10 codes

#SBATCH --partition=compute
#SBATCH --job-name=parcor_icd10
#SBATCH --account=sc-users
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-3
#SBATCH --output="<path_to_file>"%x-%A-%2a.out

## change directory
cd "<path_to_file>"

## let Rfast / BLAS use the cores Slurm granted
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

## Use Array Index to select the run
echo "Job ID: $SLURM_ARRAY_TASK_ID"
## column 1 = sex; columns 2/3 = optional test sizes (empty on a real run)
sex="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $1}' input/Input.partial.correlation.txt)"
nrow="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $2}' input/Input.partial.correlation.txt)"
ncol="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $3}' input/Input.partial.correlation.txt)"

echo "Node ID: $SLURM_NODELIST"

## Do some logging
echo "sex=${sex} nrow=${nrow} ncol=${ncol}"
printf -v date '%(%Y-%m-%d %H:%M:%S)T\n' -1
echo "Starting script at: $date"

## This is the container to be used
R_CONTAINER='"<path_to_container>"'

# This is the script that is executed
# Get with rstudioapi::getSourceEditorContext()$path
R_SCRIPT='05_compute_partial_correlation_disease.R'

# Enter all directories you need, simply in a comma-separated list
BIND_DIR=""<path_to_file>","<path_to_file>""

## The container
singularity exec \
--bind $BIND_DIR \
$R_CONTAINER Rscript $R_SCRIPT "${sex}" "${nrow}" "${ncol}"

## Some more logging
printf -v date '%(%Y-%m-%d %H:%M:%S)T\n' -1
echo "Finishing script at: $date"
echo "Done!"
