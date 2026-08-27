#!/bin/sh

## run delta-AUC (exposure gain over the Delphi LP)

#SBATCH --partition=compute
#SBATCH --job-name=delta_auc
#SBATCH --account=sc-users
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-1260
#SBATCH --output="<path_to_file>"%x-%A-%2a.out

## NB: set --array=1-XXX to the row count 04 printed for input/Input.delta.auc.txt,
## or override on submit:  sbatch --array=1-$(wc -l < input/Input.delta.auc.txt) 06_submit_delta_auc.sh

## change directory
cd "<path_to_file>"


## Use Array Index to select the task (outcome x population)
echo "Job ID: $SLURM_ARRAY_TASK_ID"
target_pop="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $1}' input/Input.delta.auc.txt)"
target_ehr="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $2}' input/Input.delta.auc.txt)"
outc="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $3}' input/Input.delta.auc.txt)"
token="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $4}' input/Input.delta.auc.txt)"

echo "Node ID: $SLURM_NODELIST"

## Do some logging
echo "${target_pop} ${target_ehr} ${outc} ${token}"
printf -v date '%(%Y-%m-%d %H:%M:%S)T\n' -1
echo "Starting script at: $date"

## This is the container to be used
R_CONTAINER='"<path_to_container>"'

# This is the script that is executed
# Get with rstudioapi::getSourceEditorContext()$path
R_SCRIPT='07_run_delta_auc_cv.R'

# Enter all directories you need, simply in a comma-separated list
BIND_DIR=""<path_to_file>","<path_to_file>""

## The container
singularity exec \
--bind $BIND_DIR \
$R_CONTAINER Rscript $R_SCRIPT "${target_pop}" "${target_ehr}" "${outc}" "${token}"

## Some more logging
printf -v date '%(%Y-%m-%d %H:%M:%S)T\n' -1
echo "Finishing script at: $date"
echo "Done!"
