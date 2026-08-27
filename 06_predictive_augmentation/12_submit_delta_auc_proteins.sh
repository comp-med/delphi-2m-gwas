#!/bin/sh

## run protein delta-AUC models

#SBATCH --partition=compute
#SBATCH --job-name=protein_delta_auc
#SBATCH --account=sc-users
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-638
#SBATCH --output="<path_to_file>"%x-%A-%2a.out

## change directory
cd "<path_to_file>"

## one row per (target_pop, target_ehr, outc) -- build once with:
##   cut -f1-3 input/Input.Cox.models.txt | sort -u > input/Input.delta.auc.proteins.txt
## then set --array=1-<number of rows> above.
echo "Job ID: $SLURM_ARRAY_TASK_ID"
target_pop="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $1}' input/Input.delta.auc.proteins.txt)"
target_ehr="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $2}' input/Input.delta.auc.proteins.txt)"
outc="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $3}' input/Input.delta.auc.proteins.txt)"

echo "Node ID: $SLURM_NODELIST"
echo "${target_pop} ${target_ehr} ${outc}"
printf -v date '%(%Y-%m-%d %H:%M:%S)T\n' -1
echo "Starting script at: $date"

R_CONTAINER='"<path_to_container>"'
R_SCRIPT='13_run_delta_auc_cv_proteins.R'
BIND_DIR=""<path_to_file>","<path_to_file>""

singularity exec \
--bind $BIND_DIR \
$R_CONTAINER Rscript $R_SCRIPT "${target_pop}" "${target_ehr}" "${outc}"

printf -v date '%(%Y-%m-%d %H:%M:%S)T\n' -1
echo "Finishing script at: $date"
echo "Done!"