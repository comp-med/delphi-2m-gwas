#!/bin/sh
#SBATCH --partition=compute
#SBATCH --account=sc-users
#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=20G
#SBATCH --cpus-per-task=20
#SBATCH --array=1-23%3
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-%x-%j.out

cd "<path_to_file>"

echo "Job ID: $SLURM_ARRAY_TASK_ID"

PLINK="<path_to_plink2>"
# PLINK="<path_to_plink1.9>"
DIR="<path_to_ukb_plink2>"/
KEEP="<path_to_file>"

# chromosome handling
if [ "$SLURM_ARRAY_TASK_ID" -eq 23 ]; then
chr="X"
pfile="${DIR}/plink2_chr${chr}"
snpfile="gwas.catalog.snps.23.20260317.txt"
else
chr="$SLURM_ARRAY_TASK_ID"
pfile="${DIR}/plink2_chr${chr}"
snpfile="gwas.catalog.snps.${chr}.20260317.txt"
fi

# 1. find duplicates
$PLINK \
--pfile $pfile \
--keep $KEEP \
--clump $snpfile \
--clump-p1 1 \
--clump-p2 1 \
--clump-r2 0.1 \
--clump-kb 1000000 \
--rm-dup \
--threads 20 \
--out duplicates.${chr}

# 2. detect missing IDs
$PLINK \
--pfile $pfile \
--keep $KEEP \
--clump $snpfile \
--exclude duplicates.${chr}.rmdup.mismatch \
--clump-p1 1 \
--clump-p2 1 \
--clump-r2 0.1 \
--clump-kb 1000000 \
--threads 20 \
--out clumbed.${chr}

# 3. combine exclusions
cat duplicates.${chr}.rmdup.mismatch clumbed.${chr}.clumps.missing_id > exclude.${chr}.txt

# 4. final clumping
$PLINK \
--pfile $pfile \
--keep $KEEP \
--clump $snpfile \
--exclude exclude.${chr}.txt \
--clump-p1 1 \
--clump-p2 1 \
--clump-r2 0.1 \
--clump-kb 1000000 \
--threads 20 \
--out clumbed.${chr}