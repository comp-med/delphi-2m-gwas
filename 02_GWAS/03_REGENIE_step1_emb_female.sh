#!/bin/sh

## script to run the first step of REGENIE on the female-only embeddings

#SBATCH --job-name=emb_female

#! select partition
#SBATCH --partition=compute

#! select account
#SBATCH --account=sc-users

#! Specify required run time
#SBATCH --time=48:00:00

#! how many nodes
#SBATCH --nodes=1

#! how much memory per cpu
#SBATCH --mem-per-cpu=5G

#! how many tasks
#SBATCH --ntasks=1

#! how many cpus per task
#SBATCH --cpus-per-task=64

#! no other jobs can be run on the same node
##SBATCH --exclusive

#! What types of email messages do you wish to receive?
#SBATCH --mail-type=FAIL

#! set name
#SBATCH --output=slurm-step1-emb-female.out

export LC_ALL=C \
export LANGUAGE=

## export location of genotype files to be used
export dir="<path_to_file>"
  
## change to relevant directory
cd "<path_to_file>"

#pheno=phenotype.t2d.txt

## Do some logging
echo "Run REGENIE step1 on female embeddings phenotype file"
printf -v date '%(%Y-%m-%d %H:%M:%S)T\n' -1
echo "Starting script at: $date"

## run REGENIE (needs a remove flag for samples to be exclude)
"<path_to_regenie>" \
--step 1 \
--bed ${dir}/ukb22418_allChrs \
--extract "<path_to_file>" \
--keep "<path_to_file>" \
--phenoFile "<path_to_file>" \
--covarFile "<path_to_file>" \
--catCovarList sex,batch \
--threads 64 \
--bsize 1000 \
--lowmem \
--lowmem-prefix gwas/tmpdir/regenie_step1_emb120_female \
--out gwas/input/emb120_female/ukb_emb_female_step1 \
--gz

## Some more logging
printf -v date '%(%Y-%m-%d %H:%M:%S)T\n' -1
echo "Finishing script at: $date"
echo "Done!"
