#!/bin/sh

## Script to run the second step of REGENIE on the female-only embeddings
#SBATCH --job-name=emb120_female_s2

#! select partition
#SBATCH --partition=compute

#! select account
#SBATCH --account=sc-users

#! Specify required run time
#SBATCH --time=48:00:00

#! how many nodes
#SBATCH --nodes=1

#! how much memory per cpu
#SBATCH --mem-per-cpu=2G

#! how many tasks
#SBATCH --ntasks=1

#! how many cpus per task
#SBATCH --cpus-per-task=16

#! no other jobs can be run on the same node
##SBATCH --exclusive

#! run as job array
#SBATCH --array=1-23%8

#! What types of email messages do you wish to receive?
#SBATCH --mail-type=FAIL

#! set name
#SBATCH --output=slurm-%x-%j-step2-emb-female.out

source "<path_to_file>"
source "<path_to_file>"

## directory of variant inclusion files
export var=${ukb_imp}/variant_qc/EUR

export input="<path_to_file>"
export current_proj=${ppl}/Wenhuan/project/10_Delphi
export step1_output=${current_proj}/gwas/input/emb120_female
export output=${current_proj}/gwas/output_sex_stratified/emb120_female
export tmpdir=${current_proj}/gwas/tmpdir/emb120_female
export sample_qc="<path_to_file>"

cd ${current_proj}

## Use Array Index to select chromosome
echo "Chromosome: $SLURM_ARRAY_TASK_ID"
export chr=$SLURM_ARRAY_TASK_ID

## Do some logging
pheno=female_norm_embeds_phenotype.txt
echo "Run REGENIE step2 for ${pheno}"
printf -v date '%(%Y-%m-%d %H:%M:%S)T\n' -1
echo "Starting script at: $date"

##change to X for 23
if [[ $chr -eq 23 ]]; then

export chr="X"

## create variant inclusion list for the respective chromosome

awk -v chr=${chr} '{if(NR != 1 && $3 == chr && (($17+0 > 0.001 && $20+0 > 0.8 && $17+0 < 0.999) || ($17+0 > 0.01 && $20+0 > 0.3 && $17+0 < .99))) print $2}' ${var}/ukb_imp_EUR_chr${chr}_snpstat.out > ${tmpdir}/tmp.ex.chr${chr}.list

## run REGENIE
${progs_dir}/regenie/regenie_v4.0.gz_x86_64_ubuntu20_mkl \
--step 2 \
--bgen ${bgen_imp}/ukb22828_c${chr}_b0_v3.bgen \
--ref-first \
--extract ${tmpdir}/tmp.ex.chr${chr}.list \
--sample ${bgen_imp}/ukb22828_cX_b0_v3.sample \
--keep ${sample_qc} \
--phenoFile ${current_proj}/gwas/data/female_norm_embeds_phenotype.txt \
--covarFile ${input}/covariates.txt \
--catCovarList sex,batch \
--threads 16 \
--par-region 'hg19' \
--qt \
--gz \
--pred ${step1_output}/ukb_emb_female_step1_pred.list \
--bsize 400 \
--out ${output}/gwas_emb120_female_chr${chr}

else

## create variant inclusion list for the respective chromosome

awk -v chr=${chr} '{if(NR != 1 && $3 == chr && (($14+0 > 0.001 && $17+0 > 0.8 && $14+0 < 0.999) || ($14+0 > 0.01 && $17+0 > 0.3 &&  $14+0 < 0.99)) && $8+0 > 1e-15) print $2}' ${var}/ukb_imp_EUR_chr${chr}_snpstat.out > ${tmpdir}/tmp.ex.chr${chr}.list

## run REGENIE 
${progs_dir}/regenie/regenie_v4.0.gz_x86_64_ubuntu20_mkl \
--step 2 \
--bgen ${bgen_imp}/ukb22828_c${chr}_b0_v3.bgen \
--ref-first \
--extract ${tmpdir}/tmp.ex.chr${chr}.list \
--sample ${bgen_imp}/ukb22828_c1_b0_v3.sample \
--keep ${sample_qc} \
--phenoFile ${current_proj}/gwas/data/female_norm_embeds_phenotype.txt \
--covarFile ${input}/covariates.txt \
--catCovarList sex,batch \
--threads 16 \
--qt \
--gz \
--pred ${step1_output}/ukb_emb_female_step1_pred.list \
--bsize 400 \
--out ${output}/gwas_emb120_female_chr${chr}

fi

## Some more logging
printf -v date '%(%Y-%m-%d %H:%M:%S)T\n' -1
echo "Finishing script at: $date"
echo "Done!"
