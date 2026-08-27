#!/bin/bash

## script to obtain LD among SNPs present in UKB BGEN files
## Maik Pietzner 28/04/2026
## needs to be supplied with 
## ${1} - file with SNPs
## ${2} - chromosome to be used

## export location of files
export dir="<path_to_ukb_bgen>"/
export out="<path_to_file>"

## move to location
cd ${out} 
echo ${1}
echo ${2}

## reassign
snp_list="${1}"
chr="${2}"

# ## activate REGENIE
# source /etc/profile.d/conda.sh
# conda activate regenie_env

##change to X for 23
if [[ $chr -eq 23 ]]; then
  
  export chr="X"
  
fi

## run REGENIE (needs a remove flag for samples to exclude)
# "<path_to_regenie>" \
"<path_to_regenie>" \
--step 2 \
--bgen ${dir}/ukb22828_c${chr}_b0_v3.bgen \
--ref-first \
--extract ${snp_list} \
--minMAC 10 \
--sample ${dir}/ukb22828_c${chr}_b0_v3.sample \
--keep "<path_to_file>" \
--threads 20 \
--compute-corr \
--ignore-pred \
--bsize 10000 \
--out ${snp_list}
