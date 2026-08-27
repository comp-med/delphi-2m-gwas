#!/bin/bash

## script to obtain SNP dosages for protein scores
## Maik Pietzner 05/01/2023

## export location of files
export dir="<path_to_ukb_bgen>"/
export out="<path_to_file>"

## move to location
cd ${out} 
echo ${1}

## loop through all chromosomes
for i in {1..22}
do
echo "Chromosome ${i}"

## create subset bgen file
"<path_to_bgenix>" \
-g ${dir}/ukb22828_c${i}_b0_v3.bgen \
-incl-rsids ../input/snp.list.embeddings.txt > tmp.${i}.bgen

## create dosage file
"<path_to_qctool>" \
-g tmp.${i}.bgen \
-s ${dir}/ukb22828_c${i}_b0_v3.sample  \
-og - \
-ofiletype dosage > tmp.${i}.dosage

## clear up
rm tmp.${i}.bgen

done

## delete all empty files
find -L -maxdepth 1  -type f -size 0 | rm

## concatenate the rest of it (https://stackoverflow.com/questions/16890582/unixmerge-multiple-csv-files-with-same-header-by-keeping-the-header-of-the-firs/16890695#16890695)
awk 'FNR==1 && NR!=1{next;}{print}' tmp.*.dosage > snp.dosage

## get SNP info
cut -f 1-6 -d ' ' snp.dosage > snp.info

## transpose dosage matrix (https://stackoverflow.com/questions/1729824/an-efficient-way-to-transpose-a-file-in-bash)
cut -f 1-6 -d ' ' --complement snp.dosage | awk '
{
  for (i=1; i<=NF; i++)  {
      a[NR,i] = $i
  }
}
NF>p { p = NF }
END {
  for(j=1; j<=p; j++) {
      str=a[1,j]
      for(i=2; i<=NR; i++){
          str=str" "a[i,j];
      }
      print str
  }
}' - > snp.dosage.transpose
  
  
## delete files not longer needed
rm tmp*.dosage
rm snp.dosage
  
  