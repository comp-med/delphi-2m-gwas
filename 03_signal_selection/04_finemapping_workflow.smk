import yaml
import os
import glob
from itertools import product
import json
import sys
import re

# Snakemake workflow to perform fine mapping with SusieR
# Requires regional data to be ready (e.g. 04_process_regions.smk)

try:
    REGIONS_FILE = config["regions_file"]
except KeyError:
    sys.exit("ERROR: You must pass --config regions_file=... to Snakemake")

print(f"Using regions file: {REGIONS_FILE}")

SUBFOLDER = os.path.basename(REGIONS_FILE).split('_')[0]
OUT = f"gwas/output/finemapping/{SUBFOLDER}"

#REGIONS_FILE = "gwas/output/emb120_regions/with_meta/emb103_regions.bed" # The full thing


info = {}
indices = []
with open(REGIONS_FILE,'r') as infile:
    for line in infile:
        idx, phenotype, region, chrom, start, end = line.strip().split()

        info[idx] = {'idx': idx, 'phenotype': phenotype, 'region': region, 'chrom': chrom, 'start': start, 'end': end}
        indices.append(idx)

        # Output directory for every region separately
        os.makedirs(f"{OUT}/{idx}", exist_ok=True)


rule all:
    input:
        expand(f"{OUT}/{{idx}}/finemapped_results.txt.gz", idx = indices)

def get_mem(wildcards, attempt):
    # print(f"This is attempt: {attempt}")
    if attempt == 1:
        return '75G'
    elif attempt == 2:
        return '150G'
    elif attempt == 3:
        return '250G'

rule runFinemapping:
    input:
    output: f"{OUT}/{{idx}}/finemapped_results.txt.gz"
    params:
        index = lambda wcs: info[wcs.idx]['idx'],
        phenotype = lambda wcs: info[wcs.idx]['phenotype'],
        region = lambda wcs: info[wcs.idx]['region'],
        chrom = lambda wcs: info[wcs.idx]['chrom'],
        start = lambda wcs: info[wcs.idx]['start'],
        end = lambda wcs: info[wcs.idx]['end']
    resources:
        mem = get_mem,
        runtime = "5h",
        cpus_per_task = 5
    shell:
        """
            # Variables needed for singularity
            R_CONTAINER='"<path_to_container>"'
            R_SCRIPT='05_run_susie.R'
            BIND_DIR=""<path_to_file>","<path_to_file>","<path_to_file>""
            R_LIBS_USER="<path_to_file>"

            R_LIBS_USER="$R_LIBS_USER" \
            singularity exec \
              --bind $BIND_DIR \
              $R_CONTAINER Rscript $R_SCRIPT {params.index} {params.phenotype} {params.region} {params.chrom} {params.start} {params.end}
        """
