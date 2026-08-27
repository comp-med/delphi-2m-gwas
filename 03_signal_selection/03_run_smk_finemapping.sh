# Get the conda env that contains snakemake
source "<path_to_conda>"/etc/profile.d/conda.sh
conda activate snakemake

cd "<path_to_file>"

echo 'Changing to: '
echo $PWD

# To get the accessories from the Snakemake config
export PATH="<path_to_snakemake_profile>":$PATH

echo 'Snakemake version: '
snakemake -v

SNAKEFILE="04_finemapping_workflow.smk"
REGIONS_DIR="gwas/output/emb120_regions/with_meta"

# Unlock once, just in case
# Get one example file to satisfy --config during --unlock
first_file=$(ls gwas/output/emb120_regions/with_meta/emb*_regions.bed | head -n1)
first_emb=$(basename "$first_file")
first_emb=${first_emb%%_*}

# Unlock once with config
snakemake \
  --unlock \
  --executor slurm \
  --snakefile "$SNAKEFILE" \
  --profile official_slurm_github \
  --config regions_file="$first_file" embedding="$first_emb" || true

for regions_file in "$REGIONS_DIR"/emb*_regions.bed; do
  [[ -e "$regions_file" ]] || continue

  name=$(basename "$regions_file")       # e.g. emb37_regions.bed
  emb=${name%%_*}                        # -> emb37

  echo "==============================="
  echo "Running fine-mapping for $emb"
  echo "Regions file: $regions_file"
  echo "==============================="

  snakemake \
    --executor slurm \
    --snakefile "$SNAKEFILE" \
    --jobs 100 \
    --keep-going \
    --profile official_slurm_github \
    --use-conda \
    --rerun-triggers mtime \
    --config regions_file="$regions_file" embedding="$emb" \

done
