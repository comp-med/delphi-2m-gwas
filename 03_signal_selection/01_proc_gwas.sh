#!/usr/bin/env bash

source "<path_to_conda>"/etc/profile.d/conda.sh

conda activate gwas

cd "<path_to_file>"


# 1. group files according to the embedings
echo "Start grouping regenie results according to embedding dimension"
IN_DIR="output/emb120"
OUT_DIR="output/emb120_group"
mkdir -p "$OUT_DIR"

for f in "$IN_DIR"/*.regenie.gz; do
  tmp_part=${f#"$IN_DIR"/gwas_emb120_}
  part=${tmp_part%.regenie.gz}
  chrom=${part%%_*}
  emb=${part##*_}

  dest_dir="$OUT_DIR/$emb"
  mkdir -p "$dest_dir"

  echo "grouping $emb $chrom ..."

  cp -v "$f" "$dest_dir/"

done




# 2. intersect own result with UKB WGS results (all 763 studies) according to rsid (could also have done at first)

# --- INPUTS (edit these) ---
MY_DIR="output/emb120_group"            # folder with the chr*.regenie.gz
STUDY_DIR="<path_to_file>"  # either a single file OR a directory
OUT_DIR="output/emb120_intersect"       # where to write filtered chr files
mkdir -p "$OUT_DIR"

TMPD="output/misc"


# Space-delimited (.regenie.gz)
extract_ids_space() {
  local infile="$1" header_name="$2" out="$3"
  zcat "$infile" \
  | awk -F' +' -v want="$header_name" -v fname="$infile" '
      NR==1{
        for(i=1;i<=NF;i++){ s=$i; sub(/\r$/,"",s); if(tolower(s)==tolower(want)) c=i }
        if(!c){ print "WARN: no column \"" want "\" in " fname; exit 0 }
        next
      }
      { print $c }
    ' >> "$out"
}

# Tab-delimited (.tsv.gz)
extract_ids_tsv() {
  local infile="$1" header_name="$2" out="$3"
  zcat "$infile" \
  | awk -F'\t' -v want="$header_name" -v fname="$infile" '
      NR==1{
        for(i=1;i<=NF;i++){ s=$i; sub(/\r$/,"",s); if(tolower(s)==tolower(want)) c=i }
        if(!c){ print "WARN: no column \"" want "\" in " fname; exit 0 }
        next
      }
      { print $c }
    ' >> "$out"
}


# 1) Collect RSIDs for the files using column "ID"
: > "$TMPD/emb120.ids"
shopt -s nullglob
# loop through subfolders (emb1 ... emb120)
for d in "$MY_DIR"/*/; do
  [[ -d "$d" ]] || continue

  # all .regenie.gz files in this subfolder
  files=( "$d"/*.regenie.gz )
  if (( ${#files[@]} == 0 )); then
    echo "Skip $(basename "$d"): no .regenie.gz files."
    continue
  fi

  echo "Processing $(basename "$d") with ${#files[@]} files..."
  for f1 in "${files[@]}"; do
    extract_ids_space "$f1" "ID" "$TMPD/emb120.ids"
  done
done

# unique+sorted (write back to same file or change output name if you prefer)
LC_ALL=C sort -u "$TMPD/emb120.ids" -o "$TMPD/emb120.ids"
echo "Wrote $(wc -l < "$TMPD/emb120.ids") unique IDs to $TMPD/emb120.ids"



# 2) Collect RSIDs for STUDY using column "rsid" (file)
: > "$TMPD/study_all.ids"

# collect all study files (adjust patterns if needed)
files=( "$STUDY_DIR"/*.tsv.gz )
if (( ${#files[@]} == 0 )); then
  echo "No study files found under: $STUDY_DIR"
else
  echo "Found ${#files[@]} study files."
  for f in "${files[@]}"; do
    # study files are TSV (tabs); rsid column name is "rsid"
    extract_ids_tsv "$f" "rsid" "$TMPD/study_all.ids"
  done
fi
LC_ALL=C sort -u "$TMPD/study_all.ids" -o "$TMPD/study_all.ids"


# 3) Intersection of RSIDs ---
LC_ALL=C comm -12 "$TMPD/emb120.ids" "$TMPD/study_all.ids" > "$TMPD/emb120_intersect.ids"
echo "Intersect RSIDs: $(wc -l < "$TMPD/emb120_intersect.ids")"



# 4) Filter EACH of the chromosome files by the intersected RSIDs ---
# Recursively process emb subfolders: emb120_group/emb*/*.regenie.gz
shopt -s nullglob

for d in "$MY_DIR"/emb*/; do
  [[ -d "$d" ]] || continue

  # where to write (mirror tree under OUT_DIR)
  emb_name=$(basename "$d")                 # e.g. emb7
  dest_dir="$OUT_DIR/$emb_name"
  mkdir -p "$dest_dir"

  files=( "$d"/*.regenie.gz )
  if (( ${#files[@]} == 0 )); then
    echo "Skip $emb_name: no *.regenie.gz"
    continue
  fi

  for f in "${files[@]}"; do
    base=$(basename "$f" .regenie.gz)       # e.g. chr1_emb7.log10p....sorted
    out="$dest_dir/${base}.intersect.regenie.gz"
    tmp="$dest_dir/tmp.${base}.intersect.regenie.gz"

    echo "Filtering $emb_name: $(basename "$f") -> $(basename "$out")"

    # Filter by intersect IDs, keep SPACE delimiter
    zcat "$f" \
    | awk -F' +' -v OFS=' ' '
        NR==FNR { keep[$1]=1; next }              # read intersect IDs
        FNR==1 {
          for(i=1;i<=NF;i++) if(tolower($i)=="id") idc=i
          print $0; next                           # header
        }
        (idc && ($idc in keep))                    # keep matching rows
      ' "$TMPD/emb120_intersect.ids" - \
    | bgzip > "$tmp"

    # count non-empty lines; require >1 to keep (header + at least 1 data row)
    n=$(zcat -f -- "$tmp" 2>/dev/null | awk 'NF>0{c++} END{print c+0}')
    if [ "${n:-0}" -le 1 ]; then
      echo "Removing (no data rows): $tmp"
      rm -f -- "$tmp"
      continue
    fi    

    # Sort by LOG10P numeric DESC (space in/out)
    csvtk sort -d ' ' -D ' ' -k LOG10P:nr "$tmp" | bgzip > "$out"
    rm -f "$tmp"

    # Row counts (header-aware)
    prev_total=$(zcat "$f"  | wc -l)
    new_total=$(zcat "$out" | wc -l)
    prev_data=$(( prev_total>0 ? prev_total-1 : 0 ))
    new_data=$(( new_total>0 ? new_total-1 : 0 ))
    printf "Rows (incl header): prev=%d, new=%d | data rows: prev=%d, new=%d | kept=%.2f%%\n" \
      "$prev_total" "$new_total" "$prev_data" "$new_data" \
      "$(awk -v a="$prev_data" -v b="$new_data" 'BEGIN{print (a?100*b/a:0)}')"
  done
done


# 3. filter and sort
LOG10P="9.3802"
IN_DIR="output/emb120_intersect"
OUT_DIR="output/emb120_intersect_proc"
mkdir -p "$OUT_DIR"

for d in "$IN_DIR"/emb*/; do
  [[ -d "$d" ]] || continue

  emb_name=$(basename "$d")                 
  dest_dir="$OUT_DIR/$emb_name"
  mkdir -p "$dest_dir"

  # all .regenie.gz files in this subfolder
  files=( "$d"/*.regenie.gz )
  if (( ${#files[@]} == 0 )); then
    echo "Skip $emb_name: no .regenie.gz files."
    continue
  fi

  echo "Filtering $emb_name (LOG10P >= $LOG10P)..."
  for f in "${files[@]}"; do
    base=$(basename "$f" .regenie.gz)       
    tmp="$dest_dir/tmp.${base}.log10p${LOG10P}.regenie.gz"
    
    # Read space-delimited, write space-delimited
    csvtk filter2 -d ' ' -D ' ' -f "\$LOG10P >= ${LOG10P}" "$f" | bgzip > "$tmp"
    # If only header remains (<=1 line), skip
    if [ "$(wc -l < "$tmp")" -le 1 ]; then
      echo " $base No rows passed filter (only header). Skipping."
      rm -f "$tmp"
      continue
    fi

    out="$dest_dir/${base}.log10p${LOG10P}.regenie.gz"

    echo "  Sorting and writing -> $out"
    csvtk sort -d ' ' -D ' ' -k LOG10P:nr "$tmp" | bgzip > "$out"
    rm -f "$tmp"
  done
done


# 4. merged files across all embs and chromosome to one file
IN_DIR="output/emb120_intersect_proc"     # parent containing emb1..emb120 subfolders
OUT_DIR="output/merged_emb120_intersect_proc"   # where to write per-emb merges
mkdir -p "$OUT_DIR"

for emb in {1..120}; do
  emb_name="emb${emb}"
  emb_dir="$IN_DIR/$emb_name"
  [[ -d "$emb_dir" ]] || { echo "Skip $emb_name: no folder"; continue; }

  # any chr files?
  any=( "$emb_dir"/*.regenie.gz )
  if (( ${#any[@]} == 0 )); then
    echo "Skip $emb_name: no regenie.gz files"
    continue
  fi

  out_plain="$OUT_DIR/${emb_name}.merged.regenie"
  : > "$out_plain"
  header_written=0

  echo "Merging $emb_name in chr order (1..22)..."
  for chr in {1..22}; do
    # Support e.g. chr7_embX... or chr7... (both patterns)
    files=( "$emb_dir"/gwas_emb120_chr${chr}_${emb_name}*.regenie.gz)
    for f in "${files[@]}"; do
      [ -e "$f" ] || continue
      total=$(zcat "$f" | wc -l)
      (( total > 0 )) || continue

      if (( header_written == 0 )); then
        # header with 'emb' inserted as first column
        zcat "$f" | head -1 | awk 'BEGIN{OFS=" "} {print "emb", $0}' > "$out_plain"
        header_written=1
      fi

      # append data rows, prefixing with embX
      if (( total > 1 )); then
        zcat "$f" | tail -n +2 | awk -v e="$emb_name" 'BEGIN{OFS=" "} {print e, $0}' >> "$out_plain"
      fi
    done
  done

  if (( header_written == 1 )); then
    bgzip -f "$out_plain"
    echo "  → $OUT_DIR/${emb_name}.merged.regenie"
  else
    rm -rf "$out_plain"
    echo "  No non-empty chr files for $emb_name; nothing written."
  fi
done



# 5. obtain regional lead signals
IN_DIR="output/merged_emb120_intersect_proc"
OUT_DIR="output/emb120_regions"
PEAK_DIR="output/emb120_sentinels"
mkdir -p "$OUT_DIR"
mkdir -p "$PEAK_DIR"

## get all signals of interest (treat extended MHC region as one): normal gw-threshold
LOG10P="9.3802"
WINDOW=500000
MHC_CHR=6
MHC_START=25500000
MHC_END=34000000

BEDTOOLS="<path_to_bedtools>"  # or just "bedtools" if on PATH


for f in "$IN_DIR"/*.regenie.gz; do
  [[ -e "$f" ]] || continue
  emb_name=$(basename "$f")
  emb_name=${emb_name%%.*}                       # e.g. emb7 or emb7.merged_with_emb -> emb7
  out="$OUT_DIR/${emb_name}_regions.bed"

  echo "Building regions for $emb_name from $(basename "$f") ..."

  # Stream file, filter LOG10P > threshold, convert to windows with MHC special rules,
  # sort and merge.
  zcat "$f" \
  | awk -F' +' -v OFS='\t' -v T="$LOG10P" -v W="$WINDOW" -v M1="$MHC_START" -v M2="$MHC_END" '
      # Skip header (first line or any line whose first field is CHROM)
      NR==1 || $1=="CHROM" { next }

      {
        chrom = $2
        pos   = $3 + 0
        logp  = ($14 + 0)        # LOG10P assumed in column 14 (emb is the 1st column)
        if (logp <= T) next

        # normalize numeric chromosome for comparisons; keep original label for output
        cn = chrom; sub(/^chr/,"",cn)

        # default ±W window
        start = pos - W; if (start < 0) start = 0
        end   = pos + W

        # extended MHC handling on chr6
        if (cn + 0 == 6) {
          if (pos >= M1+W && pos <= M2-W) {
            start = M1; end = M2
          } else if (pos >= M1 && pos <= M1+W) {
            start = (pos-W < 0 ? 0 : pos-W); end = M2
          } else if (pos >= M2-W && pos <= M2) {
            start = M1; end = pos+W
          }
        }

        print chrom, start, end
      }
    ' \
  | sort -k1,1 -k2,2n \
  | "$BEDTOOLS" merge -i - > "$out"

  echo "  -> wrote $out"
done

echo "Done: regions per embedding in $OUT_DIR"


## now add the strongest signal for each region

for f in "$IN_DIR"/*.regenie.gz; do
  emb=$(basename "$f")
  emb=${emb%%.*}                                   # e.g., emb7
  regions="$OUT_DIR/${emb}_regions.bed"
  [[ -s "$regions" ]] || { echo "Skip $emb: no regions"; continue; }

  out="$PEAK_DIR/${emb}_regional_sentinels.txt"
  echo "Picking strongest signals for $emb ..."
  : > "$out"

  # For each region: pick the row with max |BETA/SE|
  # Notes: CHROM=$2, POS=$3, LOG10P=$14, BETA=$11, SE=$12
  while read -r reg_chr reg_start reg_end; do
    zcat "$f" \
    | awk -F' +' -v OFS=' ' \
        -v T="$LOG10P" \
        -v rchr="$reg_chr" -v rs="$reg_start" -v re="$reg_end" '
        BEGIN{
          # normalize region chr to numeric if possible (strip leading "chr")
          rc=rchr; sub(/^chr/,"",rc)
          max=-1; best=""
        }
        NR==1 { next }                        # skip header
        {
          chrom=$2; pos=$3+0; logp=$14+0; b=$11+0; se=$12+0
          if (logp <= T) next

          # normalize data chrom (strip leading "chr")
          cn=chrom; sub(/^chr/,"",cn)

          # must be within this region
          if (cn==rc && pos>=rs && pos<=re) {
            if (se!=0) {
              z = (b/se); if (z<0) z=-z
              if (z>max) { max=z; best=$0 }
            }
          }
        }
        END{
          if (best!="") print best, rs, re    # append region bounds to the line
        }' >> "$out"
  done < "$regions"

  echo "  -> wrote $out"
done


# 6. edit region file, 3 columns -> 6 columns for finemapping (exclude MHC region)
#!/usr/bin/env bash

REGIONS_DIR="output/emb120_regions"   # adjust if needed
OUT_DIR="${REGIONS_DIR}/with_meta"
mkdir -p "$OUT_DIR"
echo "make dir"

# MHC definition
MHC_CHR=6
MHC_START=25500000
MHC_END=34000000

for f in "$REGIONS_DIR"/*_regions.bed; do
  [[ -e "$f" ]] || continue

  base=$(basename "$f")
  phenotype=${base%_regions.bed}
  out="${OUT_DIR}/${phenotype}_regions.bed"
  tmp="${out}.tmp"

  awk -v PHENO="$phenotype" \
      -v OFS='\t' \
      -v MHC_CHR="$MHC_CHR" \
      -v MHC_START="$MHC_START" \
      -v MHC_END="$MHC_END" '
    # skip blank lines or comments
    /^\s*$/ || /^#/ { next }

    {
      chr=$1; start=$2; end=$3

      # normalize "chr6" -> "6"
      sub(/^chr/, "", chr)

      # ---- skip regions that overlap MHC ----
      # overlap if not (end < MHC_START or start > MHC_END)
      if (chr + 0 == MHC_CHR && !(end < MHC_START || start > MHC_END))
        next

      # increment global index only for kept regions
      idx++

      # per-chromosome counter starting at 0 (only for kept regions)
      c = count[chr] + 0
      region_id = chr "_" c
      count[chr] = c + 1

      print idx, PHENO, region_id, chr, start, end
    }
  ' "$f" > "$tmp"

  if [[ -s "$tmp" ]]; then
    mv "$tmp" "$out"
    echo "Wrote: $out"
  else
    rm -f "$tmp"
    echo "Skipped ${phenotype}: only MHC regions (no non-MHC regions left)."
  fi
done



