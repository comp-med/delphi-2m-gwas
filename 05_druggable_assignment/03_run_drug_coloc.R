#!/usr/bin/env Rscript

## script to run coloc for embeddings and drug gwas loci - one REGION per task
## Maik Pietzner 04/06/2026
rm(list = ls())

## get the arguments from the command line
args <- commandArgs(trailingOnly = T)

## little options
options(stringsAsFactors = F)

## set working directory
setwd("<path_to_file>")

## --> packages required <-- ##

require(data.table)
require(coloc)
require(doMC)
require(arrow)
library(polars)
library(fs)
library(checkmate)
# require(smap)

## --> import parameters <-- ##

## one task now handles one genomic REGION (not one region x disease x embedding)
l.num     <- as.numeric(args[1])
# l.num   <- 1

## how many cores for the per-embedding extraction (matches --cpus-per-task)
n.cores   <- 3
registerDoMC(n.cores)

## the cheap embedding signal gate: skip coloc unless the embedding has at least
## suggestive signal in the region (max -log10 p); set to 0 to test everything
emb.log10p.min <- 5

print(l.num)

#-----------------------------------------#
##--   resolve the region for this task --##
#-----------------------------------------#

## full input table (region x disease x embedding) - same file as before
tmp        <- fread("input/input.drug.colocalisation.txt")
## the unique regions, in a stable order; row l.num is this task's region
regions    <- unique(tmp[, .(CHR_numeric, region.start.third.merge, region.end.third.merge)])
setorder(regions, CHR_numeric, region.start.third.merge, region.end.third.merge)

## N.B.: coordinates are in build 38(!)
chr.s      <- regions$CHR_numeric[l.num]
pos.s      <- regions$region.start.third.merge[l.num]
pos.e      <- regions$region.end.third.merge[l.num]

## all drug indications mapped to this region
diseases   <- unique(tmp[ CHR_numeric == chr.s &
                            region.start.third.merge == pos.s &
                            region.end.third.merge == pos.e]$file)
## all embeddings to be tested (constant 120, but read from the file to be safe)
embeddings <- sort(unique(tmp[ CHR_numeric == chr.s &
                                 region.start.third.merge == pos.s &
                                 region.end.third.merge == pos.e]$embedding))

cat("Region", chr.s, pos.s, pos.e, "->", length(diseases), "drug GWAS x",
    length(embeddings), "embeddings\n")
cat("--------------------------------------------------\n")

## functions: region query + the shared per-triplet coloc routine
source("../functions/extract_genomic_region.R")
source("../functions/extract_genomic_data.R")
source("../functions/fn_drug_coloc.R")

#-----------------------------------------#
##-- drug GWAS slices: one per disease  --##
#-----------------------------------------#

## queried ONCE per disease for the whole region (was repeated 120x before)
cat("\n importing", length(diseases), "drug GWAS slices\n")
drug.slices <- lapply(diseases, function(d){
  extract_genomic_data(
    paste0("<path_to_file>",
           d, "/grch38/formatted_sumstats_grch38.parquet"),
    chr = chr.s, bp_min = pos.s, bp_max = pos.e)
})
names(drug.slices) <- diseases

#-----------------------------------------#
##-- liftover for the region (once)     --##
#-----------------------------------------#

cat("\n importing lift over file\n")
lift.coords <- smap::extract_genomic_data(
  parquet_location = "<path_to_file>",
  chr        = paste0("chr", chr.s),
  bp_min     = pos.s,
  bp_max     = pos.e,
  chr_column = "chr_hg38",
  bp_column  = "pos_hg38"
)
## hg19 window that the embedding files (build 37) need to be queried on
low.hg19   <- min(lift.coords$pos_hg19, na.rm = T)
upp.hg19   <- max(lift.coords$pos_hg19, na.rm = T)

#-----------------------------------------#
##-- loop embeddings: extract, gate, run --#
#-----------------------------------------#

## extract every embedding's regional slice in parallel (the slow zcat|awk step,
## now done once per embedding for the whole region rather than once per triplet)
cat("\n extracting", length(embeddings), "embedding slices\n")
emb.slices <- mclapply(embeddings, function(emb){
  fread(cmd = paste0(
    "zcat <path_to_file>",
    emb, ".allchr.results.gz | awk -v chr=", chr.s,
    " -v low=", low.hg19, " -v upp=", upp.hg19,
    " '{if(($1 == chr && $2 >= low && $2 <= upp) || NR == 1) print $0}'"))
}, mc.cores = n.cores)
names(emb.slices) <- embeddings

## test only the embeddings that carry signal in the region (the gate)
emb.test   <- embeddings[ sapply(emb.slices, embedding.has.signal, log10p.min = emb.log10p.min)]
cat("\n", length(emb.test), "of", length(embeddings),
    "embeddings pass the signal gate (-log10p >=", emb.log10p.min, ")\n")
cat("--------------------------------------------------\n")

## run coloc for every surviving embedding x disease pair in the region
for(emb in emb.test){
  for(d in diseases){
    cat("Run coloc with", d, chr.s, pos.s, pos.e, "on", emb, "\n")
    ## copy() so the in-place keys added inside the function do not leak between pairs
    run.drug.coloc(res.gwas    = drug.slices[[d]],
                   res.emb     = copy(emb.slices[[emb]]),
                   lift.coords = lift.coords,
                   embedding   = emb,
                   disease     = d,
                   chr.s       = chr.s,
                   pos.s       = pos.s,
                   pos.e       = pos.e)
  }
}

## --> completion sentinel: reconcile success against the 460 region tasks <-- ##
dir.create("logs/done", showWarnings = F, recursive = T)
fwrite(data.table(region_id   = l.num,
                  chr         = chr.s,
                  pos.s       = pos.s,
                  pos.e       = pos.e,
                  n.diseases  = length(diseases),
                  n.emb.gate  = length(emb.test),
                  n.pairs.run = length(emb.test) * length(diseases)),
       paste0("logs/done/region.", l.num, ".done"),
       sep = "\t", row.names = F, quote = F, na = NA)

cat("\nFinished region", chr.s, pos.s, pos.e, "\n")
cat("--------------------------------------------------\n")