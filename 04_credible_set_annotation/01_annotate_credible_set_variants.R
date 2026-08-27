#######################################################
#### Process credible set variants from embeddings ####
#### Maik Pietzner                      26/02/2026 ####
#######################################################

rm(list=ls())
setwd("<path_to_file>")
options(stringsAsFactors = F)
load(".RData")

## --> packages needed <-- ##
require(data.table)
require(doMC)
require(igraph)

###########################################
####          import results           ####
###########################################

## import lead credible set variants
res.credible <- fread("<path_to_file>")

###########################################
####            LD-clumping            ####
###########################################

## export list of SNPs
write.table(unique(res.credible$id), "snp.list.embeddings.txt", col.names = F, row.names = F, quote = F)

## run
system("./obtain_snps.sh")

## import
snp.dat        <- fread("../tmpdir/snp.dosage.transpose")
snp.info       <- fread("../tmpdir/snp.info")
## assign names
names(snp.dat) <- c("f.eid", snp.info$rsid)

## subset to unrelated EUR population
tmp.dat        <- fread("<path_to_file>")
## apply
snp.dat        <- snp.dat[ f.eid %in% tmp.dat$V1]

## compute LD matrix
snp.ld         <- Rfast::cora(snp.dat[,-1])
## convert to data table (keep only entries from upper triangle)
snp.ld         <- melt(
  as.data.table(snp.ld, keep.rownames="var1"),
  id.vars="var1"
)[, `:=`(var1=as.character(var1), var2=as.character(variable))
][var1 <= var2]
## convert to r2
snp.ld[, value := value^2]

## convert to graph (use max to allow for edge weights)
ld.sub        <- graph_from_data_frame(snp.ld[ value >= .7])
## get all separate components
ld.sub        <- components(ld.sub)$membership
## convert to data frame
ld.sub        <- data.table(ID=names(ld.sub), R2.group=ld.sub)

## add to credible set results
res.credible  <- merge(res.credible, ld.sub, by.x = "id", by.y = "ID")

###########################################
####        create complete set        ####
###########################################

## --> call in regional sentinels <-- ##

## import regional sentinel results to include findings from the MHC region
res.regional        <- fread("<path_to_file>")
## adopt names
jj                  <- fread("zcat <path_to_file> | head -2")
## assign names
names(res.regional) <- c("emb", names(jj), "region.start", "region.end")

## --> separate MHC region file <-- ##

## create MHC regions file to the credible set variant file
tmp.mhc             <- res.regional[ CHROM == 6 & GENPOS >= 25500000+500000 & GENPOS <= 34000000-500000]
names(tmp.mhc)      <- gsub("ID", "id", names(tmp.mhc))

## --> pull stats from REGENIE files for credible set variants <-- ##

## write file with all relevant IDs (credible sets and MHC)
write.table(unique(c("ID", res.credible$id, res.regional$ID)), "snp.list.embeddings.complete.txt", col.names = F, row.names = F, quote = F)

## iterate through all embeddings
snp.lookup.emb      <- rbindlist(mclapply(1:120, function(x){
  ## grep the relevant results
  tmp <- fread(cmd=paste0("zgrep -wf snp.list.embeddings.complete.txt <path_to_file>", x,".allchr.results.gz"))
  ## add embedding and return
  return(data.table(emb = paste0("emb", x), tmp))
}, mc.cores = 10))

## --> create harmonized results set <-- ##

## combine
res.complete        <- merge(res.credible, snp.lookup.emb, by.x = c("id", "emb"), by.y = c("ID", "emb"))
## add MHC region
res.complete        <- rbind(res.complete, tmp.mhc, fill = T)

## delete some columns
res.complete[, c("chrom", "genpos", "allele0", "allele1", "pval_marginal", "beta_marginal", "se_marginal", "R2_leadvariant", "pheno") := NULL]
## harmonize columns
res.complete[, R2.group := ifelse(!is.na(R2.group), R2.group, 0)]
res.complete[, startpos_region := ifelse(!is.na(startpos_region), startpos_region, region.start)]
res.complete[, endpos_region := ifelse(!is.na(endpos_region), endpos_region, region.end)]
## edit some columns
res.complete[, c("region.start", "region.end") := NULL]

## export for gene assignment
write.table(res.complete, "Lead.credible.set.variants.incl.MHC.Delphi.embeddings.UKB.20260415.txt", sep = "\t", row.names = F)

###########################################
####          gene assignment          ####
###########################################

#----------------------------------------#
##-- import human genes for reference --##
#----------------------------------------#

## import gene list: build37 coordinates from here https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/
# human.genes     <- rtracklayer::readGFF("<path_to_file>",
#                                         filter = list(type=c("gene")))
## https://ftp.ensembl.org/pub/grch37/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz
human.genes     <- rtracklayer::readGFF("Homo_sapiens.GRCh37.87.gtf.gz",
                                        filter = list(type=c("gene")))
human.genes     <- as.data.table(human.genes)
## parse some chromosome names: careful chromosome is stored as factor!!
human.genes[, CHROM := ifelse(seqid == "X", 23, as.numeric(as.character(seqid)))] ## assigns some chromosomes wrongly!!
human.genes     <- human.genes[ !is.na(gene_name)]
## drop some unspecific genes
sort(table(human.genes$gene_biotype))
human.genes     <- human.genes[ gene_biotype %in% c("protein_coding", "processed_transcript", "IG_J_gene", "IG_C_gene", "IG_V_gene", "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene")]

#----------------------------------------#
##--       import V2G assignments     --##
#----------------------------------------#

## N.B.: candidate effector genes were assigned with the ensemble variant-to-gene
## classifier developed and described separately; that code is not part of this
## repository. The table imported here carries one row per locus with the columns
## consumed downstream: 'locus_id', 'ensembl_gene_id', 'hgnc_symbol',
## 'maximum.prob', 'median.prob', 'index.max.prob', 'final.variant.gene.rank' and
## 'distance_to_footprint'. Genes were retained where the maximum classifier
## probability exceeded 0.25, otherwise the nearest protein-coding gene was used.
res.v2g          <- fread("<path_to_file>")

## --> annotate closet gene <-- ##

## import function to do so
source("../functions/map_closest_gene.R")
## apply
tmp              <- map_nearest_gene(unique(res.complete[ , .(id, CHROM, GENPOS)]), human.genes)

#----------------------------------------#
##--           add information        --##
#----------------------------------------#

## closest gene
res.complete     <- merge(res.complete, tmp, by = c("id", "CHROM", "GENPOS"))
## V2G assignment
res.complete     <- merge(res.complete, res.v2g, by.x = "id", by.y = "locus_id", all.x = T)

## how often not the nearest gene
res.complete[ !is.na(hgnc_symbol) & hgnc_symbol != nearest_gene_name]

## store
write.table(res.complete, "Lead.credible.set.variants.incl.MHC.Delphi.embeddings.UKB.GWAS.plus.catalog.gene.annotation.20260421.txt", sep = "\t", row.names = F)

###########################################
####       GWAS catalog overlap        ####
###########################################

#---------------------------------#
##--      obtain proxies       --##
#---------------------------------#

## export to generate a list of proxy variables
for(j in 1:23){
  ## write all rsIDs into file
  write.table(unique(res.complete[ CHROM == j]$id), paste("snp.list", j, "txt", sep="."), col.names = F, row.names = F, quote = F)
}
## needs to run script '02_extract_proxies.sh' as batch command in the cluster

## import: N.B.: SNP pairs are only listed once!
r2.proxies <- rbindlist(lapply(1:23, function(x) fread(paste0("ld.proxies.", x,".vcor"))))
## edit naming
names(r2.proxies) <- gsub("#", "", names(r2.proxies))

#-------------------------------------#
##--      lift over to build 38    --##
#-------------------------------------#

## import file
lift.over  <- fread("<path_to_file>",
                    select = c("chr_hg19", "pos_hg19", "pos_hg38"))
lift.over  <- unique(lift.over[, .(chr_hg19, pos_hg19, pos_hg38)])
## align chromosome naming
lift.over[, chr_hg19 := gsub("chr", "", chr_hg19)]

## align coding
r2.proxies[, CHROM_A := as.character(CHROM_A)]
r2.proxies[, CHROM_B := as.character(CHROM_B)]
## merge with proxies
r2.proxies <- merge(r2.proxies, lift.over, by.x = c("CHROM_A", "POS_A"), by.y = c("chr_hg19", "pos_hg19"), all.x=T) 
r2.proxies <- merge(r2.proxies, lift.over, by.x = c("CHROM_B", "POS_B"), by.y = c("chr_hg19", "pos_hg19"), all.x=T, suffixes = c("_A", "_B")) 

## drop
rm(lift.over); gc(reset = T)

#-------------------------------------#
##--       import GWAS catalog     --##
#-------------------------------------#

## import latest release (26/02/2026)
gwas.catalogue           <- fread("gwas_catalog_download_20260317.tsv")

## rename some columns
names(gwas.catalogue)[8] <- "TRAIT"

## prune GWAS catalogue data
gwas.catalogue           <- gwas.catalogue[ !is.na(`OR or BETA`) & is.finite(`OR or BETA`) & CHR_ID != ""]
## generate risk allele and drop everything w/o this information
gwas.catalogue[, riskA := sapply(`STRONGEST SNP-RISK ALLELE`, function(x) strsplit(x,"-")[[1]][2])] 
gwas.catalogue[, riskA := trimws(riskA, which = "b")] 
## drop interaction entries
ii                       <- grep("[0-9]", gwas.catalogue$riskA)
gwas.catalogue           <- gwas.catalogue[-ii,]
## only genome-wide significant ones
gwas.catalogue           <- gwas.catalogue[ PVALUE_MLOG > 7.3 ]
## N = 625,252entries

## create another entry to possible merge on (careful, genome build 38 mapping)
gwas.catalogue[, snp.id := paste0(ifelse(CHR_ID == "X", 23, CHR_ID), ":", CHR_POS)]

#-------------------------------------#
##-- map variants to the GWAS cat. --##
#-------------------------------------#

## get all variants that need to be annotated
emb.var  <- unique(res.complete[, .(id, CHROM, GENPOS, ALLELE0, ALLELE1)])
## n = 318 variants

## import function to do so
source("../functions/annotate_gwas_catalog.R")

## lenient threshold
gwas.r1      <- annotate.variants.gwas(emb.var, gwas.catalogue, r2.proxies, .1)

## more stringent threshold
gwas.r8      <- annotate.variants.gwas(emb.var, gwas.catalogue, r2.proxies, .8)

## combine both
tmp.gwas     <- merge(gwas.r1, gwas.r8, by = c("id", "CHROM", "GENPOS", "ALLELE0", "ALLELE1"), suffixes = c(".r1", ".r8"))

## add to results data set
res.complete <- merge(res.complete, tmp.gwas)

## store
write.table(res.complete, "Lead.credible.set.variants.incl.MHC.Delphi.embeddings.UKB.GWAS.plus.catalog.gene.annotation.20260421.txt", sep = "\t", row.names = F)

###########################################
####        EFO term enrichment        ####
###########################################

## create position
gwas.catalogue[, GENPOS.numeric.hg38 := as.numeric(CHR_POS)]

#--------------------------------#
##--   LD-clump GWAS catalog  --##
#--------------------------------#

## export all snps
gwas.snps         <- unique(gwas.catalogue[, .(CHR_ID, SNPS, PVALUE_MLOG)])
## make chromosome numeric
gwas.snps[, CHR_ID := as.numeric(ifelse(CHR_ID == "X", 23, CHR_ID))]
## to enable clumping based on pseudo p-values
gwas.snps[, P := 10^(-PVALUE_MLOG/100)]
gwas.snps[, SNP := SNPS]
## prune redundant SNPs
gwas.snps         <- gwas.snps[ order(SNP, P)]
gwas.snps[, ind := 1:.N, by = "SNP"]
gwas.snps <- gwas.snps[ ind == 1]

## write to file for clumping
for(j in 1:23){
  ## store by chromosome
  write.table(gwas.snps[ CHR_ID == j, .(SNP, P)], paste0("../tmpdir/gwas.catalog.snps.", j, ".20260317.txt"), sep = "\t", row.names = F, quote = F)
}

## --> import results from clumping <-- ##

## collate results
gwas.clump        <- dir("../tmpdir/")
gwas.clump        <- grep("clumps$", gwas.clump, value = T)
gwas.clump        <- rbindlist(lapply(gwas.clump, function(x) fread(paste0("../tmpdir/", x))))
## edit names
names(gwas.clump) <- gsub("\\#", "", names(gwas.clump))

## collate SNPs missing
miss.clump        <- dir("../tmpdir/")
miss.clump        <- grep("missing", miss.clump, value = T)
miss.clump        <- rbindlist(lapply(miss.clump, function(x) fread(paste0("../tmpdir/", x), header = F)), fill = T)
## add information to gwas catalog file
gwas.catalogue[, miss.clumping := SNPS %in% miss.clump$V1]
table(gwas.catalogue$miss.clumping)
#  FALSE   TRUE 
# 610672  14580 

## create clump groups
gwas.clump        <- as.data.table(gwas.clump)
clump.tmp         <- gwas.clump[, .(ld.buddies = unlist(strsplit(SP2, ","))), by = "ID"]
names(clump.tmp)  <- c("lead.snp.locus.r1", "ld.buddy")

## amend to GWAS catalog: two-step process
gwas.catalogue    <- merge(gwas.catalogue, clump.tmp, by.x = "SNPS", by.y = "ld.buddy", all.x = T)
gwas.catalogue[, lead.snp.locus.lead := SNPS %in% clump.tmp$lead.snp.locus.r1]
## harmonize
gwas.catalogue[, lead.snp.locus.r1 := ifelse(lead.snp.locus.lead == T, SNPS, lead.snp.locus.r1)]

## add information for embedding: r2 > 0.1
tmp.emb           <- res.complete[, .(tmp.r1 = unlist(strsplit(id.gwas.r1, "\\|"))), by = "emb"]
tmp.emb           <- tmp.emb[, .(embedding.r1 = paste(sort(unique(emb)), collapse = "|")), by = "tmp.r1"]
gwas.catalogue    <- merge(gwas.catalogue, unique(tmp.emb[ embedding.r1 != ""]), by.x = "SNPS", by.y = "tmp.r1", all.x = T)
## add information for embedding: r2 > 0.8
tmp.emb           <- res.complete[, .(tmp.r8 = unlist(strsplit(id.gwas.r8, "\\|"))), by = "emb"]
tmp.emb           <- tmp.emb[, .(embedding.r8 = paste(sort(unique(emb)), collapse = "|")), by = "tmp.r8"]
gwas.catalogue    <- merge(gwas.catalogue, unique(tmp.emb[ embedding.r8 != ""]), by.x = "SNPS", by.y = "tmp.r8", all.x = T)

#-----------------------------------#
##-- enrichment by mapped trait  --##
#-----------------------------------#

## --> convert for enrichment <-- ##

## adopt ontology mapping
gwas.catalogue[, term := sub(".*/([^/]+)$","\\1",sub(",.*","",MAPPED_TRAIT_URI))]
gwas.catalogue[, prefix := toupper(sub("_.*","",term))]

## prune accordingly: include only outcomes listed at least 5 times in the GWAS catalog
gwas.enr          <- table(gwas.catalogue$MAPPED_TRAIT)

## Define a helper to clean and collapse pipe-separated strings
clean_pipe        <- function(x) {
  # 1. Split by pipe
  # 2. Unlist into a single vector
  # 3. Remove NAs and empty strings
  # 4. Get unique values
  # 5. Paste back together
  vals <- unlist(strsplit(as.character(x), "\\|"))
  vals <- vals[!is.na(vals) & vals != "NA" & vals != ""]
  if(length(vals) == 0) return(NA_character_)
  paste(sort(unique(vals)), collapse = "|")
}

## aggregate accordingly: exclude MHC region!
gwas.enr          <- gwas.catalogue[ 
  !is.na(lead.snp.locus.r1) & MAPPED_TRAIT %in% names(gwas.enr[gwas.enr >= 10]) & !(CHR_ID == 6 & GENPOS.numeric.hg38 >= 28510020 - 1e6 & GENPOS.numeric.hg38 <= 33480577 + 1e6), 
  .(
    mapped_trait = paste(unique(MAPPED_TRAIT), collapse = "|"),
    mapped_term  = paste(unique(term), collapse = "|"),
    embedding.r1 = clean_pipe(embedding.r1),
    embedding.r8 = clean_pipe(embedding.r8)
  ), 
  by = "lead.snp.locus.r1"
]

## --> perform enrichment <-- ##

# 1. Get the total list of unique SNP loci (The "Universe")
# This must come from the original table before any filtering/splitting
all_loci          <- unique(gwas.enr[ !is.na(lead.snp.locus.r1)]$lead.snp.locus.r1)
total_n           <- length(all_loci)

# 2. Identify which SNPs have which Traits (Long format)
trait_map         <- gwas.enr[!is.na(lead.snp.locus.r1), .(
  trait = unlist(strsplit(as.character(mapped_trait), "\\|"))
), by = lead.snp.locus.r1]

# 3. Identify which SNPs have which Embeddings (Long format)
emb_map           <- gwas.enr[!is.na(lead.snp.locus.r1), .(
  emb = unlist(strsplit(as.character(embedding.r1), "\\|"))
), by = lead.snp.locus.r1]

# Remove empty strings/NAs from maps
trait_map         <- trait_map[!is.na(trait) & trait != ""]
emb_map           <- emb_map[!is.na(emb) & emb != ""]

# Get unique lists for looping
all_traits        <- unique(trait_map$trait)
all_embs          <- unique(emb_map$emb)

## run in parallel
registerDoMC(10)

## enrichment testing
enr.embedding     <- rbindlist(mclapply(all_traits, function(t) {
  
  # IDs that have Trait T
  ids_with_trait <- trait_map[trait == t, unique(lead.snp.locus.r1)]
  
  ## apply to each embedding
  rbindlist(lapply(all_embs, function(e) {
    
    # IDs that have Embedding E
    ids_with_emb <- emb_map[emb == e, unique(lead.snp.locus.r1)]
    
    # Calculate the 4 cells of the Fisher Matrix
    # a: Both Trait and Embedding
    a <- length(intersect(ids_with_trait, ids_with_emb))
    
    # b: Trait ONLY (Trait yes, Emb no)
    b <- length(setdiff(ids_with_trait, ids_with_emb))
    
    # c: Embedding ONLY (Emb yes, Trait no)
    c <- length(setdiff(ids_with_emb, ids_with_trait))
    
    # d: Neither (The true background)
    # Total - (Those with Trait OR those with Embedding)
    d <- total_n - length(union(ids_with_trait, ids_with_emb))
    
    # Fisher's Exact Test
    mat <- matrix(c(a, c, b, d), nrow = 2)
    ft  <- fisher.test(mat)
    
    data.table(
      trait        = t,
      embedding    = e,
      count_both   = a,
      n_trait      = length(ids_with_trait),
      n_emb        = length(ids_with_emb),
      intersection = paste(sort(intersect(ids_with_trait, ids_with_emb)), collapse = "|"),
      odds_ratio   = as.numeric(ft$estimate),
      p_val        = ft$p.value
    )
  }))
}, mc.cores = 10))

# Multiple testing correction
enr.embedding[, fdr := p.adjust(p_val, method = "BH")]

## --> prune results <-- ##

# 1. Pre-process SNP lists
enr.embedding[, snp_list := lapply(as.character(intersection), function(x) {
  if (is.na(x) || x == "" || x == "NA") return(character(0))
  unique(strsplit(x, "\\|")[[1]])
})]

## import function to do so
source("../functions/prune_consistent_by_factor.R")

# 2. Execute
enr.embedding.pruned <- prune_consistent_by_factor(enr.embedding, factor_col = "embedding", overlap_threshold = 0.1)

#-----------------------------------#
##-- enrichment by mapped trait  --##
##--       any embedding         --##
#-----------------------------------#

## enrichment testing
enr.embedding.loci <- rbindlist(mclapply(all_traits, function(t) {
  
  # IDs that have Trait T
  ids_with_trait <- trait_map[trait == t, unique(lead.snp.locus.r1)]
  
  # IDs that have Embedding E
  ids_with_emb <- emb_map[, unique(lead.snp.locus.r1)]
  
  # Calculate the 4 cells of the Fisher Matrix
  # a: Both Trait and Embedding
  a <- length(intersect(ids_with_trait, ids_with_emb))
  
  # b: Trait ONLY (Trait yes, Emb no)
  b <- length(setdiff(ids_with_trait, ids_with_emb))
  
  # c: Embedding ONLY (Emb yes, Trait no)
  c <- length(setdiff(ids_with_emb, ids_with_trait))
  
  # d: Neither (The true background)
  # Total - (Those with Trait OR those with Embedding)
  d <- total_n - length(union(ids_with_trait, ids_with_emb))
  
  # Fisher's Exact Test
  mat <- matrix(c(a, c, b, d), nrow = 2)
  ft  <- fisher.test(mat)
  
  data.table(
    trait        = t,
    count_both   = a,
    n_trait      = length(ids_with_trait),
    n_emb        = length(ids_with_emb),
    odds_ratio   = as.numeric(ft$estimate),
    intersection = paste(sort(intersect(ids_with_trait, ids_with_emb)), collapse = "|"),
    p_val        = ft$p.value
  )
  
}, mc.cores = 10))

# Multiple testing correction
enr.embedding.loci[, fdr := p.adjust(p_val, method = "BH")]

## --> prune results <-- ##

## import function to prune
source("../functions/prune_with_tracking.R")

# 1. Pre-process (Ensure snp_list is a list of character vectors)
enr.embedding.loci[, snp_list := lapply(as.character(intersection), function(x) {
  if (is.na(x) || x == "" || x == "NA") return(character(0))
  unique(strsplit(x, "\\|")[[1]])
})]

# 2. Sort by FDR and filter out empty rows
enr.embedding.loci        <- enr.embedding.loci[order(fdr)][lengths(snp_list) > 0]

# 4. Execute
enr.embedding.loci.pruned <- prune_with_tracking(enr.embedding.loci, overlap_threshold = 0.1)

###########################################
####     ICD-10 catalog overlap        ####
###########################################

#-------------------------------#
##--   prepare SNP look-up   --##
#-------------------------------#

## import SNP mapping: UKB imputed <-> WGS ICD10
snp.mapping    <- fread("<path_to_file>")

## clean data
grep("\\.x|\\.y", names(res.complete), value = T)
## clean accordingly
res.complete[, (grep("\\.y", names(res.complete), value = T)) := NULL]
## now the names
names(res.complete) <- gsub("\\.x", "", names(res.complete))

## create comparable identifier
res.complete[, marker_name_hg19 := paste(ifelse(CHROM == 23, "X", CHROM), GENPOS, pmin(ALLELE0, ALLELE1), pmax(ALLELE0, ALLELE1), sep = "_")]

## subset to SNPs included here: missing SNP rs11384540; multi-allelic allele
snp.mapping    <- snp.mapping[ marker_name_hg19 %in% res.complete$marker_name_hg19] 
## look at duplicates
jj             <- table(snp.mapping$marker_name_hg19)
snp.mapping[ marker_name_hg19 %in% names(jj[ jj > 1])]
## remove manually
snp.mapping    <- snp.mapping[ alt_hg19 == alt_hg38]
## add
res.complete   <- merge(res.complete, unique(snp.mapping[, .(marker_name_hg19, marker_name_hg38)]), by = "marker_name_hg19",
                        all.x = T)
## three not mappable

## create look-up file, including LD-proxies
r2.proxies[, snp.id_A := paste(CHROM_A, pos_hg38_A, sep = ":")]
r2.proxies[, snp.id_B := paste(CHROM_A, pos_hg38_B, sep = ":")]
## write all SNPs to file
write.table(c("ID", unique(c(r2.proxies[UNPHASED_R2 >= .5]$snp.id_A, r2.proxies[UNPHASED_R2 >= .5]$snp.id_B))), "Lookup.SNPs.WGS.stats.20260421.txt", 
            sep = "\t", row.names = F, quote = F, col.names = F)

## table of all possible outcomes
tmp <- dir("<path_to_file>")
write.table(tmp, "icd10.gwas.stats", sep = "\t", row.names = F, quote = F, col.names = F)

#-------------------------------#
##--   collate SNP look-up   --##
#-------------------------------#

## import results
snp.lookup.wgs <- dir("../wgs_lookup/")
snp.lookup.wgs <- rbindlist(mclapply(snp.lookup.wgs, function(x){
  ## read in
  tmp <- fread(paste0("../wgs_lookup/", x))
  ## add
  tmp[, icd10 := x]
  ## return results
  return(tmp)
}, mc.cores = 10), fill = T)
# ## get study accession
# snp.lookup.wgs[, study_accession := gsub("lookup\\.|\\.h*", "", icd10)]

## get the SNPs covered and decide on the proxies for each
snp.tmp <- unique(snp.lookup.wgs[, .(chromosome, base_pair_location, effect_allele, other_allele, ID)])
## create MarkerName
snp.tmp[, marker_name_hg38 := paste(ifelse(chromosome == 23, "X", chromosome), base_pair_location, pmin(effect_allele, other_allele), pmax(effect_allele, other_allele), sep = "_")]

## import mapping again
snp.mapping    <- fread("<path_to_file>")
## reduce to what is needed
snp.mapping    <- snp.mapping[marker_name_hg19 %in% res.complete$marker_name_hg19 | marker_name_hg38 %in% snp.tmp$marker_name_hg38]

## subset to SNPs also included in the imputed data
snp.tmp        <- merge(snp.tmp, unique(snp.mapping[, .(marker_name_hg19, marker_name_hg38, rsid)]), by = "marker_name_hg38")
## indicator which ones to keep: careful may still contain multi-allelic variants
snp.tmp[, keep.snp := marker_name_hg19 %in% res.complete$marker_name_hg19]

## the snps to keep
snp.keep.dir   <- snp.tmp[ keep.snp == T, .(marker_name_hg19, marker_name_hg38)]
## find proxies for the remaining ones
snp.keep.pro   <- r2.proxies[ ID_A %in% res.complete[ !(marker_name_hg19 %in% snp.keep.dir$marker_name_hg19)]$id & UNPHASED_R2 >= .5]
## add markername
snp.keep.pro[, CHROM_B := as.numeric(CHROM_B)]
snp.keep.pro   <- merge(snp.keep.pro, snp.tmp, 
                        by.x = c("CHROM_B", "pos_hg38_B"),
                        by.y = c("chromosome", "base_pair_location"),
                        all.x = T)
## order
snp.keep.pro   <- snp.keep.pro[ order(ID_A, -UNPHASED_R2)]
snp.keep.pro[, ind := 1:.N, by = "ID_A" ]
## subset accordingly
snp.keep.pro   <- snp.keep.pro[ ind == 1 & !is.na(marker_name_hg38), .(ID_A, marker_name_hg19, marker_name_hg38, UNPHASED_R2)]
snp.keep.pro   <- merge(snp.keep.pro, unique(res.complete[, .(id, marker_name_hg19)]),
                        by.x = "ID_A", by.y = "id", suffixes = c(".proxy", ".lead"))
## combine again
names(snp.keep.dir) <- c("marker_name_hg19.lead", "marker_name_hg38")
snp.keep.wgs   <- unique(rbind(snp.keep.dir, snp.keep.pro[, !"ID_A"], fill = T))
## fewer SNPs; due collapsing LD blocks on same SNP; some are simply missing

## add WGS ID
snp.keep.wgs   <- merge(snp.keep.wgs, snp.tmp[, .(marker_name_hg38, ID)])
## drop some multi-allelic variants
snp.keep.wgs   <- merge(snp.keep.wgs, snp.lookup.wgs[ ID %in% snp.lookup.wgs$ID, 
                                                      .(af.median = median(effect_allele_frequency)),
                                                      by = "ID"],
                        by = "ID")
## compute MAF
snp.keep.wgs[, maf.median := ifelse(af.median >= .5, 1-af.median, af.median)]
snp.keep.wgs    <- snp.keep.wgs[ order(marker_name_hg38, -af.median)]
snp.keep.wgs[, ind := 1:.N, by = "marker_name_hg38"]
snp.keep.wgs    <- snp.keep.wgs[ ind == 1]

## subset
snp.lookup.wgs  <- merge(snp.lookup.wgs, snp.keep.wgs[, .(ID, marker_name_hg19.lead)])
snp.lookup.wgs[, study_accession := sub(".*\\.(GCST[0-9]+)\\..*", "\\1", icd10)]
## add R2 group
snp.lookup.wgs  <- merge(snp.lookup.wgs, unique(res.complete[, .(marker_name_hg19, R2.group)]),
                         by.x = "marker_name_hg19.lead", by.y = "marker_name_hg19")
## matched 308 out of 318 SNPs

## import ICD10 annotation
icd10.map       <- fread("<path_to_file>")
## merge
snp.lookup.wgs  <- merge(snp.lookup.wgs, icd10.map[, !c("auc", "n_snps")], by.x = "study_accession", by.y = "Study Accession")

###########################################
####  cell type and tissue enrichment  ####
###########################################

#-------------------------------#
##--        HPA data         --##
#-------------------------------#

## import HPA v25 data
hpa.tissue        <- fread("<path_to_file>")
hpa.cells         <- fread("<path_to_file>")

## rename
names(hpa.tissue) <- c("ensembl_id", "gene", "tissue", "expression")
names(hpa.cells)  <- c("ensembl_id", "gene", "cell_type", "expression")

#-------------------------------#
##--    perform enrichment   --##
#-------------------------------#

## import function to do so
source("../functions/gene_enrichment_hpa.R")

## do in parallel
registerDoMC(10)

## run for one embedding: tissues
enrich.tissue <- rbindlist(mclapply(unique(res.complete$emb), function(x){
  
  ## get the genes to query
  genes <- na.omit(sapply(unique(res.complete[ emb == x]$hgnc_symbol), function(x) strsplit(x, "\\|")[[1]][1]))
  
  ## do only if enough genes
  if(length(genes) > 2){
    ## run enrichment
    enr <- test_expression_enrichment(hpa.tissue, 
                                      query_genes = genes,
                                      gene_col = "gene",
                                      group_col = "tissue",
                                      expr_col = "expression",
                                      min_query_n = 2,
                                      specificity_fold = 5,
                                      specificity_expr_min = 1,
                                      n_perm = 1e3)
    ## take only forward what is needed
    enr <- enr$summary
    ## add the embedding tested
    enr[, emb := x]
    ## return results of interest
    return(enr)
  }
  
}, mc.cores = 10))
## only liver for two embeddings

## run for one embedding: cell types
enrich.cells  <- rbindlist(mclapply(unique(res.complete$emb), function(x){
  
  ## get the genes to query
  genes <- na.omit(sapply(unique(res.complete[ emb == x]$hgnc_symbol), function(x) strsplit(x, "\\|")[[1]][1]))
  
  ## do only if enough genes
  if(length(genes) > 2){
    ## run enrichment
    enr <- test_expression_enrichment(hpa.cells, 
                                      query_genes = genes,
                                      gene_col = "gene",
                                      group_col = "cell_type",
                                      expr_col = "expression",
                                      log_fold = T,
                                      min_query_n = 2,
                                      specificity_fold = 1.5,
                                      specificity_expr_min = 5,
                                      n_perm = 1e3)
    ## take only forward what is needed
    enr <- enr$summary
    ## add the embedding tested
    enr[, emb := x]
    ## return results of interest
    return(enr)
  }
  
}, mc.cores = 10))

###########################################
####       sex-stratified results      ####
###########################################

#-------------------------------------#
##--         import results        --##
#-------------------------------------#

## import results
res.regional.sex <- dir("<path_to_file>")
res.regional.sex <- rbindlist(lapply(res.regional.sex, function(x){
  
  ## import results
  tmp <- fread(paste0("<path_to_file>", x))
  ## add embedding
  tmp[, emb := gsub("\\.sexdiff\\.final\\.tsv", "", x)]
  ## return results
  return(tmp)
  
}))

#------------------------------------#
##--  reference to prev. results  --##
#------------------------------------#

## export list of SNPs: include both, all, sex-specific and MHC results
write.table(c(unique(res.credible$id), res.regional.sex$MarkerName), "snp.list.embeddings.txt", col.names = F, row.names = F, quote = F)

## run: no X-chromosome!
system("./obtain_snps.sh")

## import
snp.dat        <- fread("../tmpdir/snp.dosage.transpose")
snp.info       <- fread("../tmpdir/snp.info")
## assign names
names(snp.dat) <- c("f.eid", snp.info$rsid)

## subset to unrelated EUR population
tmp.dat        <- fread("<path_to_file>")
## apply
snp.dat        <- snp.dat[ f.eid %in% tmp.dat$V1]

## compute LD matrix
snp.ld         <- Rfast::cora(snp.dat[,-1])
## convert to data table (keep only entries from upper triangle)
snp.ld         <- melt(
  as.data.table(snp.ld, keep.rownames="var1"),
  id.vars="var1"
)[, `:=`(var1=as.character(var1), var2=as.character(variable))
][var1 <= var2]
## convert to r2
snp.ld[, value := value^2]

## convert to graph (use max to allow for edge weights)
ld.sub           <- graph_from_data_frame(snp.ld[ value >= .7])
## get all separate components
ld.sub           <- components(ld.sub)$membership
## convert to data frame
ld.sub           <- data.table(ID=names(ld.sub), R2.group=ld.sub)
## add previous assignment
ld.sub           <- merge(ld.sub, unique(res.complete[, .(id, R2.group)]),
                          by.x = "ID", by.y = "id", all = T, suffixes = c(".updated", ".previous"))

## test whether selected SNPs are in LD clumps with overall findings
res.regional.sex <- merge(res.regional.sex, ld.sub, by.x = "MarkerName", by.y = "ID",
                          all.x = T)

## additional whether the same variant was reported
res.regional.sex[, match.complete.results := paste(MarkerName, emb) %in% paste(res.complete$id, res.complete$emb)]
## n = 3

#------------------------------------#
##--            summary           --##
#------------------------------------#

## sex-dimorphic effects
nrow(res.regional.sex[ P_male > .05 | P_female > .05])
## n = 13
nrow(res.regional.sex[ P_male > .05])
nrow(res.regional.sex[ P_female > .05])

## sex-dimorphic effects: from overall results
nrow(res.regional.sex[ (P_male > .05 | P_female > .05) & (match.complete.results == T | !is.na(R2.group.previous))])

## look into examples
View(res.regional.sex[ P_male > .05 | P_female > .05])

## write to file
write.table(res.regional.sex, "Results.sex.stratified.analysis.UKB.embeddings.20260427.txt", sep = "\t", row.names = F)

###########################################
####        SNP-based heritability     ####
###########################################

## import results
res.ldsc <- fread("ldsc_summary_sorted_pval.tsv")

## reporting for the paper
summary(res.ldsc$h2)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00670 0.01427 0.01915 0.02097 0.02473 0.06350 

## recompute p-value for sanity checking
res.ldsc[, p_h2      := pnorm(z, lower.tail = FALSE)] 

## how many significant
nrow(res.ldsc[ pval < .05/nrow(res.ldsc)])

###########################################
####     numbers for the manuscript    ####
###########################################

#------------------------------#
##--    generic reporting   --##
#------------------------------#

## create p-value column
res.complete[, PVAL := 10^-LOG10P]
## add MAF to the results
res.complete[, MAF := ifelse(A1FREQ > .5, 1 - A1FREQ, A1FREQ)]

## how many passing more stringent statistical significance: p<4.166667e-10
nrow(res.complete[ LOG10P > -log10(5e-8/120)])
## n = 419

## how many embeddings with at least one variant associated
length(unique(res.complete$emb))
## n = 100

## how many variants per embedding
summary(as.vector(table(res.complete$emb)))
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00    1.00    2.00    4.429    6.00   34.00 

## top embedding
tail(sort(table(res.complete$emb)), 10)

## how many with low PIP
nrow(res.complete[ pip > .5])
## n = 100

## how many variants associated with more than one embedding
length(unique(res.complete$R2.group))-1
## how many associate with more than one embedding
jj <- table(res.complete[ R2.group != 0]$R2.group)
sum(jj > 1)
## n = 54
sort(table(res.complete[ R2.group != 0]$R2.group))

## look at most shared signals
View(res.complete[ R2.group %in% as.numeric(names(jj[ jj >= 10]))])

## look at specific signals; strong association with one but not others
View(res.complete[ R2.group %in% as.numeric(names(jj[ jj == 1]))])

#------------------------------#
##--     gene enrichment    --##
#------------------------------#

## --> pathway enrichment <-- ##

## import results from Wenhuan
res.pathway <- dir("<path_to_file>")
res.pathway <- rbindlist(lapply(res.pathway, function(x){
  ## import results
  res <- fread(paste0("<path_to_file>", x))
  ## add embedding
  res[, emb := gsub("_pathway_fdr.tsv", "", x)]
  ## return results
  return(res)
}))
## drop ubiquitous terms
res.pathway <- res.pathway[ term_size < 500]
## add fc 
res.pathway[, fc := (intersection_size/term_size)/(query_size/effective_domain_size)]

## import function to prune findings
source("../functions/pathway_pruning.R")

## apply
res.pathway.pruned <- prune_pathways(res.pathway, jaccard_threshold = .1, rank_by = "new_genes")
View(res.pathway.pruned[ intersection_size > 2])
write.table(res.pathway[ intersection_size > 2], "Pathway.enrichment.effector.gene.by.embedding.20260422.txt", sep = "\t", row.names = F)

## how many embeddings with at least one pathway
length(unique(res.pathway[ intersection_size > 2]$emb))
## n = 15

## how many embeddings with cholesterol metabolism

## --> tissue enrichment <-- ##

## look into persistent results
enrich.tissue[ fisher_padj < .05 & fisher_OR > 1 & n_specific_query > 2]
## n = 11

## export results
write.table(enrich.tissue[ fisher_padj < .05 & fisher_OR > 1 & n_specific_query > 2], "Results.tissue.enrichment.effector.genes.UKB.embeddings.txt", 
            sep = "\t", row.names = F)

# ## import results
# res.tissue <- dir("<path_to_file>")
# res.tissue <- rbindlist(lapply(res.tissue, function(x){
#   ## import 
#   res <- fread(paste0("<path_to_file>", x))
#   ## add embedding
#   res[, emb := gsub("_hpa_tissue_enrichment.tsv", "", x)]
#   ## return results
#   return(res)
# }))
# 
# ## sig results
# res.tissue[ or > 1 & fdr < .05 & d1 > 2]
# ## only liver for embedding 92

## --> cell type enrichment <-- ##

## look into persistent results
enrich.cells[ fisher_padj < .05 & fisher_OR > 1 & n_specific_query > 2]
## n = 10

## export results
write.table(enrich.cells[ fisher_padj < .05 & fisher_OR > 1 & n_specific_query > 2], "Results.cell.type.enrichment.effector.genes.UKB.embeddings.txt", 
            sep = "\t", row.names = F)


# ## import results
# res.cell.type <- dir("<path_to_file>")
# res.cell.type <- rbindlist(lapply(res.cell.type, function(x){
#   ## import 
#   res <- fread(paste0("<path_to_file>", x))
#   ## add embedding
#   res[, emb := gsub("_hpa_cell_type_enrichment.tsv", "", x)]
#   ## return results
#   return(res)
# }))
# 
# ## sig results
# res.cell.type[ or > 1 & fdr < .05 & d1 > 2]
# ## only liver for embedding 92



#------------------------------#
##--     GWAS enrichment    --##
#------------------------------#

## reference to GWAS catalog
nrow(res.complete[ num_reported.r1 > 0])
## all

## high confidence
length(unique(res.complete[ num_reported.r8 > 0]$R2.group))

## how many embeddings have at least one enrichment passing significance
length(table(enr.embedding[ fdr < .05 & count_both >= 3]$embedding))
## n = 51

## write results to file
write.table(enr.embedding[ count_both >= 3 & fdr < .05, !"snp_list"], "Results.Embedding.GWAS.loci.enrichment.20260421.txt", sep = "\t", row.names = F)

## look up examples: most promiscuous embedding after pruning
tail(sort(table(enr.embedding.pruned[ count_both >= 3 & fdr < .05]$embedding)))

## how many embeddings enriched for 'Hypercholesterolemia'
nrow(enr.embedding[ trait == "Hypercholesterolemia" & fdr < .05 & odds_ratio > 1])
## how many embeddings enriched for 'asthma'
nrow(enr.embedding[ trait == "asthma" & fdr < .05 & odds_ratio > 1])
## how many embeddings enriched for 'educational attainment'
nrow(enr.embedding[ trait == "educational attainment" & fdr < .05 & odds_ratio > 1])
## how many embeddings enriched for 'educational attainment'
nrow(enr.embedding[ trait == "body mass index" & fdr < .05 & odds_ratio > 1])

#-------------------------------#
##--      WGS SNP look-up    --##
#-------------------------------#

## how many r2 groups with at least on sig. finding (excluding MHC region)
uniqueN(snp.lookup.wgs[ R2.group != 0 & p_value < 5e-8]$R2.group)
## n = 123 (p<5e-8); n = 137 (p<1e-3)
uniqueN(snp.lookup.wgs[ R2.group != 0 ]$R2.group)
## n = 123/144

## snps with no finding
tmp.r2group <- snp.lookup.wgs[, .(num.gws.sig = uniqueN(study_accession[p_value < 5e-8]),
                                  num.int.sig = uniqueN(study_accession[p_value < 1e-3]),
                                  num.nom.sig = uniqueN(study_accession[p_value < .05])),
                              by = "R2.group"]
## strong effects
uniqueN(snp.lookup.wgs[ p_value < 5e-8 & R2.group != 0 & abs(beta) > log(1.5)]$R2.group)
## n = 11
View(snp.lookup.wgs[ p_value < 5e-8 & R2.group != 0 & abs(beta) > log(1.5)])

## no effects: look into r2 groups
res.complete[ R2.group %in% tmp.r2group[ num.nom.sig == 0]$R2.group, .(R2.group, id, emb, BETA, SE, LOG10P, pip, num_reported.r1, num_reported.r8, trait_reported.r8, nearest_gene_name, hgnc_symbol)]

## write to file
fwrite(snp.lookup.wgs[ p_value < 5e-8], "Results.WGS.look.up.ICD10.codes.20260630.txt", sep = "\t", row.names = F, na = NA)

## import results from association testing
res.assoc <- fread("../../02_association_analysis/input/Results.Embedding.associations.UKB.minimal.extensive.20260319.txt.gz")

#------------------------------#
##--  export for reporting  --##
#------------------------------#

## write to file
write.table(res.complete, "Lead.credible.set.variants.incl.MHC.Delphi.embeddings.UKB.GWAS.plus.catalog.gene.annotation.20260421.txt", sep = "\t", row.names = F)

#------------------------------#
##--    BioRender figure    --##
#------------------------------#

## number embeddings with GSDMB as putative effector gene
uniqueN(res.complete[ grep("GSDMB|ORMDL3", hgnc_symbol)]$emb)
## number embeddings with GSDMB as putative effector gene
uniqueN(res.complete[ grep("IL33", hgnc_symbol)]$emb)
## number embeddings with TSPL as putative effector gene
uniqueN(res.complete[ grep("TSLP", hgnc_symbol)]$emb)
## number embeddings with TSPL as putative effector gene
uniqueN(res.complete[ grep("IL1RL1", hgnc_symbol)]$emb)
## number embeddings with TSPL as putative effector gene
uniqueN(res.complete[ grep("IL18R1", hgnc_symbol)]$emb)
## number embeddings with LRCC32 as putative effector gene
uniqueN(res.complete[ grep("LRRC32", hgnc_symbol)]$emb)
## number embeddings with SMAD3 as putative effector gene
uniqueN(res.complete[ grep("SMAD3", hgnc_symbol)]$emb)

## number embeddings with LPA as putative effector gene
uniqueN(res.complete[ grep("LPA", hgnc_symbol)]$emb)
## number embeddings with APOB as putative effector gene
uniqueN(res.complete[ grep("APOB", hgnc_symbol)]$emb)
## number embeddings with PCK9 as putative effector gene
uniqueN(res.complete[ grep("PCSK9", hgnc_symbol)]$emb)
## number embeddings with LDLR as putative effector gene
uniqueN(res.complete[ grep("LDLR", hgnc_symbol)]$emb)
## number embeddings with NPC1L1 as putative effector gene
uniqueN(res.complete[ grep("NPC1L1", hgnc_symbol)]$emb)
## number embeddings with LRCC32 as putative effector gene
uniqueN(res.complete[ grep("ABCG5|ABCG8", hgnc_symbol)]$emb)
## number embeddings with HMGCR as putative effector gene
uniqueN(res.complete[ grep("HMGCR", hgnc_symbol)]$emb)


###########################################
####     figures for the manuscript    ####
###########################################

#------------------------------#
##--     embQTL discovery   --##
#------------------------------#

## --> create chromosome mapping <-- ##

## import build 37: from wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
chr.dat        <- fread("hg19.chrom.sizes")
## change names
names(chr.dat) <- c("chromosome", "length")
## create chromosome
chr.dat[, CHR_numeric := as.numeric(ifelse(chromosome == "chrX", 23, gsub("chr", "", chromosome)))]
## drop NA
chr.dat        <- chr.dat[ !is.na(CHR_numeric)]
## order
chr.dat        <- chr.dat[ order(CHR_numeric)]
## add position in the plot
chr.dat[, chr_offset := cumsum(as.numeric(shift(length, fill = 0)))]

## add to credible set results
res.complete   <- merge(res.complete, chr.dat[, .(CHR_numeric, chr_offset)], 
                        by.x = "CHROM",
                        by.y = "CHR_numeric")
## position in the plot
res.complete[, plt.pos := chr_offset + GENPOS]

## --> embedding cluster <-- ##

## all pairs of emb entries sharing at least one R2.group
dt.edges     <- res.complete[, {
  if(.N > 1) {
    tmp <- combn(emb, 2, simplify = F)
    data.table(
      from = sapply(tmp, `[`, 1),
      to   = sapply(tmp, `[`, 2)
    )
  }
}, by = R2.group]

## count shared R2.groups per pair as edge weight and collapse
dt.edges[, weight := .N, by = .(from, to)]
dt.edges     <- unique(dt.edges[, .(from, to, weight)])

## build undirected graph from weighted edge list
g            <- graph_from_data_frame(dt.edges, directed = F, vertices = unique(res.complete[, .(emb)]))
cl           <- cluster_louvain(g, weights = E(g)$weight)

## get the clustering
emb.clusters <- data.table(
  emb         = cl$names,
  emb.cluster = cl$membership
)
## add number of embedding associated loci
emb.clusters[, num.embeddings := sapply(emb, function(x) nrow(res.complete[ emb == x]))]

## summary
table(emb.clusters$emb.cluster)
#  1  2  3  4  5 
# 24 54 16  1  3 

## pathway enrichment by cluster (based on included loci)
cluster.gene.enr <- rbindlist(lapply(1:5, function(x){
  
  ## get all possible genes by r2 group; drop MHC
  tmp <- res.complete[ emb %in% emb.clusters[ emb.cluster == x]$emb & R2.group != 0, .(R2.group, hgnc_symbol) ]
  
  ## parse accordingly
  tmp <- lapply(unique(tmp$R2.group), function(k){
    ## get all genes
    genes <- unlist(lapply(tmp[ R2.group == k]$hgnc_symbol, function(c) strsplit(c, "\\|")[[1]]))
    ## count
    genes <- sort(table(genes), decreasing = T)
    ## return top
    return(names(genes[1]))
  })
  ## get gene list
  tmp <- unique(unlist(tmp))
  
  ## only if sufficient number of genes
  if(length(tmp) > 0){
    ## enrichment testing
    enr <- gprofiler2::gost(query = tmp, 
                            organism = "hsapiens", ordered_query = FALSE, 
                            multi_query = FALSE, significant = TRUE, exclude_iea = TRUE,
                            measure_underrepresentation = FALSE, evcodes = TRUE, 
                            user_threshold = 0.05, correction_method = "fdr", 
                            domain_scope = "annotated",
                            numeric_ns = "", sources = c("KEGG", "REAC"), as_short_link = FALSE)
    
    ## return results
    return(data.table(emb.cluster = x,
                      enr$result))
  }else{
    return(data.table(emb.cluster = x))
  }
}))
## how many findings
table(cluster.gene.enr$emb.cluster)

## get top pick for each
cluster.gene.enr <- cluster.gene.enr[ order(emb.cluster, p_value)]
## add
cluster.gene.enr[, ind := 1:.N, by = emb.cluster]
cluster.gene.enr[ ind == 1, .(emb.cluster, p_value, term_name, intersection)]

## add sensible annotations
tmp.anno         <- data.table(emb.cluster      = 1:5,
                               emb.anno         = c("Cholesterol metabolism", "MHC region", "Interleukin-33 signaling", "LDL remodelling",
                                                    "ROBO signalling"),
                               emb.cluster.sort = c(1,5,3,2,4),
                               emb.num          = sapply(1:5, function(x) nrow(emb.clusters[ emb.cluster == x])))
## add
emb.clusters     <- merge(emb.clusters, tmp.anno, by = intersect(names(emb.clusters),
                                                                 names(tmp.anno)))
## define order for plotting
emb.clusters     <- emb.clusters[ order(emb.cluster.sort, -num.embeddings)]
emb.clusters[, emb.pos := 1:.N]

## add to results
res.complete     <- merge(res.complete, emb.clusters, by = intersect(names(res.complete), 
                                                                     names(emb.clusters)))

## open PDF
pdf("../graphics/Summary.credible.sets.2D.Manhattan.20260626.pdf", width = 6.3, height = 6.3)
## define plotting parameters
par(mar=c(.2,2.5,1.5,.2), tck = -.01, mgp = c(.8,.2,0), cex.axis = .7, cex.lab = .7, xaxs = "i", yaxs = "i",
    bty = "l", lwd = .5)

## establish layout
layout(matrix(1:4, 2, 2, byrow = T), heights = c(.15,.85), widths = c(.85,.15))

## --> locus count <-- ##

## what to plot: count and gene
tmp <- res.complete[, .(locus.count = uniqueN(emb),
                        plt.pos = min(plt.pos),
                        gene = paste(unique(hgnc_symbol), collapse = "|")),
                    by = R2.group]

## empty plot
plot(c(0, chr.dat$chr_offset[ nrow(chr.dat)] + chr.dat$length[nrow(chr.dat)]), c(0, log10(85)),
     type = "n", ylab = "#Embeddings", xlab = "", xaxt = "n", yaxt = "n")
## add axis
axis(2, lwd=.5, at = log10(c(1,5,10,50)), labels = c(1,5,10,50))

## add hits
arrows(tmp$plt.pos, 0, tmp$plt.pos, log10(tmp$locus.count), lwd=.5, length = 0)
## add pins
points(tmp$plt.pos,  log10(tmp$locus.count), pch=21, cex=.6, lwd=.5, bg = "#E8601C", xpd = NA)

## add top genes
tmp.label <- tmp[locus.count > 5]

## shorten gene label: keep only unique gene symbols, drop NAs,
## cap at 3 genes to avoid overly long strings at dense loci;
## skip cleaning for manually set labels (e.g. MHC)
tmp.label[, gene.label := sapply(gene, function(x) {
  if(x %in% c("MHC")) return(x)
  g <- unique(unlist(strsplit(x, "\\|")))
  g <- g[g != "NA"]
  if(length(g) > 3) g <- c(g[1:3], "...")
  paste(g, collapse = "|")
})]

## add labels with white background box for legibility
text(tmp.label$plt.pos, log10(tmp.label$locus.count) + 0.05,
     labels = tmp.label$gene.label,
     cex    = .6,
     srt    = 45,
     xpd    = NA,
     adj    = c(0, 0.5),
     font   = 2)

## --> empty plot <-- ## 

par(bty = "n")

plot(1, 1, type = "n", ylab = "", xlab = "", xaxt = "n", yaxt = "n")


## --> 2D Manhattan <-- ##

## adopt plotting margins
par(mar=c(2.5,2.5,.2,.2))

## empty plot
plot(c(0, chr.dat$chr_offset[ nrow(chr.dat)] + chr.dat$length[nrow(chr.dat)]), c(0,99),
     type = "n", ylab = "", xlab = "Chromosomal position", xaxt = "n", yaxt = "n", ylim = rev(c(0,99)))
## add axis
axis(1, lwd=.5, at=chr.dat$chr_offset + (chr.dat$length)/2, labels = chr.dat$CHR_numeric)

## get plotting coordinates
pm <- par("usr")

## stratify by chromosome
rect(chr.dat$chr_offset, pm[3], chr.dat$chr_offset + chr.dat$length, pm[4], border = NA,
     col = c("white", "grey90"))

## stratify embeddings
# abline(h=1:max(res.complete$emb.pos)-.5, lwd=.5, lty = 2, col = "grey30")
abline(h=res.complete[, .( min.pos = min(emb.pos)), by = emb.cluster]$min.pos-.5, lwd=.5, lty = 2, col = "grey30")

## add points
points(res.complete$plt.pos, res.complete$emb.pos, cex = .6, pch=21, lwd=.3, bg = "grey50",
       col="white", xpd = NA)

## add names
text(pm[1], emb.clusters$emb.pos, pos = 2, xpd = NA,
     cex = .5, labels = emb.clusters$emb)

## --> Number of loci <-- ##

## adopt plotting margins
par(mar=c(2.5,.2,.2,.5))

## empty plot
plot(c(0, max(emb.clusters$num.embeddings)), c(0,99),
     type = "n", ylab = "", xlab = "#Loci", 
     xaxt = "n", yaxt = "n", ylim = rev(c(0,99)))
## add axis
axis(1, lwd=.5)

## add counts
rect(0, emb.clusters$emb.pos-.4, emb.clusters$num.embeddings, emb.clusters$emb.pos+.4,
     border = NA, col = "#4EB3D3")

## close device
dev.off()


#------------------------------#
##--     GWAS enrichment    --##
#------------------------------#

## --> prep for plotting <-- ##

## omit associations with poor overlap
tmp.plot   <- unique(enr.embedding.pruned[ count_both >= 5 & fdr < .05]$trait)
## keep only those traits for visualization
tmp.plot   <- enr.embedding[ trait %in% tmp.plot]

## create numbering for traits
num.traits <- data.table(trait = unique(tmp.plot$trait), 
                         num.trait = 1:length(unique(tmp.plot$trait)))
## add
tmp.plot   <- merge(tmp.plot, num.traits)

## --> create ordering of embeddings <-- ##

# 1. Transform FDR to -log10 for better clustering contrast
tmp.plot[, log_fdr := -log10(fdr)]

# 2. Pivot to Wide Format (Traits x Embeddings)
# Fill missing associations with 0
tmp.mat           <- dcast(tmp.plot, trait ~ embedding, 
                           value.var = "log_fdr", 
                           fill = 0)

# 3. Convert to a standard matrix for clustering
row_names         <- tmp.mat$trait
tmp.mat           <- as.matrix(tmp.mat [, -1, with = FALSE])
rownames(tmp.mat) <- row_names

# 4. Extract Embedding Clusters
# Compute the distance matrix for columns (embeddings)
# emb_hclust        <- hclust(as.dist(1 - cor(tmp.mat, method = "pearson")), method = "ward.D2")
emb_hclust        <- hclust(proxy::dist(t(tmp.mat), method = "cosine"), method = "ward.D2")

## 5. map to a file to annotate the plotting
num.emb    <- data.table(embedding = emb_hclust$labels,
                         num.emb = emb_hclust$order)
## add
tmp.plot   <- merge(tmp.plot, num.emb, by = "embedding")

## --> actual plot <-- ##

## create PDF
pdf("../graphics/Heatmap.GWAS.locus.enrichment.by.embedding.20260421.pdf", width = 6.3, height = 2)
## graphical parameters
par(mar=c(2,4,1.5,.5), mgp = c(.6,0,0), cex.axis = .5, cex.lab = .5, tck=-.01, lwd=.5, xaxs = "i", yaxs = "i")

## colour gradient
col.vec <- colorRampPalette(c("white", "red3"))(round(max(tmp.plot$log_fdr)+2))

## --> heatmap <-- ##

## emty plot
plot(c(.5, max(num.emb$num.emb)+.5), c(.5, max(num.traits$num.trait)+.5), type = "n", xlab = "", ylab = "",
     xaxt = "n", yaxt = "n")

## add p-value heatmap from enrichment
rect(tmp.plot$num.emb -.5, tmp.plot$num.trait -.5, tmp.plot$num.emb + .5, tmp.plot$num.trait + .5,
     col = col.vec[round(tmp.plot$log_fdr)+1], border = NA)

## add rectangles for those meeting significance
jj <- which(tmp.plot$fdr < .05)
rect(tmp.plot$num.emb[jj] -.5, tmp.plot$num.trait[jj] -.5, tmp.plot$num.emb[jj] + .5, tmp.plot$num.trait[jj] + .5,
     col = NA, border = "grey20", lwd=.5)

## plotting coordinates
pm <- par("usr")

## add names: embeddings
text(num.emb$num.emb, pm[3]-(pm[4]-pm[3])*.05, cex=.4, xpd=NA, labels = num.emb$embedding, offset = 0,
     srt=90, pos = 2)

## add names: traits
text(pm[1], num.traits$num.trait, cex=.5, xpd=NA, labels = stringr::str_to_title(num.traits$trait), offset = .1, pos = 2)

## box
box(lwd=.5)

## --> add color gradient <-- ##

## length
l <- seq(pm[1]+(pm[2]-pm[1])*.05, pm[1]+(pm[2]-pm[1])*.35, length.out = length(col.vec))
## rectangle for the colours
rect(l-(l[2]-l[1])/2, pm[4]+(pm[4]-pm[3])*.1, l+(l[2]-l[1])/2, pm[4]+(pm[4]-pm[3])*.15, border=NA, col=col.vec, xpd=NA)
## box
rect(l[1]-(l[2]-l[1])/2, pm[4]+(pm[4]-pm[3])*.1, l[length(l)]+(l[2]-l[1])/2, pm[4]+(pm[4]-pm[3])*.15, border="black", col=NA, lwd=.3, xpd=NA)
## add header
text(pm[1]+(pm[2]-pm[1])*.1, pm[4]+(pm[4]-pm[3])*.18, cex=.4, labels = "-log10(FDR)", pos=4,
     offset = .2, xpd=NA)
## simple axis
text(l[round(c(1, c(.2, .4, .6, .8, 1)*length(col.vec)))], pm[4]+(pm[4]-pm[3])*.08, 
     labels=sprintf("%.f", c(0, .2, .4, .6, .8, 1)*length(col.vec)), pos=1, cex=.3, offset = .1, xpd=NA)

dev.off()