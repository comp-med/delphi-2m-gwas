## shared functions for embedding x drug colocalisation
## Maik Pietzner 04/06/2026

## --> packages assumed loaded by the caller: data.table, coloc, arrow <-- ##

#-----------------------------------------#
##--   quick signal check (the gate)   --##
#-----------------------------------------#

## cheap test: does the embedding GWAS carry any signal in this region at all?
## if not, there is nothing for coloc to pick up against the drug instrument and
## we skip the (comparatively expensive) merge / coloc / LD / plotting entirely.
embedding.has.signal <- function(res.emb, log10p.min = 6){
  ## res.emb is the raw regional REGENIE slice for one embedding
  if(is.null(res.emb) || nrow(res.emb) == 0) return(FALSE)
  return(max(res.emb$LOG10P, na.rm = T) >= log10p.min)
}

#-----------------------------------------#
##--   one embedding x one drug locus  --##
#-----------------------------------------#

## does everything the original 03_run_drug_coloc.R body did for a single
## (embedding, disease, region) triplet: align to the drug GWAS, run coloc,
## add LD-based sensitivity + locus-compare plot where there is H4 evidence,
## write the per-triplet output file, and return the summary row invisibly.
run.drug.coloc <- function(res.gwas,      ## drug GWAS regional slice (hg38, A2 = effect)
                           res.emb,       ## embedding regional slice (raw REGENIE, hg19)
                           lift.coords,   ## hg19<->hg38 mapping for this region
                           embedding,
                           disease,
                           chr.s, pos.s, pos.e,
                           p12        = 5e-6,
                           h4.ld.thr  = .5,    ## run LD + plot only above this PP.H4
                           do.plot    = T){
  
  #-- map embedding variants to hg38 so they can meet the drug GWAS --#
  
  ## skip the liftover merge when the caller already attached the hg38 id
  ## (the workflow pre-lifts once in 05_prep_query_gwas.R); otherwise do it here
  if(!"marker_name_hg38_ucsc" %in% names(res.emb)){
    ## create hg19 identifier
    res.emb[, marker_name_hg19_ucsc := paste0("chr", CHROM, ":", GENPOS, "_", pmin(ALLELE0, ALLELE1), "_", pmax(ALLELE0, ALLELE1))]
    ## attach hg38 identifier
    res.emb     <- merge(res.emb, lift.coords[, .(marker_name_hg19_ucsc, marker_name_hg38_ucsc)],
                         by = "marker_name_hg19_ucsc")
    ## drop multi-allelic variants
    jj          <- table(res.emb$marker_name_hg19_ucsc)
    res.emb     <- res.emb[ names(jj[ jj == 1 ])]
  }
  
  #-- combine the two GWAS on their shared variant set --#
  
  res.comb    <- merge(res.gwas[, .(ID_UCSC, SNP, CHR_ENSEMBL, BP, A1, A2, BETA, SE, P)],
                       res.emb[, .(marker_name_hg38_ucsc, ID, GENPOS, ALLELE0, ALLELE1, A1FREQ, BETA, SE, LOG10P)],
                       by.x       = "ID_UCSC",
                       by.y       = "marker_name_hg38_ucsc",
                       suffixes   = c(".drug", ".embedding"))
  
  ## nothing in common -> bail out gracefully
  if(nrow(res.comb) == 0){
    cat("  no overlapping variants for", embedding, "vs", disease, "- skipping\n")
    return(invisible(NULL))
  }
  
  ## align effect direction + allele frequency for the embedding
  res.comb[, BETA.embedding := ifelse(ALLELE1 == A2, BETA.embedding, -BETA.embedding)]
  res.comb[, A2FREQ         := ifelse(ALLELE1 == A2, A1FREQ, 1 - A1FREQ)]
  
  #-- run colocalisation --#
  
  ## embedding (quantitative)
  D1            <- list(
    beta     = res.comb$BETA.embedding,
    varbeta  = res.comb$SE.embedding^2,
    type     = "quant",
    N        = max(res.emb$N, na.rm = T),
    sdY      = 1,
    snp      = res.comb$ID,
    position = res.comb$GENPOS
  )
  
  ## drug indication (binary; dummy s / N do not affect the single-variant abf)
  D2            <- list(
    beta     = res.comb$BETA.drug,
    varbeta  = res.comb$SE.drug^2,
    s        = .01,
    N        = 1e5,
    type     = "cc",
    snp      = res.comb$ID,
    position = res.comb$GENPOS
  )
  
  ## run coloc
  res.coloc     <- coloc.signals(D1, D2, method = "single", p12 = p12)
  
  ## SNPs of interest: tweak if a lead happens to fall outside the overlap
  snps.interest <- data.table(
    type.snp = c("embQTL", "regional.lead.drug", "top.shared"),
    ID       = c(res.coloc$summary$best1,
                 res.coloc$summary$best2,
                 res.coloc$summary$best4)
  )
  ## add effect estimates of interest
  snps.interest <- merge(snps.interest,
                         res.comb[, .(ID, ID_UCSC, CHR_ENSEMBL, GENPOS,
                                      BP,            ## position in build 38
                                      A1, A2,        ## A2 = effect allele
                                      A2FREQ, BETA.embedding, SE.embedding, LOG10P,
                                      BETA.drug, SE.drug, P)])
  
  #-- sensitivity (LD) only where there is coloc evidence --#
  
  if(res.coloc$summary$PP.H4.abf > h4.ld.thr){
    
    ## --> LD among respective lead signals <-- ##
    
    ## stem for the temporary SNP-list / LD files
    ld.stem   <- paste("tmpdir/snp.list", embedding, chr.s, pos.s, pos.e, disease, "txt", sep = ".")
    ## export SNP list and obtain the LD matrix via REGENIE
    write.table(res.comb$ID, ld.stem, sep = "\t", row.names = F, col.names = F, quote = F)
    system(paste("./obtain_ld_matrix.sh", ld.stem, chr.s))
    ## importer for the binary LD matrix
    source("../functions/import_ld_matrix_regenie.R")
    ## read (LD matrix - not for fine-mapping!)
    snplist   <- readLines(paste(ld.stem, "corr.snplist", sep = "."))
    ld_mat    <- get.corr.sq.matrix(paste(ld.stem, "corr", sep = "."))
    ## clean up
    system(paste("rm", paste("tmpdir/snp.list", embedding, chr.s, pos.s, pos.e, disease, "*", sep = ".")))
    ## LD with the lead embQTL
    snps.interest[, R2.lead.embQTL := ld_mat[snps.interest$ID, res.comb[ which.max(abs(BETA.embedding/SE.embedding))]$ID]]
    
  }
  
  #-- assemble + write the per-triplet output --#
  
  res.out   <- data.table(
    embedding         = embedding,
    region_start.hg38 = pos.s,
    region_end.hg38   = pos.e,
    disease           = disease,
    res.coloc$summary,
    snps.interest
  )
  write.table(res.out,
              paste("output/drug.coloc", embedding, disease, chr.s, pos.s, pos.e, "txt", sep = "."),
              sep = "\t", row.names = F)
  
  ## locus-compare plot, again only where there is evidence
  if(do.plot & res.coloc$summary$PP.H4.abf > h4.ld.thr){
    source("../functions/plot_locus_compare.R")
    png(paste("graphics_coloc/coloc", embedding, disease, chr.s, pos.s, pos.e, "png", sep = "."),
        width = 16, height = 8, units = "cm", res = 300)
    par(mar = c(1.5, 1.5, 1, .5), mgp = c(.6, 0, 0), cex.axis = .5,
        cex.lab = .5, tck = .01, cex.main = .6, font.main = 2)
    layout(matrix(c(1, 1, 1, 2, 3, 4), 3, 2), heights = c(.43, .37, .2))
    plot.locus.compare(res.comb, res.out, res.out[ type.snp == "embQTL"]$ID, ld_mat)
    dev.off()
  }
  
  return(invisible(res.out))
}