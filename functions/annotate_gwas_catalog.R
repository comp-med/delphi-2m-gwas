###################################################
## function to annotate a set of variants with 
## results from the GWAS catalog

annotate.variants.gwas <- function(gen.var, gwas.catalogue, r2.proxies, r2.thr){
  
  ## 'gen.var'        -- data set for genetic variants to be annotated
  ## 'gwas.catalogue' -- parsed release of the GWAS catalogue (assumed in build38)
  ## 'r2.proxies'     -- lsit of proxy variants for the ones in 'gen.var'
  ## 'r2.thr'         -- r2 threshold to declare linkage for proxa variants
  
  ## do in parallel
  registerDoMC(10)
  
  ## go through each variant and 1) identify all proxies, 2) map to GWAS catalog findings, and 3) reduce redundancy
  res.gwas <- mclapply(1:nrow(gen.var), function(x){
    
    print(x)
    
    ## get all possible proxies (careful; does not include the SNP itself)
    snp               <- r2.proxies[ (ID_A == gen.var$id[x] | ID_B == gen.var$id[x]) & UNPHASED_R2 >= r2.thr]
    
    ## convert back to easy format
    snp[, lead.id := gen.var$id[x]]
    snp[, proxy.id := ifelse(ID_A == gen.var$id[x], ID_B, ID_A)]
    snp[, pos.proxy := ifelse(ID_A == gen.var$id[x], pos_hg38_B, pos_hg38_A)]
    ## subset to minimum needed
    snp              <- snp[, .(lead.id, proxy.id, CHROM_A, pos.proxy, UNPHASED_R2)]
    
    ## create snp id to optimize merging multiple mappings
    snp[, snp.id := paste0(CHROM_A, ":", pos.proxy)]
    
    ## create two versions of mapping
    snp.id          <- merge(snp, gwas.catalogue, by.x="proxy.id", by.y="SNPS")
    snp.pos           <- merge(snp, gwas.catalogue, by = "snp.id")
    ## edit
    snp.id$snp.id   <- snp.id$snp.id.x
    snp.id$snp.id.x <- snp.id$snp.id.y <- NULL
    snp.pos$SNPS      <- snp.pos$proxy.id
    snp.id$SNPS     <- snp.id$proxy.id
    ## combine
    snp               <- unique(rbind(snp.id, snp.pos))
    
    ## prepare return
    if(nrow(snp) > 0){
      ## sort 
      snp              <- snp[order(TRAIT, lead.id, -UNPHASED_R2)]
      ## create indicator
      snp[, ind := 1:.N, by=c("TRAIT", "lead.id")]
      ## keep only one finding per trait
      snp              <- snp[ind == 1]
      
      ## report summary back
      snp              <- data.table(id.gwas=paste(sort(unique(snp$proxy.id)), collapse = "||"),
                                     trait_reported=paste(sort(unique(snp$TRAIT)), collapse = "||"),
                                     mapped_trait=paste(sort(unique(snp$MAPPED_TRAIT)), collapse = "||"),
                                     mapped_trait_efo=paste(sort(unique(snp$MAPPED_TRAIT_URI)), collapse = "||"),
                                     study_id=paste(sort(unique(snp$`STUDY ACCESSION`)), collapse = "||"),
                                     source_gwas=paste(sort(unique(snp$PUBMEDID)), collapse = "||"),
                                     num_reported=nrow(snp))
      
    }else{
      
      ## report summary back
      snp              <- data.table(id.gwas="",
                                     trait_reported="",
                                     mapped_trait="",
                                     mapped_trait_efo="",
                                     study_id="",
                                     source_gwas="",
                                     num_reported=0)
    }
    
    ## return data set
    return(data.table(gen.var[x,], snp))
  }, mc.cores=10)
  # })
  ## combine 
  res.gwas <- rbindlist(res.gwas)
  
  return(res.gwas)
  
}