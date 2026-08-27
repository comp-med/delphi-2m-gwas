doWakefield <- function(pheno, chr.s, pos.s, pos.e, region = region) {
  print('SuSie failed, using Wakefield approximation')

  res <- fread(
    cmd = paste0(
      "zcat <path_to_file>",
      pheno,
      "/gwas_emb120_chr",
      chr.s,
      "_",
      pheno,
      ".regenie.gz | awk -v chr=",
      chr.s,
      " -v low=",
      pos.s,
      " -v upp=",
      pos.e,
      " '{if(($1 == chr && $2 >= low && $2 <= upp) || NR == 1) print $0}'"
    ),
    sep = " ",
    header = T,
    data.table = F
  )

  ## drop SNPs that have possibly failed in REGENIE
  res <- subset(res, is.na(EXTRA))

  ## create MarkerName to enable mapping to the LD file
  res$MarkerName <- apply(res, 1, function(x) {
    paste0(
      "chr",
      as.numeric(x[1]),
      ":",
      as.numeric(x[2]),
      "_",
      paste(sort(x[4:5]), collapse = "_")
    )
  })
  res$id <- 1:nrow(res)

  ## read in LD matrix
  ld <- fread(
    paste0(
      'gwas/output/proc_emb120_region_data/',
      pheno,
      '/',
      region,
      '/ldmat.ld'
    ),
    data.table = F
  ) %>%
    as.matrix()
  ld[is.na(ld)] <- 0
  diag(ld) <- 1
  ld <- data.frame(ld)
  colnames(ld) = 1:ncol(ld)

  ## Implement Wakefield approximation for when Susie fails to converge
  z <- res$BETA / res$SE
  v <- res$SE^2
  res$pip <- ppfunc(z, v, W = .2)

  ## compute the 95%-credible set
  cred.set <- credset(res$pip, thr = .95)

  ## add to the data
  res$cs <- NA
  res$cs[cred.set$credset] <- 1

  ## add LD (this time only one causal variant)
  ii <- res[which.max(res$pip), 'id']

  ## add LD of the remaining variants in the credset to the lead variant
  ii <- data.frame(snp.id = 1:nrow(ld), R2.1 = ld[, ii]^2)
  res <- merge(res, ii, by.x = 'id', by.y = 'snp.id')
  res[which(is.na(res$cs)), 'R2.1'] <- NA # We don't need the LD to variants that are not inside credible sets

  # Take only the variant that are in an actual credible set
  # res %<>% filter(!is.na(cs))

  #####
  # Pseudo joint model
  #####

  # Check whether lead variant reaches genome-wide thresholds. Only one variant
  top.snp <- res[which.max(res$pip), 'id']

  # Import necessary data to run the models, filter on the right individuals
  phen <- fread(paste0(
    '<path_to_file>'
  )) %>%
    select(all_of(c('FID', pheno))) %>%
    set_colnames(c('eid', 'phenotype'))
  covs <- fread(
    '<path_to_file>'
  )
  inc.list <- fread(
    '<path_to_file>'
  )
  phen <- merge(phen, covs, by.x = 'eid', by.y = 'FID') %>%
    filter(eid %in% inc.list$V1)

  # Write a list of snp ids that we need the dosages for to run the model
  tmp <- res[which.max(res$pip), 'ID']
  cat(
    tmp,
    file = paste0(
      'gwas/output/finemapping/',
      pheno,
      '/',
      index,
      '/variants_dosages.txt'
    ),
    sep = '\n'
  )

  foo <- get_dosages(index, chr.s, pos.s, pos.e, pheno)
  ## separate out into two data sets to ease downstream computation
  snp.dat <- foo[[1]]
  snp.info <- foo[[2]]

  ## add results variant ID to snp.info as well as credible set information
  snp.info <- merge(
    snp.info,
    res[, c('MarkerName', 'pip', 'cs', 'id')],
    by = 'MarkerName'
  )

  ## add SNP data
  phen <- merge(phen, snp.dat, by.x = "eid", by.y = "ID_1")

  #--------------------------------#
  ##--      run joint model     --##
  #--------------------------------#
  print('Starting generalised linear model')

  # Running a joint model on the top SNPs per credible set
  snp.info %>% filter(id %in% top.snp) %>% pull(XID) -> vars

  ## run model
  m.joint <- summary(glm(
    paste(
      'phenotype',
      " ~",
      paste(c(vars, "age", "sex", paste0("pc", 1:10)), collapse = " + ")
    ),
    data = phen,
    family = gaussian
  ))$coefficients[-1, , drop = F]
  ## look at SNPs only
  m.joint <- as.data.frame(m.joint[vars, , drop = F])

  # Normal colnames
  colnames(m.joint) <- c(
    'estimate_joint',
    'se_joint',
    'tstat_joint',
    'pval_joint'
  )

  ## add LOG10P from the GWAS
  m.joint$XID <- row.names(m.joint)
  m.joint <- merge(m.joint, snp.info, by = 'XID')

  m.joint <- merge(
    m.joint,
    res %>% select(-pip, -cs, -MarkerName, -EXTRA),
    by = 'id'
  )

  print('These are the aggregated statistics for GWAS and joint model:')
  print(m.joint)

  ## keep only what passes the genome-wide threshold in the joint model
  print(
    'Keeping these lead variants that pass both joint and marginal statistics: '
  )

  m.joint %<>%
    mutate(
      sig = ifelse(pval_joint < 5e-8 & LOG10P > 7.3, T, F),
      dir_concordance = ifelse(sign(estimate_joint) == sign(BETA), T, F),
      lim = ifelse(
        (abs(estimate_joint) < abs(BETA) + (0.25 * abs(BETA))) &
          (abs(estimate_joint) > abs(BETA) - (0.25 * abs(BETA))),
        T,
        F
      )
    )

  m.joint$keep <- ifelse(
    m.joint$sig & m.joint$dir_concordance & m.joint$lim,
    T,
    F
  )

  print('Keeping only the variants that pass all filters: ')
  print(table(m.joint$keep))

  res <- merge(res, m.joint, all.x = T)

  # return(res)
  # rsID  CHROM GENPOS  ALLELE0 ALLELE1 pval_marginal beta_marginal se_marginal pval_joint  beta_joint  se_joint  R2_leadvariant  cs  pip method  pheno startpos_region endpos_region index region

  # We are now satisfied, create the final output
  out <- data.frame(
    # General info
    id = res$ID,
    chrom = res$CHROM,
    genpos = res$GENPOS,
    allele0 = res$ALLELE0,
    allele1 = res$ALLELE1,

    # Marginal and joint statistics
    pval_marginal = res$LOG10P,
    beta_marginal = res$BETA,
    se_marginal = res$SE,
    pval_joint = -log10(res$pval_joint),
    beta_joint = res$estimate_joint,
    se_joint = res$se_joint,

    # Fine mapping statistics
    R2_leadvariant = res$R2.1,
    cs = res$cs,
    pip = res$pip,
    method = 'Wakefield',

    # General info on region
    pheno = pheno,
    startpos_region = pos.s,
    endpos_region = pos.e,
    index = index,
    region = region
  )

  return(out)
}
