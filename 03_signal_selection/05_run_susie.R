#sessionInfo()
#Try to close a graphics device if possible, and remove all from the env
rm(list = ls())
try(dev.off())

# Load required packages
library(dplyr)
library(magrittr)
library(reshape2)
library(susieR)
library(data.table)
library(corrcoverage)

set.seed(42)

setwd(
  '<path_to_file>'
)

# Get the LD function
source('../functions/get_dosages.R')
source('../functions/function_wakefield.R')


# Output of the final table:
# rsID  CHROM GENPOS  ALLELE0 ALLELE1 pval_marginal beta_marginal se_marginal pval_joint  beta_joint  se_joint  R2_leadvariant  cs  pip method  pheno startpos_region endpos_region index region

###### Command line arguments
args <- commandArgs(trailingOnly = T)
index <- args[1]
pheno <- args[2] # this is your embedding: emb103, emb100, ...
region <- args[3]
chr.s <- args[4]
pos.s <- args[5]
pos.e <- args[6]
#
# # Test cases
# index <- '14900'
# pheno <- 'Pyruvate'
# region <- '10_1'
# chr.s <- 10
# pos.s <- 2608283
# pos.e <- 3527533

finemap_dir <- file.path("gwas/output/finemapping", pheno, index)


print(paste0(
  'Variables: index: ',
  index,
  ' | pheno: ',
  pheno,
  ' | region: ',
  region,
  ' | chromosome: ',
  chr.s,
  ' | start: ',
  pos.s,
  ' | end: ',
  pos.e
))

# Loading the association stats from REGENIE
## import depending on phenotype
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

## read in LD matrix (embedding-specific)
ld_path <- file.path(
  "gwas/output/proc_emb120_region_data",
  pheno, # e.g. emb103
  region, # e.g. 19_2
  "ldmat.ld"
)

if (!file.exists(ld_path)) {
  stop(paste0("LD matrix not found at: ", ld_path))
}

ld <- fread(ld_path, data.table = FALSE) %>% as.matrix()
ld[is.na(ld)] <- 0
diag(ld) <- 1
ld <- data.frame(ld)
colnames(ld) <- 1:ncol(ld)

## create identifier column in results to keep the mapping to the LD matrix
res$id <- 1:nrow(res)

#-----------------------------------------#
##--          run fine-mapping         --##
#-----------------------------------------#

# Run fine-mapping to obtain 95%-credible sets.
# This throws an errors due to rounding errors in the LD matrix

## Iterate over different values of L in SuSie
## Don't try to do this in parallel since that leads to very big memory usage
res.fine <- lapply(2:10, function(x) {
  print(paste0('L: ', x))
  tmp <- tryCatch(
    {
      susie_rss(
        z = res$BETA / res$SE,
        R = as.matrix(ld),
        n = 434646,
        L = x,
        coverage = .95,
        min_abs_corr = .1,
        max_iter = 10000
      )
    },
    error = function(e) {
      return(list(
        pip = rep(NA, nrow(res)),
        sets = list(cs = NA),
        converged = F
      ))
    }
  )
  ## add L
  tmp$L <- x
  ## return
  return(tmp)
})

## reduce to entries with at least one credible set
jj <- unlist(lapply(res.fine, function(x) length(na.omit(x$sets$cs))))
## delete
res.fine <- res.fine[which(jj > 0)]
print(paste0(
  'We have ',
  sum(jj > 0),
  '/9 iterations with at least one credible set'
))

# If we have results left (that is, at least one of our iterations found at least one credible set)
if (length(res.fine) > 0) {
  #------------------------------------------#
  ##--        prune sets if needed        --##
  #------------------------------------------#
  print('Pruning the sets.')

  ## do LD assessment
  ld.top <- lapply(res.fine, function(x) {
    ## add PIPs to ease selection of SNPs
    tmp <- summary(x)$vars
    ## subset to those in credible sets
    tmp <- as.data.table(subset(tmp, cs > 0))
    ## only top SNPs
    tmp <- tmp[order(cs, -variable_prob)]
    ## get only the top
    tmp[, ind := 1:.N, by = "cs"]
    tmp <- tmp[ind == 1]
    ## get the names
    tmp[, id := sapply(tmp$variable, function(x) rownames(ld)[x])]
    ## generate LD
    top.ld <- ld[tmp$id, tmp$id, drop = F]^2
    ## identify possible sets in LD (r2>0.25); set diagonal to zero to ease downstream analysis
    diag(top.ld) <- 0
    top.ld <- reshape2::melt(top.ld)
    ## subset to possible problematic candidates
    top.ld <- subset(top.ld, value >= .25)
    ## return
    return(top.ld)
  })

  ## find the maximum set of unrelated variants
  jj <- unlist(lapply(ld.top, nrow))

  ## subset accordingly
  if (sum(jj > 0) > 0) {
    res.fine <- res.fine[-which(jj > 0)]
  }

  print(paste0(
    'After LD pruning, we are left with: ',
    length(res.fine),
    ' results'
  ))

  ## add to the results, if any
  if (length(res.fine) > 0) {
    print('Continuing since there are results left after pruning')

    ## If there are more than one results, take the one with the highest
    ## number of credible sets to not miss any secondary signals.
    res.fine <- res.fine[[length(res.fine)]]

    ## add the information on pip and cs to the summary statistics
    tmp <- summary(res.fine)$vars
    colnames(tmp) <- c("variable", "pip", "cs")

    # All credible sets in here at the moment
    # Later on, ones that dont pass additional quality checks are set to -1
    tmp <- merge(res, tmp, by.x = "id", by.y = 'variable') #%>% filter(cs > 0)

    #------------------------------------------------#
    ##-- run joint model to ensure GWAS threshold --##
    #------------------------------------------------#

    # We have candidate SNPs that are in the credible set, run joint models to ensure GWAS threshold
    # We need only the top SNPs per cs

    ## identify top SNP for each credible set
    ## with_ties is important here as there are edge cases where variants are in perfect LD with each other, causing
    ## dplyr to return more than one variant per cs leading to collinearity in the glm below.
    top.snp <- tmp %>%
      filter(cs > 0) %>%
      group_by(cs) %>%
      slice_max(pip, with_ties = F) %>%
      pull(id)

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

    # Write a list of snp ids for which we need the dosages
    cat(
      tmp %>%
        filter(cs > 0) %>%
        group_by(cs) %>%
        slice_max(pip, with_ties = F) %>%
        pull(ID),
      #file = paste0('gwas/output/finemapping/', pheno, '/', index, '/variants_dosages.txt'),
      file = file.path(finemap_dir, "variants_dosages.txt"),
      sep = '\n'
    )

    foo <- get_dosages(index, chr.s, pos.s, pos.e, pheno)
    ## separate out into two data sets to ease downstream computation
    snp.dat <- foo[[1]]
    snp.info <- foo[[2]]

    ## delete and clear up
    rm(foo)
    gc()

    ## add results variant ID to snp.info as well as credible set information
    snp.info <- merge(snp.info, tmp, by = 'MarkerName')

    ## add SNP data
    phen <- merge(phen, snp.dat, by.x = "eid", by.y = "ID_1")

    # Running a joint model on the top SNP per credible set
    snp.info %>% filter(id %in% top.snp) %>% pull(XID) -> vars

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

    ## change one name
    colnames(m.joint) <- gsub(
      "Pr\\(>\\|t\\|\\)",
      "PVALUE.joint",
      colnames(m.joint)
    )

    ## add LOG10P from the GWAS
    m.joint$XID <- row.names(m.joint)
    m.joint <- merge(m.joint, snp.info, by = 'XID')

    # m.joint             <- merge(m.joint, res %>% select(-pip, -cs, -MarkerName, -EXTRA), by='id')

    # Implement the following filters:
    # Variants need genome-wide significance in joint model and marginal model both
    # Beta needs same sign
    # Beta can change no more than 25%. Very arbitrary

    m.joint %<>%
      mutate(
        sig = ifelse(PVALUE.joint < 5e-8 & LOG10P > 7.3, T, F),
        dir_concordance = ifelse(sign(Estimate) == sign(BETA), T, F),
        lim = ifelse(
          (abs(Estimate) < abs(BETA) + (0.25 * abs(BETA))) &
            (abs(Estimate) > abs(BETA) - (0.25 * abs(BETA))),
          T,
          F
        )
      )

    print('These are the aggregated statistics for GWAS and joint model:')
    print(m.joint)

    m.joint$keep <- ifelse(
      m.joint$sig & m.joint$dir_concordance & m.joint$lim,
      T,
      F
    )

    print('Keeping only the variants that pass all filters: ')
    print(table(m.joint$keep))

    # Now only contains the credible sets we want to keep
    m.joint %<>% filter(keep)

    # Set the other ones in tmp to -1 if we don't want to keep these
    # Don't throw them out since we still output the pips
    tmp %<>% mutate(cs = ifelse(cs %in% m.joint$cs, cs, -1))

    ## Filter the regional association statistics for the credible sets to keep
    # res %<>% filter(cs %in% cs.keep)

    # If any results left, ..
    if (nrow(m.joint) > 0) {
      ## 1. We run another linear model with the lead variants that remain after filtering to check genome-wide threshold

      # 1. Running a joint model on the top SNPs per credible set
      unique(m.joint$XID) -> vars

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
      # Double check if we have all the variants
      stopifnot(nrow(m.joint) == length(vars))

      ## change one name
      colnames(m.joint) <- gsub(
        "Pr\\(>\\|t\\|\\)",
        "PVALUE.joint",
        colnames(m.joint)
      )

      ## add LOG10P from the GWAS
      m.joint$XID <- row.names(m.joint)
      m.joint <- merge(
        m.joint,
        snp.info %>% select(all_of(c('MarkerName', 'XID'))),
        by = 'XID'
      )

      colnames(m.joint) <- c(
        'XID',
        'estimate_joint',
        'se_joint',
        'tstat_joint',
        'pval_joint',
        'MarkerName'
      )

      # res now contains ALL variants contains marginal as well as joint statistics, credible set (-1 if not in CS) and pip
      res <- merge(tmp, m.joint, by = 'MarkerName', all.x = T)

      # 2. Add the LD to the top variant
      res$R2_leadvariant <- NA
      for (credset in unique(res$cs)) {
        if (credset == -1) {
          next
        }
        res %>%
          filter(cs == credset) %>%
          slice_max(order_by = pip, with_ties = F) %>%
          pull(id) -> top
        res %>% filter(cs == credset) %>% pull(id) -> other

        # Mutate every variant
        for (o in other) {
          res %<>%
            mutate(
              R2_leadvariant = unlist(ifelse(
                id == o,
                ld[o, top]^2,
                R2_leadvariant
              ))
            )
        }
      }

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
        R2_leadvariant = res$R2_leadvariant,
        cs = res$cs,
        pip = res$pip,
        method = 'SuSie',

        # General info on region
        pheno = pheno,
        startpos_region = pos.s,
        endpos_region = pos.e,
        index = index,
        region = region
      )

      # Export tables
      fwrite(
        out,
        file = file.path(finemap_dir, "finemapped_results.txt.gz"),
        sep = '\t',
        row.names = F,
        quote = F
      )
    } else {
      # Calculate Wakefield approximation
      out <- doWakefield(
        pheno = pheno,
        chr.s = chr.s,
        pos.s = pos.s,
        pos.e = pos.e,
        region = region
      )
      fwrite(
        out,
        file = file.path(finemap_dir, "finemapped_results.txt.gz"),
        sep = '\t',
        row.names = F,
        quote = F
      )
    }
  } else {
    # Calculate Wakefield approximation
    out <- doWakefield(
      pheno = pheno,
      chr.s = chr.s,
      pos.s = pos.s,
      pos.e = pos.e,
      region = region
    )
    fwrite(
      out,
      file = file.path(finemap_dir, "finemapped_results.txt.gz"),
      sep = '\t',
      row.names = F,
      quote = F
    )
  }
} else {
  # Calculate Wakefield approximation
  out <- doWakefield(
    pheno = pheno,
    chr.s = chr.s,
    pos.s = pos.s,
    pos.e = pos.e,
    region = region
  )
  fwrite(
    out,
    file = file.path(finemap_dir, "finemapped_results.txt.gz"),
    sep = '\t',
    row.names = F,
    quote = F
  )
}

# Clean up some files that we no longer need
# system(paste0('rm output/finemapping/', index, '/*.bgen'))
# system(paste0('rm output/finemapping/', index, '/*dosage*'))
