#!/usr/bin/env Rscript
## script to compute partial correlation network among lifetime ICD-10 codes
## Maik Pietzner 16/06/2026  (optimised)
##
## method: each disease pair is conditioned on age, sex, AND all other diseases
## via the precision matrix (inverse of the correlation matrix). near-duplicate
## diseases are pruned first (|Pearson| > 0.7) so the matrix is invertible.
##
## called from 04_submit_parcor.sh; arguments are taken column-wise from
## input/Input.partial.correlation.txt, selected by the SLURM array index:
##   args[1] = sex.c   : "all" | "Male" | "Female"
##   args[2] = n.row   : optional, subsample this many people  (test mode)
##   args[3] = n.col   : optional, subsample this many diseases (test mode)
## a real run only carries column 1; a test run carries all three.
rm(list = ls())

## get the arguments from the command line
args <- commandArgs(trailingOnly = T)
## little options
options(stringsAsFactors = F)
## avoid conversion of numbers
options(scipen = 1)

## correct directory
setwd("<path_to_file>")

## packages needed
require(data.table)
require(Rfast)

#----------------------------#
##-- run / test controls  --##
#----------------------------#

## import sex in which to compute partial correlations
sex.c      <- args[1]
## optional test sizes (NA when the column is absent or empty -> full run)
n.row.test <- suppressWarnings(as.integer(args[2]))
n.col.test <- suppressWarnings(as.integer(args[3]))
test.mode  <- !is.na(n.row.test)
## cap I/O in test mode so reading the big files is fast, too
read.cap   <- if (test.mode) max(n.row.test * 4L, 20000L) else Inf
## minimum cases for a disease to enter the network (rare codes are unstable)
min.cases  <- 100L
set.seed(42)
cat("run partial correlations among", sex.c,
    if (test.mode) sprintf("[TEST: <=%s rows, <=%s cols]", n.row.test, n.col.test) else "",
    "\n")

#----------------------------#
##-- import relevant data --##
#----------------------------#

## import first occurrence: filtered by Wenhuan
ukb.first.occurrence <- fread("<path_to_file>",
                              nrows = read.cap)
## import labels
lab.first.occurrence <- fread("<path_to_file>")
## subset to what is included in data
lab.first.occurrence <- lab.first.occurrence[ id.ukbb %in% names(ukb.first.occurrence)]

## the dictionary ships its OWN icd10.code column and it is unreliable
## (e.g. "Date A00 first reported (cholera)" is tagged A27) -> drop it and
## re-derive the code ourselves from the description.
if ("icd10.code" %in% names(lab.first.occurrence))
  lab.first.occurrence[, icd10.code := NULL]

## flag the two field flavours by their exact description structure
lab.first.occurrence[, date.column   := grepl("^Date .* first reported", description)]
lab.first.occurrence[, source.column := grepl("^Source of report of",    description)]

## extract the ICD-10 code anchored BETWEEN "Date " and " first reported",
## so it can only ever capture the real code and only on genuine
## first-occurrence date fields (source / other date fields stay NA).
lab.first.occurrence[, icd10.code := {
  m <- regexpr("[A-Z][0-9]{2}", description)
  fifelse(m == -1L, NA_character_, substr(description, m, m + 2L))
}]

## phenotype label table = variables that genuinely carry an ICD-10 code
lab.phe              <- lab.first.occurrence[ date.column == TRUE & !is.na(icd10.code)]
cat("retained", nrow(lab.phe), "first-occurrence ICD-10 date fields\n")
## drop 'source' columns
ukb.first.occurrence <- ukb.first.occurrence[, c("f.eid", lab.phe$id.ukbb), with = F]
## recode to binary
ukb.first.occurrence[, (lab.phe$id.ukbb) := lapply(.SD, \(x) as.integer(!is.na(x))),
                     .SDcols = lab.phe$id.ukbb]
## import age and sex, and align numbers with the main experiment
ukb.dat              <- fread("../01_phenotype_preparation//data/UKB.data.prep.Delphi.embeddings.20260129.txt",
                              select = c("f.eid", "sex", "age"), nrows = read.cap)
## combine the two
ukb.dat              <- merge(ukb.dat, ukb.first.occurrence)
## delete what is no longer needed
rm(ukb.first.occurrence); gc(reset = T)

#----------------------------#
##-- apply sex-stratification#
#----------------------------#

## which cohort to run partial correlation on
if (sex.c != "all") {
  ## reduce to one sex
  ukb.dat <- ukb.dat[ sex == sex.c]
  ## covariate adjustment: age only
  cov.form <- ~ age
} else {
  ## recode for partial correlation network
  ukb.dat[, sex := ifelse(sex == "Male", 0, 1)]
  ## covariate adjustment: age + sex
  cov.form <- ~ age + sex
}

## drop rows with missing covariates (events are already 0/1, no NA)
ukb.dat              <- ukb.dat[ complete.cases(ukb.dat[, all.vars(cov.form), with = F])]

## test mode: subsample people
if (test.mode && nrow(ukb.dat) > n.row.test) {
  ukb.dat <- ukb.dat[sample(.N, n.row.test)]
}

#----------------------------#
##-- build matrices        --##
#----------------------------#

## binary event matrix
X                    <- as.matrix(ukb.dat[, lab.phe$id.ukbb, with = F])
storage.mode(X)      <- "double"                       # Rfast/Armadillo needs double, not integer

## test mode: subsample diseases
if (!is.na(n.col.test) && ncol(X) > n.col.test) {
  sel     <- sample(ncol(X), n.col.test)
  X       <- X[, sel, drop = FALSE]
  lab.phe <- lab.phe[sel]
}

## drop constant columns (e.g. sex-specific codes that are all-0 in a stratum)
cs                   <- colSums(X)
keep                 <- cs > 0 & cs < nrow(X)
if (any(!keep)) {
  cat("dropping", sum(!keep), "constant columns (no variance in this stratum)\n")
  X       <- X[, keep, drop = FALSE]
  lab.phe <- lab.phe[keep]
}
## case counts per phenotype (used for pruning)
lab.phe[, cases.w.primary := colSums(X)]

## drop rare diseases: filter BOTH lab.phe and X so it carries forward
keep.cc              <- lab.phe$cases.w.primary >= min.cases
if (any(!keep.cc)) {
  cat("dropping", sum(!keep.cc), "diseases with <", min.cases, "cases\n")
  lab.phe <- lab.phe[keep.cc]
  X       <- X[, lab.phe$id.ukbb, drop = FALSE]
}
cat("retained", nrow(lab.phe), "diseases for the network\n")

#----------------------------#
##-- unadjusted Pearson    --##
#----------------------------#

## compute correlation matrix as backbone
phe.cor              <- Rfast::cora(X)
## get names
dimnames(phe.cor)    <- list(lab.phe$id.ukbb, lab.phe$id.ukbb)

## long form (upper triangle only)
ut                   <- which(upper.tri(phe.cor), arr.ind = TRUE)
cor.dt               <- data.table(Var1    = lab.phe$id.ukbb[ut[, 1]],
                                   Var2    = lab.phe$id.ukbb[ut[, 2]],
                                   pearson = phe.cor[ut])

## identify highly correlated pairs and drop the member with FEWER cases
hi                   <- cor.dt[abs(pearson) > 0.7]
hi[lab.phe, n1 := i.cases.w.primary, on = c("Var1" = "id.ukbb")]
hi[lab.phe, n2 := i.cases.w.primary, on = c("Var2" = "id.ukbb")]
hi[, drop := fifelse(n1 >= n2, Var2, Var1)]
## reduce label set
lab.phe              <- lab.phe[!id.ukbb %in% unique(hi$drop)]
X                    <- X[, lab.phe$id.ukbb, drop = FALSE]

#----------------------------#
##-- partial correlation   --##
#----------------------------#

## conditioning covariates (age [+ sex]); raw columns, correlation centres them
cov.cols             <- all.vars(cov.form)             # "age" or c("age","sex")
Zc                   <- as.matrix(ukb.dat[, ..cov.cols])
storage.mode(Zc)     <- "double"

## full matrix to invert: [covariates, pruned diseases]
M                    <- cbind(Zc, X)
cn                   <- c(cov.cols, lab.phe$id.ukbb)
colnames(M)          <- cn

## correlation matrix -> precision matrix (inverse of correlation)
C                    <- Rfast::cora(M)
P                    <- tryCatch(solve(C),
                                 error = function(e) {
                                   cat("solve() singular; adding ridge 1e-4 to diagonal\n")
                                   solve(C + diag(1e-4, ncol(C)))
                                 })

## partial correlation = -P_ij / sqrt(P_ii * P_jj), conditioning on ALL others
dd                   <- sqrt(diag(P))
pc                   <- -P / outer(dd, dd)
diag(pc)             <- 1
pc                   <- pmax(pmin(pc, 1), -1)          # clamp float overshoot to [-1, 1]
dimnames(pc)         <- list(cn, cn)

## keep the disease-by-disease block only (drop age/sex rows & columns)
pe                   <- pc[lab.phe$id.ukbb, lab.phe$id.ukbb]

## vectorised p-values: df = n - (#variables conditioned on) - 2 = n - ncol(M)
df                   <- nrow(M) - ncol(M)
tst                  <- pe * sqrt(df / (1 - pe^2))
pv                   <- 2 * pt(-abs(tst), df)
pv[!is.finite(pv)]   <- NA_real_                       # guard r = 1 / numeric edge cases

#----------------------------#
##-- assemble + write      --##
#----------------------------#
ut                   <- which(upper.tri(pe), arr.ind = TRUE)
out                  <- data.table(Var1     = lab.phe$id.ukbb[ut[, 1]],
                                   Var2     = lab.phe$id.ukbb[ut[, 2]],
                                   estimate = pe[ut],
                                   pval     = pv[ut])
## add unadjusted Pearson r
out[cor.dt, pearson := i.pearson, on = c("Var1", "Var2")]
## add ICD-10 labels for both members
out[lab.phe, icd10.Var1 := i.icd10.code, on = c("Var1" = "id.ukbb")]
out[lab.phe, icd10.Var2 := i.icd10.code, on = c("Var2" = "id.ukbb")]

## write results to file (separate name in test mode, never overwrites a real run)
fout                 <- paste("output/partial.correlation", sex.c,
                              if (test.mode) "TEST" else NULL, "txt", sep = ".")
fwrite(out, fout, sep = "\t", na = "NA", row.names = F)

cat("done:", nrow(out), "pairs across", nrow(lab.phe), "phenotypes ->", fout, "\n")