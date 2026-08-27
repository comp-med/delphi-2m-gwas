###############################################################
#### reconstruction of exposures from the Delphi embeddings ####
#### within sex (+ pooled), run via submit_reconstruction_r2 ####
#### Maik Pietzner                              17/06/2026  ####
###############################################################

rm(list=ls())
setwd("<path_to_file>")
options(stringsAsFactors = F)

## --> packages needed <-- ##
require(data.table)
require(doMC)

## cores allocated by SLURM (falls back to 1 outside the scheduler)
n.cores   <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
cat("Using", n.cores, "cores\n")

## --> person-level sex variable in ukb.dat -- ADJUST to your column / coding <-- ##
## expected coding after recode below: 1 = male, 0 = female (UKB f.31.0.0 is 1/0)
sex.col   <- "sex"

###############################################
####                import                 ####
###############################################

## per-person mean-pooled embeddings: f.eid + emb1..emb120 (no missingness)
emb       <- fread("<path_to_file>")
## exposures + labels (lab.set carries the analysis sex: all / men / women)
ukb.dat   <- fread("input/UKB.data.prep.ICD10.code.Cox.models.20260626.txt")
lab.set   <- fread("input/UKB.labels.prep.ICD10.code.Cox.models.20260626.txt")

## embedding columns
emb.cols  <- grep("^emb", names(emb), value = T)
## n = 120
length(emb.cols)

## reconstruct continuous + binary risk factors only (skip the disease/drug
## tokens and multi-level categoricals; the latter would need a redundancy R2)
lab.use   <- lab.set[ !(category %in% c("Diseases", "Drugs")) &
                        type %in% c("numeric", "integer", "logical")]
expo.set  <- intersect(lab.use$short_name, names(ukb.dat))
## n exposures
length(expo.set)

## rank-based inverse-normal transform (as in run.cox; applied WITHIN stratum)
INT       <- function(x) qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))

## one merge; carry sex; embeddings as a numeric matrix built once
dat       <- merge(ukb.dat[, c("f.eid", sex.col, expo.set), with = F], emb, by = "f.eid")
## recode sex to 1 = male / 0 = female (handles UKB 1/0 or "Male"/"Female")
dat[, male := as.integer(get(sex.col) %in% c(1, "1", "Male", "male", "M"))]
X         <- as.matrix(dat[, ..emb.cols])
ok.x      <- complete.cases(X)

###############################################
####     per-stratum reconstruction task    ####
###############################################

## 'all'        -> reconstruct pooled AND within men AND within women
## 'men'/'women'-> only that stratum (sex-specific exposures)
tasks     <- rbindlist(lapply(expo.set, function(x){
  sx <- lab.set[short_name == x, sex][1]
  st <- if(is.na(sx) | sx == "Both") c("all", "men", "women") else sx
  data.table(short_name = x, sex = st)
}))
## n tasks
nrow(tasks)

###############################################
####          reconstruction R^2           ####
###############################################

## each fork copies the stratum block of X -- lower --cpus-per-task if memory-bound
registerDoMC(n.cores)
res.recon <- rbindlist(mclapply(1:nrow(tasks), function(j){
  
  x   <- tasks$short_name[j]
  st  <- tasks$sex[j]
  
  ## stratum mask (within sex removes the constant sex dimension of the embeddings)
  sm  <- switch(st, all = rep(T, nrow(dat)), men = dat$male == 1L, women = dat$male == 0L)
  
  ## response within the stratum: 0/1 if binary, rank-INT if continuous
  ys  <- dat[[x]][sm]
  u   <- unique(na.omit(ys))
  if(length(u) < 2) return(NULL)                         ## no variation in stratum
  bin <- is.logical(dat[[x]]) | length(u) == 2
  y   <- if(bin) as.integer(as.factor(ys)) - 1L else INT(as.numeric(ys))
  
  ## complete cases within the stratum (aligned to the sm subset)
  ok  <- !is.na(y) & ok.x[sm]
  if(sum(ok) < 1000) return(NULL)
  
  ## linear reconstruction from the embeddings (.lm.fit is light-weight)
  xs  <- X[sm, , drop = F]
  fit <- .lm.fit(cbind(1, xs[ok, , drop = F]), y[ok])
  rss <- sum(fit$residuals^2)
  tss <- sum((y[ok] - mean(y[ok]))^2)
  
  data.table(short_name = x,
             sex        = st,
             type       = if(bin) "binary" else "continuous",
             n          = sum(ok),
             recon.r2   = 1 - rss / tss)
  
}, mc.cores = n.cores), fill = T)

## attach labels + category; tidy ordering (exposure, then all/men/women)
res.recon <- merge(res.recon, lab.set[, .(short_name, label, category)],
                   by = "short_name", all.x = T)
res.recon[, sex := factor(sex, levels = c("all", "men", "women"))]
res.recon <- res.recon[ order(short_name, sex)]

## write out, datestamped -- this aggregated table is safe to share
fwrite(res.recon, paste0("output/Reconstruction_R2_exposures_by_sex_Delphi_embeddings.",
                         format(Sys.Date(), "%Y%m%d"), ".txt"),
       sep = "\t", row.names = F, quote = F, na = NA)

cat("Done -", nrow(res.recon), "exposure x stratum reconstructions\n")