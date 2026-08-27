####################################################################
#### cross-fitted delta-AUC: exposure gain over the Delphi LP     ##
#### stratified (age x sex) where populated; global age/sex-       ##
#### adjusted out-of-sample fallback when strata collapse.         ##
#### Returns the component AUCs (auc.base.cv -> auc.comb.cv) so     ##
#### the increment can be plotted; `mode` records which path ran. ##
#### Maik Pietzner                                   19/06/2026   ##
####################################################################

## --> packages needed <-- ##
require(data.table)

#-------------------------------------------------#
##--   helper: rank-INT                          --##
#-------------------------------------------------#
INT          <- function(x) qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))

#-------------------------------------------------#
##--   helper: DeLong delta-AUC + paired SE      --##
#-------------------------------------------------#
## delta = AUC(m2) - AUC(m1) on the same samples (placement method); also returns auc1, auc2
dl.cv        <- function(label, m1, m2){
  label      <- as.integer(label); ca <- label == 1L; co <- label == 0L
  m          <- sum(ca);           n  <- sum(co)
  if(m < 2L | n < 2L) return(list(delta = NA_real_, se = NA_real_, auc1 = NA_real_, auc2 = NA_real_))
  comp       <- function(s){
    X <- s[ca]; Y <- s[co]; ys <- sort(Y); xs <- sort(X)
    nlt.y <- findInterval(X, ys, left.open = T); nle.y <- findInterval(X, ys)
    V10   <- (nlt.y + 0.5 * (nle.y - nlt.y)) / n
    nlt.x <- findInterval(Y, xs, left.open = T); nle.x <- findInterval(Y, xs)
    V01   <- ((m - nle.x) + 0.5 * (nle.x - nlt.x)) / m
    list(V10 = V10, V01 = V01, auc = mean(V10))
  }
  a <- comp(m1); b <- comp(m2)
  vd <- var(b$V10 - a$V10) / m + var(b$V01 - a$V01) / n
  list(delta = b$auc - a$auc, se = sqrt(vd), auc1 = a$auc, auc2 = b$auc)
}

#-------------------------------------------------#
##--   cross-fitted delta-AUC for one exposure   --##
#-------------------------------------------------#
## The combined score glm(lab ~ lp + expo) is fit on K-1 folds and scored on the held-out fold
## (out-of-sample throughout); the Delphi-only baseline is the raw LP (no free parameters, so its
## AUC needs no cross-fitting). The gain is summarised within age x sex strata, case-weighted.
## When fewer than `min.strata` strata clear `min.cell`, we fall back to a SINGLE global out-of-sample
## DeLong on ALL cases, adjusting for age (and sex) inside both models -- this keeps every case and
## stays adjusted and out-of-sample, at the cost of a wider (honest) interval and a different reference
## frame (the baseline becomes the age/sex model, not the within-stratum LP). `mode` records which ran.
## auc.base.cv -> auc.comb.cv equals delta.auc.cv exactly in both paths.
delphi.delta.auc.cv <- function(dat, exposure,
                                dlp.col = "delphi.lp", time.col = "t.follow",
                                event.col = "event.occurred", age.col = "age", sex.col = "sex",
                                inv.transform = T, horizon.years = 10, age.bin = 5,
                                n.folds = 5L, min.cell = 25L, min.strata = 2L, seed = 42L){
  
  hz         <- horizon.years * 365.25
  d          <- dat[, .(lp   = get(dlp.col),  agev = get(age.col), sexv = as.character(get(sex.col)),
                        tf   = get(time.col), ev   = get(event.col), x = get(exposure))]
  
  ## horizon binarisation (case: event <= horizon; control: followed >= horizon; else drop)
  d[, lab := NA_integer_]
  d[ev == 1 & tf <= hz, lab := 1L]
  d[tf >= hz,           lab := 0L]
  
  ## single-df continuous encoding
  d[, expo := if(inv.transform) INT(as.numeric(x)) else as.numeric(x)]
  d          <- d[!is.na(lab) & !is.na(lp) & !is.na(expo) & !is.na(agev)]
  if(d[lab == 1L, .N] < min.cell) return(NULL)
  
  d[, strat := paste0(floor(agev / age.bin) * age.bin, "_", sexv)]
  
  ## case-balanced folds
  set.seed(seed)
  d[, fold := sample(rep_len(1:n.folds, .N)), by = lab]
  
  ## sex term only if more than one level present (drops out for men/women runs)
  sx         <- if(uniqueN(d$sexv) > 1) " + sexv" else ""
  f.comb     <- stats::as.formula("lab ~ lp + expo")
  f.base.a   <- stats::as.formula(paste0("lab ~ lp + agev", sx))
  f.comb.a   <- stats::as.formula(paste0("lab ~ lp + expo + agev", sx))
  
  ## cross-fit: comb (for the stratified path, vs raw lp) and the age/sex-adjusted pair (for the fallback)
  d[, `:=`(comb = NA_real_, base.a = NA_real_, comb.a = NA_real_)]
  for(k in 1:n.folds){
    tr <- d[fold != k]; te <- d[fold == k]
    d[fold == k, comb   := tryCatch(predict(glm(f.comb,   tr, family = binomial()), te, type = "link"), error = function(e) NA_real_)]
    d[fold == k, base.a := tryCatch(predict(glm(f.base.a, tr, family = binomial()), te, type = "link"), error = function(e) NA_real_)]
    d[fold == k, comb.a := tryCatch(predict(glm(f.comb.a, tr, family = binomial()), te, type = "link"), error = function(e) NA_real_)]
  }
  
  ## stratified DeLong (baseline = raw Delphi LP, combined = OOS lp + expo), case-weighted
  rs         <- rbindlist(lapply(unique(d$strat), function(s){
    ds <- d[strat == s & !is.na(comb)]
    if(ds[lab == 1L, .N] < min.cell | ds[lab == 0L, .N] < min.cell) return(NULL)
    dd <- dl.cv(ds$lab, ds$lp, ds$comb)
    data.table(n.case = ds[lab == 1L, .N], est = dd$delta, se = dd$se, a1 = dd$auc1, a2 = dd$auc2)
  }), fill = T)
  
  if(nrow(rs) >= min.strata){
    w        <- rs$n.case
    est      <- sum(w * rs$est) / sum(w)
    se       <- sqrt(sum((w / sum(w))^2 * rs$se^2))
    auc.base <- weighted.mean(rs$a1, w); auc.comb <- weighted.mean(rs$a2, w)
    mode     <- "stratified"; nstr <- nrow(rs); ncase <- d[lab == 1L, .N]
  } else {
    ## global age/sex-adjusted out-of-sample DeLong on ALL cases
    dg       <- d[!is.na(base.a) & !is.na(comb.a)]
    dd       <- dl.cv(dg$lab, dg$base.a, dg$comb.a)
    est      <- dd$delta; se <- dd$se; auc.base <- dd$auc1; auc.comb <- dd$auc2
    mode     <- "global_adjusted"; nstr <- 1L; ncase <- dg[lab == 1L, .N]
  }
  
  if(is.na(est) | is.na(se)) return(NULL)
  z          <- est / se
  data.table(exposure     = exposure,
             delta.auc.cv = est, se.cv = se,
             lci.cv       = est - 1.96 * se, uci.cv = est + 1.96 * se,
             pval.cv      = 2 * pnorm(-abs(z)),
             auc.base.cv  = auc.base, auc.comb.cv = auc.comb,
             n.case       = ncase, n.strata.cv = nstr, mode = mode)
}

#-------------------------------------------------#
##--   biomarker gain OVER DEMOGRAPHICS (g_B)    --##
#-------------------------------------------------#
## Same machinery, but the baseline is age (+ sex) instead of the Delphi LP. gain.cv =
## AUC(demographics + exposure) - AUC(demographics): the exposure's own discrimination gain over
## demographics. The exposure score is cross-fitted (sign learned out-of-fold) so PROTECTIVE
## exposures are not scored as negative discrimination -- gain.cv is direction-agnostic and >= 0 in
## expectation (the sign of the association lives in the Cox 'direction'/coef columns, not here).
## Within an age x sex cell demographics are ~constant (baseline AUC 0.5); the global fallback
## adjusts age/sex inside both models. `gain.mode` records which path ran; auc.dem.cv -> auc.demB.cv
## equals gain.cv exactly in both paths.
biomarker.gain.cv <- function(dat, exposure,
                              time.col = "t.follow", event.col = "event.occurred",
                              age.col = "age", sex.col = "sex",
                              inv.transform = T, horizon.years = 10, age.bin = 5,
                              n.folds = 5L, min.cell = 25L, min.strata = 2L, seed = 42L){
  
  hz         <- horizon.years * 365.25
  d          <- dat[, .(agev = get(age.col), sexv = as.character(get(sex.col)),
                        tf   = get(time.col), ev = get(event.col), x = get(exposure))]
  d[, lab := NA_integer_]
  d[ev == 1 & tf <= hz, lab := 1L]
  d[tf >= hz,           lab := 0L]
  d[, expo := if(inv.transform) INT(as.numeric(x)) else as.numeric(x)]
  d          <- d[!is.na(lab) & !is.na(expo) & !is.na(agev)]
  if(d[lab == 1L, .N] < min.cell) return(NULL)
  d[, strat := paste0(floor(agev / age.bin) * age.bin, "_", sexv)]
  
  ## cross-fit oriented scores: the SIGN is learned out-of-fold, so a PROTECTIVE exposure (higher
  ## value -> lower risk) is not scored as 'negative discrimination'. Scoring the raw exposure would
  ## put AUC below 0.5 for every protective marker (giving spurious negative gains); the fitted score
  ## orients it, mirroring the glm(LP + expo) used by delphi.delta.auc.cv. gain.cv is therefore a
  ## direction-agnostic discrimination gain (>= 0 in expectation); direction lives in the Cox columns.
  set.seed(seed)
  d[, fold := sample(rep_len(1:n.folds, .N)), by = lab]
  sx         <- if(uniqueN(d$sexv) > 1) " + sexv" else ""
  f.cb       <- stats::as.formula("lab ~ expo")                       ## oriented single predictor (stratified path)
  f.d        <- stats::as.formula(paste0("lab ~ agev", sx))           ## demographics            (global fallback baseline)
  f.dB       <- stats::as.formula(paste0("lab ~ agev + expo", sx))    ## demographics + exposure (global fallback combined)
  d[, `:=`(comb.dem = NA_real_, p.d = NA_real_, p.dB = NA_real_)]
  for(k in 1:n.folds){
    tr <- d[fold != k]; te <- d[fold == k]
    d[fold == k, comb.dem := tryCatch(predict(glm(f.cb, tr, family = binomial()), te, type = "link"), error = function(e) NA_real_)]
    d[fold == k, p.d      := tryCatch(predict(glm(f.d,  tr, family = binomial()), te, type = "link"), error = function(e) NA_real_)]
    d[fold == k, p.dB     := tryCatch(predict(glm(f.dB, tr, family = binomial()), te, type = "link"), error = function(e) NA_real_)]
  }
  
  ## stratified: demographics constant within an age x sex cell (AUC 0.5) vs the oriented OOS exposure score
  rs         <- rbindlist(lapply(unique(d$strat), function(s){
    ds <- d[strat == s & !is.na(comb.dem)]
    if(ds[lab == 1L, .N] < min.cell | ds[lab == 0L, .N] < min.cell) return(NULL)
    dd <- dl.cv(ds$lab, rep(0, nrow(ds)), ds$comb.dem)
    data.table(n.case = ds[lab == 1L, .N], est = dd$delta, se = dd$se, a1 = dd$auc1, a2 = dd$auc2)
  }), fill = T)
  
  if(nrow(rs) >= min.strata){
    w        <- rs$n.case
    est      <- sum(w * rs$est) / sum(w)
    se       <- sqrt(sum((w / sum(w))^2 * rs$se^2))
    auc.dem  <- weighted.mean(rs$a1, w); auc.demB <- weighted.mean(rs$a2, w)
    mode     <- "stratified"; nstr <- nrow(rs)
  } else {
    ## global age/sex-adjusted out-of-sample DeLong on ALL cases (already sign-oriented by the glm)
    dg       <- d[!is.na(p.d) & !is.na(p.dB)]
    dd       <- dl.cv(dg$lab, dg$p.d, dg$p.dB)
    est      <- dd$delta; se <- dd$se; auc.dem <- dd$auc1; auc.demB <- dd$auc2
    mode     <- "global_adjusted"; nstr <- 1L
  }
  
  if(is.na(est) | is.na(se)) return(NULL)
  z          <- est / se
  data.table(exposure         = exposure,
             gain.cv          = est, gain.se.cv = se,
             gain.lci.cv      = est - 1.96 * se, gain.uci.cv = est + 1.96 * se,
             gain.pval.cv     = 2 * pnorm(-abs(z)),
             auc.dem.cv       = auc.dem, auc.demB.cv = auc.demB,
             gain.n.strata.cv = nstr, gain.mode = mode)
}