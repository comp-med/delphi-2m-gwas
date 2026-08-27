####################################################################
#### "what is worth adding to Delphi-2M" -- points 1-3           ##
####  (1) the embeddings reconstruct SOME, not all, continuous    ##
####      risk factors;                                           ##
####  (2) part of Delphi-2M's accuracy reflects that             ##
####      reconstruction, so adding those markers adds little;    ##
####  (3) the residual value concentrates at the disease each     ##
####      marker DIAGNOSTICALLY DEFINES, not at low-pleiotropy or  ##
####      poorly-reconstructed markers (GDF15 = absorbed foil).   ##
####  Runs jointly over biomarkers AND proteins.                  ##
#### Maik Pietzner                                   22/06/2026   ##
####################################################################

## --> packages needed <-- ##
require(data.table)

#-------------------------------------------------#
##--   per-feature builder (one feature set)     --##
#-------------------------------------------------#
## dt.auc/dt.cox/dt.recon as in the head script; *.id give the feature key column in each.
feature.specificity <- function(dt.auc, dt.cox, dt.recon, auc.id, cox.id, recon.id,
                                feature.type, p.thresh, case.floor){
  
  A          <- copy(as.data.table(dt.auc));   setnames(A, auc.id,   "feature")
  C          <- copy(as.data.table(dt.cox));   setnames(C, cox.id,   "feature")
  R          <- copy(as.data.table(dt.recon)); setnames(R, recon.id, "feature")
  if(!"mode"    %in% names(A)) A[, mode := "stratified"]
  if(!"gain.cv" %in% names(A)) A[, gain.cv := NA_real_]
  
  A          <- A[sex %in% c("men", "women")]
  C          <- C[sex %in% c("men", "women")]
  if(!"category" %in% names(R)) R[, category := NA_character_]
  R          <- R[sex %in% c("men", "women"), .(feature, sex, recon.r2, category)]
  
  ## pleiotropy: distinct diseases the feature is associated with (demographics-adjusted Cox),
  ## restricted to the features actually analysed (present in the delta-AUC table)
  feats.keep <- unique(A$feature)
  keepcols   <- c("feature", "sex", "event", "c.delta")
  if("description" %in% names(C)) keepcols <- c(keepcols, "description")
  sig        <- C[feature %in% feats.keep & !is.na(`Pr(>|z|)`) & `Pr(>|z|)` < p.thresh, ..keepcols]
  if(!"description" %in% names(sig)) sig[, description := NA_character_]
  plei       <- sig[, .(pleiotropy = uniqueN(event)), by = .(feature, sex)]
  ## primary disease, chosen by Cox strength (NOT by ΔAUC -> no selection on the outcome)
  primary    <- sig[order(-c.delta), .(event.primary = event[1], event.primary.desc = description[1]),
                    by = .(feature, sex)]
  
  ## residual over Delphi on robust pairs (stratified, case-floored)
  As         <- A[mode == "stratified" & n.case >= case.floor,
                  .(feature, sex, event, delta.auc.cv, gain.cv, auc.base.cv, auc.comb.cv, n.case)]
  res.prim   <- merge(primary, As, by.x = c("feature", "sex", "event.primary"),
                      by.y = c("feature", "sex", "event"), all.x = T)
  res.agg    <- As[, .(residual.max = max(delta.auc.cv), residual.med = median(delta.auc.cv),
                       gain.med = median(gain.cv), n.case.max = max(n.case)), by = .(feature, sex)]
  
  dt         <- Reduce(function(a, b) merge(a, b, by = c("feature", "sex"), all = T),
                       list(plei,
                            res.prim[, .(feature, sex, event.primary, event.primary.desc,
                                         residual.primary = delta.auc.cv,
                                         auc.base.primary = auc.base.cv,
                                         auc.comb.primary = auc.comb.cv,
                                         n.case.primary   = n.case)],
                            res.agg, R))
  dt[, feature.type := feature.type]
  dt[!is.na(pleiotropy)]
}

## --> compact disease label: "Date E80 first reported (..)" -> ICD code "E80" <-- ##
icd.code   <- function(x){
  x <- gsub("^Date\\s+", "", gsub("\\s+first reported.*$", "", x))
  sub("\\s+.*$", "", x)
}

## --> 3-class grouping of the best marker per disease: <-- ##
##  1 diagnostic analyte | 2 organ-function / physiological measure | 3 secreted / tissue-leakage protein
marker.class <- function(feature, type){
  f   <- tolower(feature)
  fun <- "fev1|forced_expiratory|blood_pressure|bone_mineral|bmd|haemoglobin|hemoglobin|haematocrit|platelet|erythrocyte|distribution_width|eosinophil|lymphocyte|monocyte|neutrophil|basophil|reticulocyte|leukocyte|white_blood|red_blood|fat_mass|lean_mass|^bmi$|waist|hip_circ|whr|grip|spiro|fvc|pulse_rate|heart_rate"
  ifelse(type == "protein", 3L, ifelse(grepl(fun, f), 2L, 1L))
}


#-------------------------------------------------#
##--   joint report over all feature sets        --##
#-------------------------------------------------#
## sets: named list; each = list(auc, cox, recon, auc.id, cox.id, recon.id, highlight)
improvement.report <- function(sets, p.thresh = 3.287959e-06, case.floor = 100L,
                               adds.threshold = 0.02, graphics.dir = "../graphics", out.dir = ".."){
  
  ds         <- format(Sys.Date(), "%Y%m%d")
  col.type   <- c(biomarker = "#4C72B0", protein = "#C44E52")
  ramp       <- colorRampPalette(c("#DD8452", "#E8E8E8", "#4C72B0"))(101)
  col.recon  <- function(r) ramp[pmin(101L, pmax(1L, round(r / 0.5 * 100) + 1L))]
  hl         <- unlist(lapply(sets, function(s) s$highlight))
  
  ## assemble the three views across all sets
  feat       <- rbindlist(lapply(names(sets), function(nm){
    s <- sets[[nm]]
    feature.specificity(s$auc, s$cox, s$recon, s$auc.id, s$cox.id, s$recon.id, nm, p.thresh, case.floor)
  }), fill = T)
  feat[, highlighted := toupper(feature) %in% toupper(hl)]
  feat[, icd := icd.code(event.primary.desc)]
  feat[is.na(icd) | icd == "", icd := event.primary]
  
  recon.all  <- rbindlist(lapply(names(sets), function(nm){
    R <- copy(as.data.table(sets[[nm]]$recon)); setnames(R, sets[[nm]]$recon.id, "feature")
    R <- R[sex %in% c("men", "women"), .(recon.r2 = mean(recon.r2, na.rm = T)), by = feature]
    R[, feature.type := nm]; R
  }), fill = T)
  
  ## within-sex reconstruction for every feature (individual-level supplement)
  recon.bysex <- rbindlist(lapply(names(sets), function(nm){
    R <- copy(as.data.table(sets[[nm]]$recon)); setnames(R, sets[[nm]]$recon.id, "feature")
    if(!"category" %in% names(R)) R[, category := NA_character_]
    R <- R[sex %in% c("men", "women"), .(feature, sex, recon.r2, category)]
    R[, feature.type := nm]; R
  }), fill = T)
  
  ## per-pair table: every analysed marker-disease pair (individual-level supplement)
  pairs.all  <- rbindlist(lapply(names(sets), function(nm){
    s <- sets[[nm]]
    A <- copy(as.data.table(s$auc)); setnames(A, s$auc.id, "feature")
    if(!"mode"        %in% names(A)) A[, mode := "stratified"]
    if(!"gain.cv"     %in% names(A)) A[, gain.cv := NA_real_]
    if(!"description" %in% names(A)) A[, description := NA_character_]
    for(cc in c("lci.cv", "uci.cv", "pval.cv")) if(!cc %in% names(A)) A[, (cc) := NA_real_]
    R <- copy(as.data.table(s$recon)); setnames(R, s$recon.id, "feature")
    R <- R[sex %in% c("men", "women"), .(feature, sex, recon.r2)]
    A <- A[sex %in% c("men", "women") & mode == "stratified" & n.case >= case.floor,
           .(feature, sex, event, description, delta.auc.cv, lci.cv, uci.cv, pval.cv, gain.cv, n.case)]
    A <- merge(A, R, by = c("feature", "sex"), all.x = T)
    A[, feature.type := nm]; A
  }), fill = T)
  
  sp         <- function(d, xc, yc) d[is.finite(get(xc)) & is.finite(get(yc)),
                                      .(r = round(cor(get(xc), get(yc), method = "spearman"), 3), n = .N), by = feature.type]
  
  w          <- function(x, nm) fwrite(x, file.path(out.dir, paste0(nm, ".", ds, ".txt")),
                                       sep = "\t", row.names = F, quote = F, na = NA)
  
  #-- POINT 1: reconstruction -- some, not all --------------------------------#
  cat("\n================ POINT 1 | reconstruction of risk factors ================\n")
  s1.recon   <- recon.all[, .(n = .N, median = round(median(recon.r2), 3),
                              iqr.lo = round(quantile(recon.r2, .25), 3), iqr.hi = round(quantile(recon.r2, .75), 3),
                              max = round(max(recon.r2), 3), pct.ge.0.2 = round(mean(recon.r2 >= 0.2), 2),
                              pct.ge.0.3 = round(mean(recon.r2 >= 0.3), 2)), by = feature.type]
  print(s1.recon)
  cat("-- best reconstructed (per type) --\n")
  print(recon.all[order(-recon.r2), head(.SD, 6), by = feature.type, .SDcols = c("feature", "recon.r2")])
  
  #-- POINT 2: absorption -- on AVERAGE a marker adds little (marker-centric) -#
  cat("\n================ POINT 2 | added value over Delphi-2M (marker-centric) ================\n")
  s2.absorb  <- pairs.all[, .(pairs = .N,
                              med.gain.over.dem = round(median(gain.cv, na.rm = T), 3),
                              med.delta.over.delphi = round(median(delta.auc.cv, na.rm = T), 3),
                              pct.absorbed = round(1 - median(delta.auc.cv, na.rm = T) / median(gain.cv, na.rm = T), 2),
                              pct.adds.gt.005 = round(mean(delta.auc.cv > 0.05, na.rm = T), 3),
                              pct.adds.gt.thr = round(mean(delta.auc.cv > adds.threshold, na.rm = T), 3)), by = feature.type]
  rd         <- sp(pairs.all, "recon.r2", "delta.auc.cv"); setnames(rd, "r", "sp.recon.delta")
  rg         <- sp(pairs.all, "recon.r2", "gain.cv");      setnames(rg, "r", "sp.recon.gain")
  s2.absorb  <- merge(merge(s2.absorb, rd[, .(feature.type, sp.recon.delta)], by = "feature.type", all.x = T),
                      rg[, .(feature.type, sp.recon.gain)], by = "feature.type", all.x = T)
  print(s2.absorb)
  
  #-- DISEASE FLOOR: best single marker per disease (disease-centric) ---------#
  ## the fair counterpart to absorption -- for each disease, how much does its BEST
  ## available marker add over Delphi-2M? (max over markers -> optimistic, winner's-cursed,
  ## but OOS per pair; this is the single-marker ceiling, NOT the average)
  cat("\n================ DISEASE FLOOR | best single marker per disease ================\n")
  db         <- pairs.all[is.finite(delta.auc.cv)]
  ev.desc    <- unique(db[!is.na(description), .(event, description)], by = "event")
  best.all   <- db[order(-delta.auc.cv), .(best.delta = delta.auc.cv[1], best.lci = lci.cv[1], best.uci = uci.cv[1],
                                           best.pval = pval.cv[1], best.marker = feature[1],
                                           best.type = feature.type[1], n.markers = .N), by = .(event, sex)]
  best.bio   <- db[feature.type == "biomarker"][order(-delta.auc.cv),
                                                .(best.delta.bio = delta.auc.cv[1], best.marker.bio = feature[1]), by = .(event, sex)]
  best.pro   <- db[feature.type == "protein"][order(-delta.auc.cv),
                                              .(best.delta.pro = delta.auc.cv[1], best.marker.pro = feature[1]), by = .(event, sex)]
  disease.best <- Reduce(function(a, b) merge(a, b, by = c("event", "sex"), all = T),
                         list(best.all, best.bio, best.pro))
  disease.best <- merge(disease.best, ev.desc, by = "event", all.x = T)
  disease.best[, icd := icd.code(description)]
  disease.best[is.na(icd) | icd == "", icd := event]
  disease.best[, marker.class := marker.class(best.marker, best.type)]
  disease.best[, class.label   := c("diagnostic analyte", "organ function", "secreted protein")[marker.class]]
  disease.best <- disease.best[order(-best.delta)]
  s4.disease <- data.table(
    n.diseases                     = nrow(disease.best),
    med.best.delta                 = round(median(disease.best$best.delta, na.rm = T), 3),
    n.gt.002                       = sum(disease.best$best.delta > 0.02, na.rm = T),
    n.gt.005                       = sum(disease.best$best.delta > 0.05, na.rm = T),
    n.gt.010                       = sum(disease.best$best.delta > 0.10, na.rm = T),
    pct.gt.005                     = round(mean(disease.best$best.delta > 0.05, na.rm = T), 3),
    best.is.protein.among.gt.005   = sum(disease.best[best.delta > 0.05]$best.type == "protein", na.rm = T),
    best.is.biomarker.among.gt.005 = sum(disease.best[best.delta > 0.05]$best.type == "biomarker", na.rm = T))
  print(s4.disease)
  cat("-- most improvable diseases (best marker -> ICD); cf. a single pleiotropic protein is rarely best --\n")
  print(head(disease.best[, .(icd, sex, best.delta, best.marker, best.type,
                              best.delta.bio, best.marker.bio, best.delta.pro, best.marker.pro, n.markers)], 25))
  
  #-- POINT 3: single-marker residual at the disease it DEFINES ---------------#
  cat("\n================ POINT 3 | single-marker residual at the defining disease ================\n")
  s3.resid   <- feat[is.finite(residual.primary),
                     .(features = .N, median.residual.primary = round(median(residual.primary), 3),
                       pct.adds.gt.005 = round(mean(residual.primary > 0.05), 3),
                       pct.adds.gt.thr = round(mean(residual.primary > adds.threshold), 3)), by = feature.type]
  print(s3.resid)
  cat("-- the exceptions: top residual at the defining disease (marker -> ICD) --\n")
  print(head(feat[is.finite(residual.primary)][order(-residual.primary),
                                               .(feature, feature.type, sex, recon.r2, pleiotropy, residual.primary, icd)], 20))
  cat("-- highlighted features --\n")
  print(feat[highlighted == T, .(feature, feature.type, sex, recon.r2, pleiotropy,
                                 residual.primary, icd)][order(-residual.primary)])
  cat("-- NB pleiotropy/reconstruction ~ residual are POSITIVE and confounded; not an explanatory axis --\n")
  cat("   pleiotropy ~ residual.primary:\n"); print(sp(feat, "pleiotropy", "residual.primary"))
  cat("   recon.r2  ~ residual.primary:\n"); print(sp(feat, "recon.r2", "residual.primary"))
  
  #-- write supplementary tables (summary + individual level) -----------------#
  w(s1.recon,                                       "Improvement.summary.reconstruction")
  w(s2.absorb,                                      "Improvement.summary.absorption")
  w(s4.disease,                                     "Improvement.summary.disease_floor")
  w(s3.resid,                                       "Improvement.summary.residual")
  w(recon.bysex[order(feature.type, -recon.r2)],    "Improvement.table.reconstruction_byfeature")
  w(pairs.all[order(feature.type, -delta.auc.cv)],  "Improvement.table.pairs")
  w(disease.best,                                   "Improvement.table.disease_best_marker")
  w(feat[order(-residual.primary)],                 "Improvement.table.feature_specificity")
  
  #-- FIGURE 1: reconstruction distribution (point 1) ------------------------#
  pdf(file.path(graphics.dir, paste0("Reconstruction.distribution.", ds, ".pdf")), width = 10, height = 5)
  par(mfrow = c(1, length(sets)), mar = c(4.6, 4.8, 2.6, 1.2), mgp = c(2.7, 0.6, 0), tcl = -0.3, las = 1)
  for(nm in names(sets)){
    d <- recon.all[feature.type == nm][order(recon.r2)]
    d[, q := (1:.N) / .N]
    plot(d$q, d$recon.r2, type = "l", lwd = 2, col = col.type[nm], ylim = c(0, max(0.5, max(d$recon.r2))),
         xlab = paste0("rank fraction of ", nm, "s"), ylab = "reconstruction R\u00B2 (within sex)", main = "")
    abline(h = median(d$recon.r2), col = "grey70", lty = 2)
    lab <- rbind(tail(d, 5), d[toupper(feature) %in% toupper(hl)])
    text(lab$q, lab$recon.r2, labels = toupper(lab$feature), pos = 2, cex = 0.6, col = "grey25")
    mtext(sprintf("%ss: median R\u00B2 = %.2f", nm, median(d$recon.r2)), side = 3, line = 0.4, adj = 0, cex = 0.9)
  }
  dev.off()
  
  #-- FIGURE 2: residual distribution, reconstruction as colour (point 3) ----#
  ## near-spike at 0 (absorbed) with a thin tail of diagnostic readouts; tail labelled marker -> ICD
  pdf(file.path(graphics.dir, paste0("Residual.distribution.", ds, ".pdf")), width = 10, height = 5)
  par(mfrow = c(1, length(sets)), mar = c(4.6, 1.2, 2.6, 1.2), mgp = c(2.7, 0.6, 0), tcl = -0.3, las = 1)
  for(nm in names(sets)){
    d <- feat[feature.type == nm & is.finite(residual.primary)]
    set.seed(1L); d[, yj := runif(.N, -1, 1)]
    plot(d$residual.primary, d$yj, pch = 21, cex = 0.9, bg = col.recon(d$recon.r2), col = "grey50",
         xlab = "residual \u0394AUC over Delphi-2M (at the defining disease)", ylab = "", yaxt = "n",
         xlim = c(min(0, min(d$residual.primary)), max(d$residual.primary) * 1.18), main = "")
    abline(v = 0, col = "grey80"); abline(v = adds.threshold, col = "grey60", lty = 3)
    tail.d <- d[residual.primary > quantile(residual.primary, 0.92) | highlighted == T]
    text(tail.d$residual.primary, tail.d$yj, labels = paste0(toupper(tail.d$feature), " \u2192 ", tail.d$icd),
         pos = 4, cex = 0.58, offset = 0.3, col = "grey20")
    mtext(sprintf("%ss: median %.3f, %d%% add > %.2f", nm, median(d$residual.primary),
                  round(100 * mean(d$residual.primary > adds.threshold)), adds.threshold),
          side = 3, line = 0.4, adj = 0, cex = 0.85)
    if(nm == names(sets)[length(sets)])
      legend("bottomright", legend = c("recon R\u00B2 ~0", "~0.25", "~0.5"), pt.bg = ramp[c(1, 51, 101)],
             pch = 21, col = "grey50", bty = "n", cex = 0.7, title = "reconstruction")
  }
  dev.off()
  
  invisible(list(features = feat[], reconstruction = recon.all[], reconstruction.bysex = recon.bysex[],
                 pairs = pairs.all[], disease.best = disease.best[],
                 summary = list(reconstruction = s1.recon[], absorption = s2.absorb[],
                                disease.floor = s4.disease[], residual = s3.resid[])))
}

#-------------------------------------------------#
##--   MAIN figure: 5a reconstruction | 5b      --##
##--   absorption | 5c pleiotropy/absorption |   --##
##--   5d forest of substantively-improved dx    --##
#-------------------------------------------------#
## res = the (invisible) list returned by improvement.report()
main.figure  <- function(res, adds.threshold = 0.02, delta.min = 0.05, p.thresh = NULL,
                         label.markers = c("gdf15", "eda2r", "tnfrsf10b", "plaur", "igfbp4", "bmi"),
                         graphics.dir = "../graphics", file = NULL){
  
  require(data.table)
  ds         <- format(Sys.Date(), "%Y%m%d")
  if(is.null(file)) file <- file.path(graphics.dir, paste0("Main.figure.Delphi.improvement.", ds, ".pdf"))
  recon      <- as.data.table(res$reconstruction); pairs <- as.data.table(res$pairs)
  feat       <- as.data.table(res$features);       db    <- as.data.table(res$disease.best)
  col.type   <- c(biomarker = "#4C72B0", protein = "#C44E52")
  col.sex    <- c(men = "#4C72B0", women = "#C44E52")
  col.cls    <- c("#5B8C5A", "#C8A24B", "#9B5DA0")            # class 1 / 2 / 3
  cls.lab    <- c("1 diagnostic\nanalyte", "2 organ\nfunction", "3 secreted\nprotein")
  
  ## per-marker pleiotropy / absorption summary (panel c)
  M          <- feat[, .(pleiotropy = mean(pleiotropy, na.rm = T), gain = mean(gain.med, na.rm = T),
                         residual = mean(residual.med, na.rm = T)), by = .(feature, feature.type)]
  
  pdf(file, width = 11, height = 11)
  layout(matrix(c(1, 2, 3, 4, 4, 4), nrow = 2, byrow = T), heights = c(1, 2.7))
  par(mar = c(4.4, 4.6, 3.0, 1.0), mgp = c(2.6, 0.6, 0), tcl = -0.3, las = 1)
  
  ## (a) reconstruction is partial --------------------------------------------#
  plot(NA, xlim = c(0, 1), ylim = c(0, max(0.5, max(recon$recon.r2, na.rm = T))),
       xlab = "rank fraction of features", ylab = "reconstruction R\u00B2 (within sex)", main = "")
  for(nm in names(col.type)){
    d <- recon[feature.type == nm][order(recon.r2)]; d[, q := (1:.N) / .N]
    lines(d$q, d$recon.r2, lwd = 2.4, col = col.type[nm])
    lab <- unique(rbind(tail(d, 3), d[toupper(feature) %in% toupper(label.markers)]), by = "feature")
    if(nrow(lab)) text(lab$q, lab$recon.r2, toupper(lab$feature), pos = 2, cex = 0.5, col = "grey30")
  }
  legend("topleft", legend = paste0(names(col.type), "s"), col = col.type, lwd = 2.4, bty = "n", cex = 0.8)
  mtext("a", 3, line = 1.4, adj = 0, font = 2, cex = 1.2); mtext("reconstruction is partial", 3, line = 0.3, adj = 0, cex = 0.78)
  
  ## (b) near-total absorption ------------------------------------------------#
  d          <- pairs[is.finite(delta.auc.cv)]
  br         <- seq(floor(min(d$delta.auc.cv) * 100) / 100, ceiling(max(d$delta.auc.cv) * 100) / 100, by = 0.0025)
  h          <- hist(d$delta.auc.cv, breaks = br, plot = FALSE)
  plot(h, col = "grey82", border = NA, xlim = c(-0.02, 0.10),
       xlab = "\u0394AUC over Delphi-2M (per pair)", ylab = "marker\u2013disease pairs", main = "")
  abline(v = median(d$delta.auc.cv), col = "#C44E52", lwd = 2); abline(v = adds.threshold, lty = 3, col = "grey50")
  mtext("b", 3, line = 1.4, adj = 0, font = 2, cex = 1.2)
  mtext(sprintf("median %.3f; %.0f%% add \u2265 %.2f", median(d$delta.auc.cv),
                100 * mean(d$delta.auc.cv >= adds.threshold), adds.threshold), 3, line = 0.3, adj = 0, cex = 0.78)
  
  ## (c) breadth of risk is absorbed -- two smooth fits, the shaded gap IS the absorbed value -#
  d          <- M[is.finite(pleiotropy) & pleiotropy > 0 & is.finite(gain) & is.finite(residual)][order(pleiotropy)]
  lx         <- log10(d$pleiotropy)
  fg         <- lowess(lx, d$gain, f = 0.6); fr <- lowess(lx, d$residual, f = 0.6)
  fr.y       <- approx(fr$x, fr$y, xout = fg$x, rule = 2)$y
  xg         <- 10^fg$x
  plot(NA, log = "x", xlim = range(d$pleiotropy), ylim = c(0, max(fg$y) * 1.15),
       xlab = "number of diseases predicted (pleiotropy)", ylab = "AUC gain", main = "")
  polygon(c(xg, rev(xg)), c(fg$y, rev(fr.y)), col = adjustcolor("#C44E52", 0.16), border = NA)
  lines(xg, fg$y, col = "#4C72B0", lwd = 2.8)
  lines(xg, fr.y, col = "grey45", lwd = 2.2, lty = 2)
  text(10^quantile(lx, 0.5), max(fg$y) * 0.42, "absorbed by\nDelphi-2M", col = "#9A3B40", cex = 0.72, font = 3)
  ann        <- M[toupper(feature) %in% c("GDF15", "LOG_CRP_CLEANED")]
  if(nrow(ann) > 0){
    lab.map <- c(GDF15 = "GDF-15", LOG_CRP_CLEANED = "CRP")
    points(ann$pleiotropy, ann$gain, pch = 21, bg = "#DD8452", col = "black", cex = 1.2)
    text(ann$pleiotropy, ann$gain, lab.map[toupper(ann$feature)], pos = 2, cex = 0.62, col = "grey20")
  }
  legend("topleft", legend = c("gain over age + sex", "residual over Delphi-2M"), lwd = c(2.8, 2.2), lty = c(1, 2),
         col = c("#4C72B0", "grey45"), bty = "n", cex = 0.7)
  mtext("c", 3, line = 1.4, adj = 0, font = 2, cex = 1.2); mtext("breadth of risk is absorbed", 3, line = 0.3, adj = 0, cex = 0.78)
  
  ## (d) forest of substantively-improved diseases (>= delta.min) -------------#
  if(!"marker.class" %in% names(db)) db[, marker.class := marker.class(best.marker, best.type)]
  for(cc in c("best.lci", "best.uci")) if(!cc %in% names(db)) db[, (cc) := NA_real_]
  fo         <- db[is.finite(best.delta) & best.delta >= delta.min]
  if(!is.null(p.thresh) && "best.pval" %in% names(fo)) fo <- fo[is.na(best.pval) | best.pval < p.thresh]
  par(mar = c(4.4, 13.5, 1.6, 1.0))
  if(nrow(fo) == 0){ plot.new(); text(0.5, 0.5, "no disease with \u0394AUC \u2265 delta.min"); dev.off(); return(invisible(file)) }
  setorder(fo, marker.class, best.delta)
  yy <- numeric(nrow(fo)); pos <- 0; prev <- fo$marker.class[1]
  for(i in seq_len(nrow(fo))){ if(fo$marker.class[i] != prev){ pos <- pos + 2; prev <- fo$marker.class[i] }; pos <- pos + 1; yy[i] <- pos }
  fo[, y := yy]
  xlim       <- c(0, max(c(fo$best.uci, fo$best.delta), na.rm = T) * 1.05)
  plot(NA, xlim = xlim, ylim = c(0, max(fo$y) + 1), yaxt = "n", main = "",
       xlab = "\u0394AUC over Delphi-2M (best marker per disease)", ylab = "")
  abline(v = delta.min, lty = 3, col = "grey70")
  segments(fo$best.lci, fo$y, fo$best.uci, fo$y, col = col.sex[fo$sex], lwd = 1.5)
  points(fo$best.delta, fo$y, pch = 21, bg = col.sex[fo$sex], col = "grey25", cex = 1.0)
  mtext(paste0(fo$icd, ": ", toupper(fo$best.marker)), side = 2, at = fo$y, las = 1, line = 0.3, cex = 0.5)
  for(cl in sort(unique(fo$marker.class))){
    yr <- fo[marker.class == cl, range(y)]
    mtext(cls.lab[cl], side = 2, at = mean(yr), las = 1, line = 10.2, cex = 0.6, font = 2, col = col.cls[cl])
    if(cl != max(fo$marker.class)) abline(h = max(yr) + 1, col = "grey90")
  }
  legend("bottomright", legend = names(col.sex), pt.bg = col.sex, pch = 21, col = "grey25", bty = "n", cex = 0.75)
  mtext("d", 3, line = 0.3, adj = 0, font = 2, cex = 1.2)
  dev.off()
  invisible(file)
}

#-------------------------------------------------#
##--  FINAL analysis: does Delphi absorb MORE   --##
##--  of a marker's value the more diseases it   --##
##--  predicts? loss = gain(dem) - residual.     --##
##--  GDF-15 = archetype of the pleiotropic,     --##
##--  highly-absorbed morbidity protein.         --##
#-------------------------------------------------#
## res = the (invisible) list returned by improvement.report()
absorption.vs.pleiotropy <- function(res, gain.floor = 0.01,
                                     label.markers = c("gdf15", "tnfrsf10b", "plaur", "adm", "tnfrsf1a",
                                                       "igfbp4", "log_crp_cleaned", "log_cysc_cleaned"),
                                     graphics.dir = "../graphics", out.dir = ".."){
  
  require(data.table)
  ds         <- format(Sys.Date(), "%Y%m%d")
  col.type   <- c(biomarker = "#4C72B0", protein = "#C44E52")
  D          <- as.data.table(res$features)
  
  ## one row per marker (mean over sexes); loss = gain over demographics - residual over Delphi
  M          <- D[, .(pleiotropy = mean(pleiotropy, na.rm = T), gain = mean(gain.med, na.rm = T),
                      residual = mean(residual.med, na.rm = T), recon = mean(recon.r2, na.rm = T)),
                  by = .(feature, feature.type)]
  M[, loss := gain - residual]
  M[, absorption.frac := ifelse(gain > 0, 1 - residual / gain, NA_real_)]
  
  scor       <- function(d, x, y){
    d <- d[is.finite(get(x)) & is.finite(get(y))]
    if(nrow(d) > 5) round(cor(d[[x]], d[[y]], method = "spearman"), 3) else NA_real_
  }
  pcor       <- function(d, x, y, z){           # partial Spearman of (x,y) given z, on ranks
    d <- d[is.finite(get(x)) & is.finite(get(y)) & is.finite(get(z))]
    if(nrow(d) < 6) return(NA_real_)
    rx <- rank(d[[x]]); ry <- rank(d[[y]]); rz <- rank(d[[z]])
    round(cor(residuals(lm(rx ~ rz)), residuals(lm(ry ~ rz))), 3)
  }
  
  cat("\n================ FINAL | does Delphi absorb more of pleiotropic markers? ================\n")
  tab        <- rbindlist(lapply(unique(M$feature.type), function(t){
    d  <- M[feature.type == t]; df <- d[gain > gain.floor]
    data.table(feature.type = t, n = nrow(d),
               rho.pleio.gain          = scor(d,  "pleiotropy", "gain"),      # pleiotropic = stronger
               rho.pleio.loss          = scor(d,  "pleiotropy", "loss"),      # absolute AUC absorbed
               rho.pleio.absorbed.frac = scor(df, "pleiotropy", "absorption.frac"),
               rho.pleio.frac.adj.gain = pcor(df, "pleiotropy", "absorption.frac", "gain"))
  }))
  print(tab)
  cat("-- GDF-15 and the 10 most pleiotropic proteins (loss = gain absorbed by Delphi) --\n")
  print(M[feature.type == "protein"][order(-pleiotropy)][1:10,
                                                         .(feature, pleiotropy = round(pleiotropy), gain = round(gain, 3), residual = round(residual, 3),
                                                           loss = round(loss, 3), absorption.frac = round(absorption.frac, 2), recon = round(recon, 2))])
  
  fwrite(M[order(-pleiotropy)], file.path(out.dir, paste0("Improvement.table.absorption_vs_pleiotropy.", ds, ".txt")),
         sep = "\t", row.names = F, quote = F, na = NA)
  
  #-- FIGURE: gain (filled) and surviving residual (grey) vs pleiotropy --------#
  ## the gap between the two curves IS the loss; it widens with pleiotropy while the residual stays ~0
  pdf(file.path(graphics.dir, paste0("Absorption.vs.pleiotropy.", ds, ".pdf")), width = 10, height = 5)
  par(mfrow = c(1, length(unique(M$feature.type))), mar = c(4.6, 4.8, 2.6, 1.2), mgp = c(2.7, 0.6, 0), tcl = -0.3, las = 1)
  for(t in unique(M$feature.type)){
    d <- M[feature.type == t & is.finite(pleiotropy) & pleiotropy > 0 & is.finite(gain)]
    plot(d$pleiotropy, d$gain, log = "x", pch = 21, bg = adjustcolor(col.type[t], 0.45), col = "grey55", cex = 0.8,
         xlab = "number of diseases predicted (pleiotropy)", ylab = "AUC gain over age + sex",
         ylim = c(0, max(d$gain, na.rm = T) * 1.05), main = "")
    points(d$pleiotropy, d$residual, pch = 20, col = "grey70", cex = 0.5)
    if(nrow(d) > 15){
      lw <- lowess(log10(d$pleiotropy), d$gain);     lines(10^lw$x, lw$y, col = col.type[t], lwd = 2.4)
      l2 <- lowess(log10(d$pleiotropy), d$residual); lines(10^l2$x, l2$y, col = "grey45", lwd = 2.2, lty = 2)
    }
    lab <- d[toupper(feature) %in% toupper(label.markers)]
    if(nrow(lab) > 0){
      points(lab$pleiotropy, lab$gain, pch = 21, bg = "#DD8452", col = "black", cex = 1.2, lwd = 1.1)
      text(lab$pleiotropy, lab$gain, toupper(lab$feature), pos = 2, cex = 0.6, col = "grey20")
    }
    mtext(sprintf("%ss: gain \u2197 with pleiotropy, residual over Delphi stays \u22480", t), side = 3, line = 0.4, adj = 0, cex = 0.82)
    if(t == unique(M$feature.type)[1])
      legend("topleft", legend = c("gain over age+sex", "residual over Delphi-2M"), lwd = c(2.4, 2.2), lty = c(1, 2),
             col = c(col.type[t], "grey45"), bty = "n", cex = 0.72)
  }
  dev.off()
  invisible(M[order(-pleiotropy)])
}