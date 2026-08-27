####################################################################
#### reconstruction vs headroom -- the per-pair comparison       ##
#### sourced by the head script; call once for exposures, once    ##
#### for proteins. Prints the correlation three ways (max-over-    ##
#### disease [old], per-pair, per-pair + case floor) and, when     ##
#### gain.cv is present, the salience (gain~recon) vs residual      ##
#### (dAUC~recon) contrast; draws the corrected scatter.          ##
#### Maik Pietzner                                   19/06/2026   ##
####################################################################

require(data.table)

## dt.pairs : per (exposure, disease, sex) delta-AUC table (res.auc.delta / res.auc.proteins / cox.auc.reporting)
## dt.recon : per (exposure, sex) reconstruction R2 (recon.exposure / recon.proteins)
## id.col   : the exposure key in dt.pairs ("short_name" for phenotypes, "exposure" for proteins)
## y.col    : headroom over Delphi (default cross-fitted delta.auc.cv)
## case.col : incident-case count used for the floor
## restrict.continuous : if a 'type' column is present in dt.recon, keep numeric/integer only
recon.vs.headroom <- function(dt.pairs, dt.recon,
                              id.col = "short_name", y.col = "delta.auc.cv",
                              case.col = "n.case", case.floor = 200L,
                              gain.col = "gain.cv", group = "exposures",
                              graphics.dir = "../graphics"){
  
  P          <- copy(as.data.table(dt.pairs))
  R          <- copy(as.data.table(dt.recon))
  
  ## harmonise the exposure key to 'short_name'
  if(id.col != "short_name") setnames(P, id.col, "short_name")
  
  ## reconstruction side: within sex, optionally continuous-only, carry category
  if(!"category" %in% names(R)) R[, category := NA_character_]
  if("type" %in% names(R))      R <- R[type %in% c("continuous")]
  R          <- R[sex %in% c("men", "women"), .(short_name, sex, recon.r2, category)]
  print(head(R))
  
  ## pairs side: within sex, keep the columns we need
  P          <- P[sex %in% c("men", "women")]
  if(!case.col %in% names(P)) P[, (case.col) := NA_integer_]
  keep       <- c("short_name", "sex", y.col, case.col, intersect(gain.col, names(P)))
  P          <- P[, ..keep]
  
  mp         <- merge(R, P, by = c("short_name", "sex"))
  if(nrow(mp) == 0){ cat("[", group, "] no overlap between pairs and reconstruction\n"); return(invisible(NULL)) }
  
  #-- correlation, three ways ------------------------------------------------#
  cor.by.sex <- function(d, xc, yc, lbl) d[is.finite(get(xc)) & is.finite(get(yc)),
                                           .(group = group, set = lbl,
                                             r = round(cor(get(xc), get(yc)), 3), n = .N), by = sex]
  
  ## (1) max over disease per exposure -- reproduces the inflated, winner's-cursed value
  mx         <- mp[is.finite(get(y.col)), .(max.dauc = max(get(y.col)), recon.r2 = recon.r2[1]), by = .(short_name, sex)]
  t.max      <- cor.by.sex(mx, "recon.r2", "max.dauc", "1 max over disease (old)")
  ## (2) per pair, all predictive
  t.pair     <- cor.by.sex(mp, "recon.r2", y.col, "2 per pair (all predictive)")
  ## (3) per pair, with the case floor
  t.flr      <- cor.by.sex(mp[get(case.col) >= case.floor], "recon.r2", y.col, paste0("3 per pair, n.case>=", case.floor))
  cor.tab    <- rbindlist(list(t.max, t.pair, t.flr), use.names = T)
  
  ## (4) salience vs residual, if the biomarker's own gain is present
  has.gain   <- gain.col %in% names(mp)
  if(has.gain){
    g.sal <- cor.by.sex(mp[get(case.col) >= case.floor], "recon.r2", gain.col, paste0("4 gain~recon, n>=", case.floor))
    cor.tab <- rbindlist(list(cor.tab, g.sal), use.names = T)
  }
  cat("\n--- reconstruction vs headroom [", group, "] ---\n", sep = "")
  print(cor.tab[order(set, sex)])
  
  #-- figure -----------------------------------------------------------------#
  pal        <- c("Biomarker"="#C44E52","Body composition"="#4C72B0","Blood cell counts"="#8172B3",
                  "Pulmonary"="#55A868","Bone"="#CCB974","Cardiovascular"="#DD8452","Proteomics"="#4C72B0")
  col.for    <- function(cat) ifelse(is.na(cat) | !(cat %in% names(pal)), "#9C9C9C", pal[cat])
  
  df         <- mp[get(case.col) >= case.floor & is.finite(recon.r2) & is.finite(get(y.col))]
  if(nrow(df) > 0){
    df[, cx := 0.6 + 1.4 * (log10(get(case.col)) - min(log10(get(case.col)))) / max(diff(range(log10(get(case.col)))), 1e-6)]
    nrow.fig <- if(has.gain) 2 else 1
    pdf(file.path(graphics.dir, paste0("Recon_vs_headroom_perpair.", group, ".", format(Sys.Date(), "%Y%m%d"), ".pdf")),
        width = 8, height = 4 * nrow.fig)
    par(mfrow = c(nrow.fig, 2), mar = c(4.6, 4.8, 2.6, 1.2), mgp = c(2.7, 0.6, 0), tcl = -0.3, las = 1)
    draw <- function(d, yy, ylab){
      plot(d$recon.r2, d[[yy]], type = "n", xlab = "within-sex reconstruction R\u00B2", ylab = ylab, main = "")
      abline(h = 0, col = "grey80")
      points(d$recon.r2, d[[yy]], pch = 21, cex = d$cx, bg = adjustcolor(col.for(d$category), 0.7), col = "grey35")
      if(nrow(d) > 5) lines(lowess(d$recon.r2, d[[yy]]), col = "#C44E52", lwd = 2)
      mtext(sprintf("%s   r = %.2f   (n = %d)", d$sex[1], cor(d$recon.r2, d[[yy]]), nrow(d)), side = 3, line = 0.5, adj = 0, cex = 0.9)
    }
    for(s in c("men", "women")) draw(df[sex == s], y.col, "\u0394AUC over Delphi (cross-fitted)")
    if(has.gain) for(s in c("men", "women")) draw(df[sex == s & is.finite(get(gain.col))], gain.col, "gain over demographics")
    dev.off()
  }
  invisible(mp[])
}