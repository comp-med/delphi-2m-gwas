library(data.table)

#' Multi-method expression enrichment testing for a gene list across
#' tissues or cell types.
#'
#' Works with any long-format expression table — bulk tissue (e.g. Human
#' Protein Atlas consensus RNA) or single-cell / single-nucleus pseudobulk
#' data — as long as it has three columns: one for gene identifiers, one for
#' the grouping label (tissue, cell type, cluster …), and one for the
#' expression value (nTPM, CPM, mean counts, …).
#'
#' Four complementary tests are run per group:
#'
#'   1. Fisher's exact test (one-sided, "greater") on specifically-expressed
#'      genes.  A gene is called "specific" to a group when its expression in
#'      that group exceeds `specificity_expr_min` AND its expression-ratio
#'      relative to the mean across all other groups exceeds
#'      `specificity_fold`.  The 2x2 contingency table is:
#'
#'                          specific  |  not specific
#'        query genes    :    a       |      b
#'        background     :    c       |      d
#'
#'      Tests whether query genes are over-represented among specifically
#'      expressed genes.
#'
#'   2. Wilcoxon rank-sum test (one-sided, "greater") comparing raw expression
#'      of query genes vs. background genes within each group.  Captures
#'      graded, non-specific upregulation across all genes without requiring
#'      a specificity threshold.
#'
#'   3. Permutation mean-expression test.  Observed mean expression of the
#'      query set is compared against a null distribution built by resampling
#'      equally-sized gene sets from the full gene universe.  Yields an
#'      empirical p-value and a z-score for effect size.
#'
#'   4. KS enrichment score (GSEA-style).  Genes are ranked by decreasing
#'      expression within each group; a running enrichment score tracks
#'      whether query genes cluster near the top of this rank list.
#'
#' @param expr_dt            data.table in long format.  Must contain the
#'                           columns named by `gene_col`, `group_col`, and
#'                           `expr_col`.
#' @param query_genes        Character vector of gene identifiers to test.
#' @param gene_col           Name of the gene-identifier column. Default "gene".
#' @param group_col          Name of the grouping column (tissue, cell type …).
#'                           Default "tissue".
#' @param expr_col           Name of the expression-value column.
#'                           Default "expression".
#'
#' @param specificity_fold   Numeric >= 1. Minimum fold-change of a gene's
#'                           expression in the focal group relative to the mean
#'                           expression across all *other* groups for the gene
#'                           to be called specifically expressed.  Default 2.
#'                           Increase (e.g. 5) for stricter tissue-specificity
#'                           comparable to HPA "tissue enhanced" calls.
#'                           Set to 1 to require only `specificity_expr_min`.
#'                           When `log_fold = TRUE` this threshold is compared
#'                           against log1p(focal) - log1p(mean_other), so it is
#'                           interpreted as a log1p fold-change (e.g. 1 ≈ a
#'                           ~2.7-fold difference on the raw scale).
#' @param specificity_expr_min  Numeric >= 0. Minimum absolute expression
#'                           value in the focal group (same units as `expr_col`)
#'                           for a gene to be called specifically expressed.
#'                           Prevents near-zero/noisy ratios from inflating
#'                           fold-change.  Default 1.
#' @param log_fold          Logical. If TRUE, fold-change for specificity is
#'                           computed in log1p space:
#'                             log1p(focal_expr) - log1p(mean_other_expr)
#'                           and compared against `specificity_fold`.
#'                           Recommended for nCPM, CPM, or any right-skewed
#'                           count-based unit where raw ratios are unstable.
#'                           If FALSE (default), raw ratio focal / mean_other
#'                           is used, which is appropriate for nTPM bulk data.
#'
#' @param n_perm             Number of permutations for the permutation mean
#'                           test and the KS permutation null.  Default 1 000.
#' @param ks_abs             Logical. If TRUE use the maximum absolute running
#'                           sum as the enrichment score (mirrors classic GSEA).
#'                           Default TRUE.
#' @param adjust_method      p-value adjustment method passed to p.adjust(),
#'                           applied across groups within each method.
#'                           Default "BH".
#' @param min_query_n        Minimum number of query genes present in a group
#'                           for that group to be tested.  Groups with fewer
#'                           query genes receive NA.  Default 3.
#' @param seed               Random seed for reproducibility.  Default 42.
#'
#' @return A named list with five elements:
#'   $fisher      — data.table: per-group Fisher test + specificity counts
#'   $wilcoxon    — data.table: per-group Wilcoxon results
#'   $permutation — data.table: per-group permutation results
#'   $ks          — data.table: per-group KS enrichment results
#'   $summary     — data.table: all adjusted p-values merged and sorted by
#'                              median adjusted p-value across all four methods
#'
test_expression_enrichment <- function(
    expr_dt,
    query_genes,
    gene_col             = "gene",
    group_col            = "tissue",
    expr_col             = "expression",
    specificity_fold     = 2,
    specificity_expr_min = 1,
    log_fold             = FALSE,
    n_perm               = 1000L,
    ks_abs               = TRUE,
    adjust_method        = "BH",
    min_query_n          = 3L,
    seed                 = 42L
) {
  
  ## ── input validation ───────────────────────────────────────────────────────
  stopifnot(is.data.table(expr_dt))
  for (col in c(gene_col, group_col, expr_col))
    if (!col %in% names(expr_dt))
      stop(sprintf("Column '%s' not found in expr_dt.", col))
  stopifnot(is.character(query_genes), length(query_genes) > 0)
  stopifnot(specificity_fold >= 1, specificity_expr_min >= 0)
  stopifnot(is.logical(log_fold), length(log_fold) == 1L)
  stopifnot(n_perm >= 1L, min_query_n >= 1L)
  
  set.seed(seed)
  
  ## work on a copy with canonical internal names to keep downstream code clean
  dt <- copy(expr_dt[, .SD, .SDcols = c(gene_col, group_col, expr_col)])
  setnames(dt, c(gene_col, group_col, expr_col), c("gene", "group", "expr"))
  
  all_genes   <- unique(dt$gene)
  query_genes <- intersect(query_genes, all_genes)
  n_query     <- length(query_genes)
  n_total     <- length(all_genes)
  
  if (n_query < min_query_n)
    stop(sprintf(
      "Only %d query gene(s) matched in the expression table (min_query_n = %d).",
      n_query, min_query_n))
  
  message(sprintf("Query genes matched : %d", n_query))
  message(sprintf("Total genes in table: %d", n_total))
  message(sprintf("Groups to test      : %d", uniqueN(dt$group)))
  
  dt[, is_query := gene %in% query_genes]
  
  ## ── pre-compute per-gene × group specificity flag ─────────────────────────
  ## Specifically expressed = expr >= specificity_expr_min  AND
  ##   expr / mean(expr across all OTHER groups) >= specificity_fold
  ##
  ## mean_other is derived from the gene's grand sum so no nested by() needed:
  ##   mean_other = (grand_sum - focal_expr) / (n_groups_observed - 1)
  
  message(sprintf(
    "Computing specificity flags (%sfold >= %.2f, min expr >= %.2f) ...",
    if (log_fold) "log1p " else "",
    specificity_fold, specificity_expr_min))
  
  gene_agg <- dt[, .(grand_sum = sum(expr, na.rm = TRUE),
                     n_obs     = .N),
                 by = gene]
  
  dt <- merge(dt, gene_agg, by = "gene", all.x = TRUE)
  
  dt[, mean_other := fcase(
    n_obs > 1L, (grand_sum - expr) / (n_obs - 1L),
    n_obs == 1L, 0                     # single group: no "other" exists
  )]
  
  ## fold-change metric: raw ratio OR log1p difference depending on log_fold.
  ## log1p space: log1p(focal) - log1p(mean_other)
  ##   - compresses the right tail of CPM/nCPM distributions
  ##   - zeros become 0 after log1p, so a gene expressed only in the focal
  ##     group gets log1p(expr) - 0 = log1p(expr), a finite and interpretable
  ##     value rather than Inf
  ##   - specificity_fold is then a log1p-scale threshold (e.g. 1 ≈ 2.7x raw)
  if (log_fold) {
    dt[, fold_over_other := log1p(expr) - log1p(mean_other)]
  } else {
    dt[, fold_over_other := fcase(
      mean_other >  0, expr / mean_other,
      mean_other == 0 & expr > 0, Inf,   # only expressed here → infinite fold
      default = 0
    )]
  }
  
  dt[, is_specific := expr >= specificity_expr_min &
       fold_over_other >= specificity_fold]
  
  ## ── helper: KS running enrichment score ───────────────────────────────────
  ## Expects is_query_sorted to be a logical vector ordered by decreasing expr.
  .ks_score <- function(is_query_sorted, nq) {
    n         <- length(is_query_sorted)
    hit_step  <-  1 / nq
    miss_step <- -1 / (n - nq)
    running   <- cumsum(ifelse(is_query_sorted, hit_step, miss_step))
    if (ks_abs) running[which.max(abs(running))] else max(running)
  }
  
  ## ── 1. Fisher's exact test ────────────────────────────────────────────────
  message("Running Fisher's exact tests ...")
  
  fisher_res <- dt[, {
    nq <- sum(is_query)
    if (nq < min_query_n) {
      .(n_query             = nq,
        n_specific_query    = NA_integer_,
        n_specific_bg       = NA_integer_,
        specific_query_genes = NA_character_,
        odds_ratio          = NA_real_,
        p_value             = NA_real_)
    } else {
      idx_a <- is_query &  is_specific    # query   & specifically expressed
      idx_b <- is_query & !is_specific    # query   & not specifically expressed
      idx_c <- !is_query &  is_specific   # background & specifically expressed
      idx_d <- !is_query & !is_specific   # background & not specifically expressed
      
      a <- sum(idx_a); b <- sum(idx_b)
      c <- sum(idx_c); d <- sum(idx_d)
      
      gene_list <- paste(sort(gene[idx_a]), collapse = "|")
      
      if ((a + c) == 0L) {
        # no specifically expressed genes in this group → p = 1, OR = NA
        .(n_query             = nq,
          n_specific_query    = 0L,
          n_specific_bg       = 0L,
          specific_query_genes = NA_character_,
          odds_ratio          = NA_real_,
          p_value             = 1)
      } else {
        ft <- fisher.test(matrix(c(a, c, b, d), nrow = 2L),
                          alternative = "greater")
        .(n_query             = nq,
          n_specific_query    = a,
          n_specific_bg       = c,
          specific_query_genes = if (a > 0L) gene_list else NA_character_,
          odds_ratio          = ft$estimate,
          p_value             = ft$p.value)
      }
    }
  }, by = group]
  
  fisher_res[!is.na(p_value),
             p_adj := p.adjust(p_value, method = adjust_method)]
  setorder(fisher_res, p_value)
  
  ## ── 2. Wilcoxon rank-sum test ─────────────────────────────────────────────
  message("Running Wilcoxon rank-sum tests ...")
  
  wilcox_res <- dt[, {
    q_expr  <- expr[ is_query]
    bg_expr <- expr[!is_query]
    nq      <- length(q_expr)
    if (nq < min_query_n) {
      .(n_query     = nq,
        mean_query  = NA_real_,
        mean_bg     = NA_real_,
        log2fc_mean = NA_real_,
        statistic   = NA_real_,
        p_value     = NA_real_)
    } else {
      wt <- wilcox.test(q_expr, bg_expr, alternative = "greater", exact = FALSE)
      mq <- mean(q_expr,  na.rm = TRUE)
      mb <- mean(bg_expr, na.rm = TRUE)
      .(n_query     = nq,
        mean_query  = mq,
        mean_bg     = mb,
        log2fc_mean = log2((mq + 1) / (mb + 1)),
        statistic   = wt$statistic,
        p_value     = wt$p.value)
    }
  }, by = group]
  
  wilcox_res[!is.na(p_value),
             p_adj := p.adjust(p_value, method = adjust_method)]
  setorder(wilcox_res, p_value)
  
  ## ── 3. Permutation mean-expression test ───────────────────────────────────
  message(sprintf("Running permutation mean tests (%d permutations) ...", n_perm))
  
  perm_res <- dt[, {
    q_expr <- expr[is_query]
    nq     <- length(q_expr)
    if (nq < min_query_n) {
      .(n_query   = nq,
        obs_mean  = NA_real_,
        perm_mean = NA_real_,
        perm_sd   = NA_real_,
        z_score   = NA_real_,
        p_value   = NA_real_)
    } else {
      obs_mean  <- mean(q_expr, na.rm = TRUE)
      all_expr  <- expr
      
      null_means <- vapply(seq_len(n_perm), function(i)
        mean(all_expr[sample.int(length(all_expr), nq)], na.rm = TRUE),
        numeric(1L))
      
      pm  <- mean(null_means)
      psd <- sd(null_means)
      z   <- if (psd > 0) (obs_mean - pm) / psd else NA_real_
      ep  <- (sum(null_means >= obs_mean) + 1L) / (n_perm + 1L)
      
      .(n_query   = nq,
        obs_mean  = obs_mean,
        perm_mean = pm,
        perm_sd   = psd,
        z_score   = z,
        p_value   = ep)
    }
  }, by = group]
  
  perm_res[!is.na(p_value),
           p_adj := p.adjust(p_value, method = adjust_method)]
  setorder(perm_res, p_value)
  
  ## ── 4. KS enrichment score ────────────────────────────────────────────────
  message("Running KS enrichment scoring ...")
  
  ks_res <- dt[, {
    ord          <- order(expr, decreasing = TRUE)
    query_sorted <- is_query[ord]
    nq           <- sum(query_sorted)
    if (nq < min_query_n) {
      .(n_query = nq, es = NA_real_, p_value = NA_real_)
    } else {
      obs_es <- .ks_score(query_sorted, nq)
      
      null_es <- vapply(seq_len(n_perm), function(i) {
        pq        <- logical(length(query_sorted))
        pq[sample.int(length(query_sorted), nq)] <- TRUE
        .ks_score(pq, nq)
      }, numeric(1L))
      
      emp_p <- (sum(null_es >= obs_es) + 1L) / (n_perm + 1L)
      .(n_query = nq, es = obs_es, p_value = emp_p)
    }
  }, by = group]
  
  ks_res[!is.na(p_value),
         p_adj := p.adjust(p_value, method = adjust_method)]
  setorder(ks_res, p_value)
  
  ## ── 5. Summary table ──────────────────────────────────────────────────────
  summary_dt <- Reduce(
    function(a, b) merge(a, b, by = c("group", "n_query"), all = TRUE),
    list(
      fisher_res[, .(group, n_query,
                     fisher_padj          = p_adj,
                     fisher_OR            = odds_ratio,
                     n_specific_query     = n_specific_query,
                     specific_query_genes = specific_query_genes)],
      wilcox_res[, .(group, n_query,
                     wilcox_padj      = p_adj,
                     wilcox_log2fc    = log2fc_mean)],
      perm_res[,   .(group, n_query,
                     perm_padj        = p_adj,
                     perm_zscore      = z_score)],
      ks_res[,     .(group, n_query,
                     ks_padj          = p_adj,
                     ks_es            = es)]
    )
  )
  
  ## restore user-supplied column name for the group variable
  setnames(summary_dt, "group", group_col)
  for (tbl in list(fisher_res, wilcox_res, perm_res, ks_res))
    setnames(tbl, "group", group_col)
  
  p_cols <- c("fisher_padj", "wilcox_padj", "perm_padj", "ks_padj")
  summary_dt[, median_padj := apply(.SD, 1L, function(x) {
    v <- as.numeric(x)
    median(v[!is.na(v)])
  }), .SDcols = p_cols]
  
  setorder(summary_dt, median_padj)
  
  list(
    fisher      = fisher_res,
    wilcoxon    = wilcox_res,
    permutation = perm_res,
    ks          = ks_res,
    summary     = summary_dt
  )
}


## ══════════════════════════════════════════════════════════════════════════════
##  Example usage — bulk tissue (Human Protein Atlas, nTPM scale)
## ══════════════════════════════════════════════════════════════════════════════
##
## Download: https://www.proteinatlas.org/about/download → "RNA consensus tissue"
##
## hpa_raw   <- fread("rna_consensus_tissue_gene_data.tsv")
## expr_bulk <- hpa_raw[, .(gene       = `Gene name`,
##                           tissue     = Tissue,
##                           expression = nTPM)]
##
## query_genes <- c("APOE", "LDLR", "PCSK9", "LPA", "ABCA1", "ABCG8",
##                  "SCARB1", "LPL", "LIPC", "APOB")
##
## res_bulk <- test_expression_enrichment(
##   expr_dt              = expr_bulk,
##   query_genes          = query_genes,
##   gene_col             = "gene",
##   group_col            = "tissue",
##   expr_col             = "expression",
##   specificity_fold     = 2,     # raw ratio: >= 2x over mean of other tissues
##   specificity_expr_min = 1,     # >= 1 nTPM in the focal tissue
##   log_fold             = FALSE, # raw ratio is stable for nTPM bulk data
##   n_perm               = 1000L,
##   adjust_method        = "BH"
## )
##
## print(res_bulk$summary[median_padj < 0.05])
##
## for (nm in names(res_bulk))
##   fwrite(res_bulk[[nm]], sprintf("enrichment_bulk_%s.tsv", nm), sep = "\t")


## ══════════════════════════════════════════════════════════════════════════════
##  Example usage — single-cell pseudobulk (nCPM scale)
## ══════════════════════════════════════════════════════════════════════════════
##
## nCPM is right-skewed with frequent zeros across cell types, so raw ratios
## are noisy.  Setting log_fold = TRUE computes fold-change in log1p space:
##   log1p(focal_nCPM) - log1p(mean_other_nCPM)
## This compresses the tail and handles zeros gracefully (log1p(0) = 0).
## The specificity_fold threshold is then on the log1p scale:
##   ~1.0  ≈  e^1 - 1 ≈ 1.7x on raw scale  (permissive)
##   ~1.5  ≈  e^1.5-1 ≈ 3.5x on raw scale  (moderate, recommended default)
##   ~2.0  ≈  e^2 - 1 ≈ 6.4x on raw scale  (strict)
##
## sc_raw  <- fread("pseudobulk_mean_expression_per_celltype.tsv")
## expr_sc <- sc_raw[, .(gene      = gene_symbol,
##                        cell_type = cell_type,
##                        nCPM      = nCPM)]
##
## res_sc <- test_expression_enrichment(
##   expr_dt              = expr_sc,
##   query_genes          = query_genes,
##   gene_col             = "gene",
##   group_col            = "cell_type",
##   expr_col             = "nCPM",
##   specificity_fold     = 1.5,   # log1p-scale threshold (~3.5x on raw nCPM)
##   specificity_expr_min = 5,     # >= 5 nCPM in focal cell type (noise floor)
##   log_fold             = TRUE,  # essential for right-skewed nCPM data
##   n_perm               = 1000L,
##   adjust_method        = "BH"
## )
##
## print(res_sc$summary[median_padj < 0.05])
##
## for (nm in names(res_sc))
##   fwrite(res_sc[[nm]], sprintf("enrichment_sc_%s.tsv", nm), sep = "\t")