# =============================================================================
# plot_regression_scatter.R
#
# Scatter plot of an exposure (feature) against a continuous outcome, with the
# fitted regression line (linear or spline) overlaid. Companion to
# regression_analysis() — re-fits the model using the SAME pipeline so the
# plot is consistent with the reported estimates.
#
# Key behaviours:
#   - Applies the same transformations used during modelling:
#       * feature RINT + z-score normalisation (if applied in the analysis)
#       * outcome RINT (if outcome_rint = TRUE for that outcome)
#   - Two display scales available for each axis:
#       * "model"    : the transformed scale the model was actually fitted on
#                      (RINT/normalised x, RINT y). The fitted line is straight
#                      for a linear model on this scale.
#       * "original" : the raw data scale. Transformations are inverted for
#                      display, so the fitted curve appears warped accordingly.
#   - Accounts for covariates: by default the scatter shows PARTIAL
#     relationships — both x and y are adjusted for covariates via partial
#     residuals (added back to the fitted component of interest), so the cloud
#     reflects the covariate-adjusted association the model estimated.
#       * covariate_handling = "partial"   partial residuals (default)
#       * covariate_handling = "marginal"  raw scatter, fit at covariate means
#       * covariate_handling = "none"      no covariates in the refit
#   - Only meaningful for continuous outcomes (linear regression). For binary,
#     count, or categorical outcomes the function stops with an informative
#     message (use plot_spline_curve() for probability curves instead).
#
# Dependencies: data.table, rms, splines  (and the helpers from
#   regression_analysis.R, which must be sourced first:
#   parse_covariates, prune_covariate_formula, rint, z_normalise, needs_rint,
#   compute_skewness, x_form_str, is_categorical)
#
# Usage:
#   source("regression_adaptive.R")
#   source("plot_regression_scatter.R")
#
#   res <- regression_analysis(dt = d, outcome = "hba1c",
#                              features = c("age","bmi"),
#                              covariates = ~ sex, outcome_rint = TRUE)
#
#   plot_regression_scatter(
#     feature      = "bmi",
#     outcome      = "hba1c",
#     results      = res,
#     dt           = d,
#     covariates   = ~ sex,
#     outcome_rint = TRUE
#   )
# =============================================================================


# -----------------------------------------------------------------------------
# Helper: resolve an axis transform spec to fn / inv / label (mirrors the
# x_transform helper in plot_spline_associations.R but applies to either axis)
# -----------------------------------------------------------------------------
..resolve_axis_transform <- function(spec) {
  if (is.null(spec) || identical(spec, "identity"))
    return(list(fn = identity, inv = identity, label = NULL, is_identity = TRUE))
  if (is.list(spec)) {
    if (!all(c("fn", "inv") %in% names(spec)))
      stop("transform list must have 'fn' and 'inv'.")
    lbl <- if (!is.null(spec$label)) spec$label else "transformed"
    return(list(fn = spec$fn, inv = spec$inv, label = lbl, is_identity = FALSE))
  }
  switch(as.character(spec),
         "log"   = list(fn = log,   inv = exp,            label = "ln",    is_identity = FALSE),
         "log2"  = list(fn = log2,  inv = function(x) 2^x, label = "log2",  is_identity = FALSE),
         "log10" = list(fn = log10, inv = function(x) 10^x,label = "log10", is_identity = FALSE),
         "sqrt"  = list(fn = sqrt,  inv = function(x) x^2, label = "sqrt",  is_identity = FALSE),
         stop(sprintf("Unknown transform '%s'.", spec))
  )
}


# -----------------------------------------------------------------------------
# Internal: re-fit the linear model for a (feature, outcome) pair and return
# everything needed to draw the scatter + fit.
# -----------------------------------------------------------------------------
.refit_for_scatter <- function(feature, outcome, results, dt,
                               covariates         = NULL,
                               spline_type        = c("rcs","ns"),
                               knots              = 4L,
                               apply_rint_feature = c("auto","always","never"),
                               skew_threshold     = 2,
                               normalise          = TRUE,
                               outcome_rint       = FALSE,
                               model_choice       = c("preferred","spline","linear"),
                               covariate_handling = c("partial","marginal","none"),
                               diag               = FALSE) {
  
  spline_type        <- match.arg(spline_type)
  apply_rint_feature <- match.arg(apply_rint_feature)
  model_choice       <- match.arg(model_choice)
  covariate_handling <- match.arg(covariate_handling)
  k <- as.integer(knots)
  
  # ---- Result rows (long format) for this (outcome, feature) ----
  # Convert to a plain data.frame FIRST so all subsequent [ , $ subsetting is
  # base-R and free of data.table's non-standard evaluation. This is critical:
  # the long-format table has multiple rows per (outcome, feature) — one per
  # term ("overall", "linear", "sp_1", ...). The model-level summary (the
  # preferred model, R2, overall p-value, transform flags) lives ONLY on the
  # term == "overall" row. Reading those fields from any other term row yields
  # wrong values (e.g. a spline-term beta mislabelled as the linear estimate).
  res_df <- as.data.frame(results)
  
  # Coerce key columns to character to avoid factor-level comparison surprises
  oc_col <- as.character(res_df$outcome)
  ft_col <- as.character(res_df$feature)
  tm_col <- as.character(res_df$term)
  
  sel <- which(oc_col == outcome & ft_col == feature)
  if (length(sel) == 0L)
    stop(sprintf("No results row for outcome='%s', feature='%s'.", outcome, feature))
  
  rr_all <- res_df[sel, , drop = FALSE]
  
  # The "overall" row DICTATES model-level behaviour (model choice, R2, p, flags)
  ov_sel <- which(as.character(rr_all$term) == "overall")
  if (length(ov_sel) == 0L)
    stop(sprintf(
      "No 'overall' term row for outcome='%s', feature='%s'. The results table must come from regression_analysis() in long format.",
      outcome, feature))
  if (length(ov_sel) > 1L)
    warning(sprintf(
      "Multiple 'overall' rows for outcome='%s', feature='%s' (%d); using the first. Check for duplicated results.",
      outcome, feature, length(ov_sel)))
  rr <- rr_all[ov_sel[1L], , drop = FALSE]
  
  if (diag) {
    message("=== .refit_for_scatter DIAGNOSTIC ===")
    message(sprintf("  matched rows for outcome='%s' feature='%s': %d",
                    outcome, feature, nrow(rr_all)))
    message(sprintf("  terms present: %s",
                    paste(as.character(rr_all$term), collapse=", ")))
    message(sprintf("  overall-row preferred_model: '%s'", rr$preferred_model[1L]))
    message(sprintf("  overall-row p_overall_association: %.4g",
                    rr$p_overall_association[1L]))
    message(sprintf("  overall-row r2_linear: %.4g  r2_nonlinear: %.4g",
                    rr$r2_linear[1L],
                    if ("r2_nonlinear" %in% names(rr)) rr$r2_nonlinear[1L] else NA_real_))
    message(sprintf("  overall-row spline_knots: %s  rint_feat: %s  outcome_rint: %s  normalised: %s",
                    rr$spline_knots[1L], rr$rint_applied_feature[1L],
                    rr$outcome_rint[1L], rr$normalised[1L]))
    lr <- rr_all[as.character(rr_all$term) == "linear", , drop = FALSE]
    if (nrow(lr) > 0L)
      message(sprintf("  linear-row beta: %.4g  p_value: %.4g",
                      lr$beta[1L], lr$p_value[1L]))
  }
  
  reg_model <- rr$regression_model[1L]
  if (!identical(reg_model, "linear"))
    stop(sprintf(
      "plot_regression_scatter() supports continuous outcomes (linear regression) only.\n  Outcome '%s' used model '%s'. Use plot_spline_curve() for probability curves.",
      outcome, reg_model
    ))
  
  ftype <- rr$feature_type[1L]
  if (identical(ftype, "categorical"))
    stop(sprintf(
      "Feature '%s' is categorical; a scatter vs a continuous exposure is not applicable.\n  Use a boxplot/stripchart of outcome by level instead.",
      feature
    ))
  
  # ---- Sync settings with the analysis (results row is authoritative) ----
  # The analysis stored the knot count and transformation flags it used.
  # If they differ from the function arguments, the refit would NOT match the
  # reported estimates, so we adopt the stored values and warn.
  if ("spline_knots" %in% names(rr) && !is.na(rr$spline_knots[1L])) {
    k_stored <- as.integer(rr$spline_knots[1L])
    if (k_stored != k) {
      warning(sprintf("knots: using %d from results (you passed %d).", k_stored, k))
      k <- k_stored
    }
  }
  if ("rint_applied_feature" %in% names(rr) && !is.na(rr$rint_applied_feature[1L])) {
    # Force feature RINT to match what the analysis actually did
    apply_rint_feature <- if (isTRUE(rr$rint_applied_feature[1L])) "always" else "never"
  }
  if ("normalised" %in% names(rr) && !is.na(rr$normalised[1L])) {
    normalise <- isTRUE(rr$normalised[1L])
  }
  if ("outcome_rint" %in% names(rr) && !is.na(rr$outcome_rint[1L])) {
    outcome_rint <- isTRUE(rr$outcome_rint[1L])
  }
  
  # ---- Covariates ----
  cov_parsed    <- parse_covariates(covariates)
  cov_form      <- if (covariate_handling == "none") "" else cov_parsed$formula_str
  cov_base_vars <- if (covariate_handling == "none") character(0) else cov_parsed$base_vars
  
  # ---- Complete cases ----
  keep_cols <- c(feature, outcome, cov_base_vars)
  keep      <- complete.cases(dt[, keep_cols, with = FALSE])
  
  x_raw <- as.numeric(dt[[feature]][keep])
  y_raw <- as.numeric(dt[[outcome]][keep])
  
  # ---- Apply outcome RINT if used in the analysis ----
  y_proc <- if (outcome_rint) rint(y_raw) else y_raw
  
  # ---- Apply feature transformation (same logic as the analysis) ----
  do_rint <- switch(apply_rint_feature,
                    "always" = TRUE, "never" = FALSE,
                    "auto"   = needs_rint(x_raw, skew_threshold))
  x_after_rint <- if (do_rint) rint(x_raw) else x_raw
  norm_mean <- mean(x_after_rint, na.rm = TRUE)
  norm_sd   <- sd(x_after_rint,   na.rm = TRUE)
  x_proc <- if (normalise && !is.na(norm_sd) && norm_sd > 0)
    (x_after_rint - norm_mean) / norm_sd else x_after_rint
  
  # ---- Covariate data frame ----
  cov_df <- if (length(cov_base_vars) > 0L)
    as.data.frame(dt[keep, cov_base_vars, with = FALSE]) else
      data.frame(row.names = seq_along(x_proc))
  adf <- cbind(data.frame(outcome = y_proc, x = x_proc), cov_df)
  
  # ---- Model form selection ----
  preferred  <- rr$preferred_model[1L]
  use_spline <- switch(model_choice,
                       "preferred" = grepl("non-linear", preferred),
                       "spline"    = TRUE, "linear" = FALSE)
  
  x_str    <- x_form_str(use_spline, spline_type, k)
  f_mod    <- stats::as.formula(paste("outcome ~", x_str, cov_form))
  fit_err  <- NULL
  fit      <- tryCatch(stats::lm(f_mod, data = adf),
                       error = function(e) { fit_err <<- conditionMessage(e); NULL })
  
  spline_fell_back <- FALSE
  if (is.null(fit) && use_spline) {
    # Before giving up on splines, try the OTHER spline basis as a fallback
    # (rcs <-> ns). rms::rcs occasionally fails inside lm() in some sessions;
    # splines::ns is part of base R and very robust.
    alt_type <- if (spline_type == "rcs") "ns" else "rcs"
    x_str_alt <- x_form_str(TRUE, alt_type, k)
    f_alt     <- stats::as.formula(paste("outcome ~", x_str_alt, cov_form))
    fit       <- tryCatch(stats::lm(f_alt, data = adf), error = function(e) NULL)
    if (!is.null(fit)) {
      warning(sprintf(
        "Spline basis '%s' failed for '%s' ~ '%s' (%s); used '%s' instead.",
        spline_type, outcome, feature,
        if (is.null(fit_err)) "unknown error" else fit_err, alt_type))
      x_str       <- x_str_alt
      spline_type <- alt_type
    } else {
      # Both spline bases failed: fall back to linear
      warning(sprintf(
        "Spline refit failed for '%s' ~ '%s' (%s); falling back to linear.",
        outcome, feature, if (is.null(fit_err)) "unknown error" else fit_err))
      x_str <- "x"
      f_mod <- stats::as.formula(paste("outcome ~ x", cov_form))
      fit   <- tryCatch(stats::lm(f_mod, data = adf), error = function(e) NULL)
      use_spline       <- FALSE
      spline_fell_back <- TRUE
    }
  }
  if (is.null(fit)) return(NULL)
  
  # Compare refit R^2 with the value stored by the analysis — a large
  # discrepancy signals a transformation/covariate mismatch.
  refit_r2 <- tryCatch({
    ss_res <- sum(stats::residuals(fit)^2)
    ss_tot <- sum((fit$model[[1L]] - mean(fit$model[[1L]]))^2)
    if (ss_tot > 0) 1 - ss_res / ss_tot else NA_real_
  }, error = function(e) NA_real_)
  stored_r2 <- if (use_spline && "r2_nonlinear" %in% names(rr))
    rr$r2_nonlinear[1L] else if ("r2_linear" %in% names(rr)) rr$r2_linear[1L] else NA_real_
  if (!is.na(refit_r2) && !is.na(stored_r2) && abs(refit_r2 - stored_r2) > 0.02) {
    warning(sprintf(
      "Refit R2 (%.3f) differs from analysis R2 (%.3f) for '%s' ~ '%s'. The covariates/transform arguments may not match the analysis call.",
      refit_r2, stored_r2, outcome, feature))
  }
  
  list(
    fit = fit, adf = adf, x_raw = x_raw, y_raw = y_raw,
    x_proc = x_proc, y_proc = y_proc,
    do_rint = do_rint, norm_mean = norm_mean, norm_sd = norm_sd,
    normalise = normalise, outcome_rint = outcome_rint,
    use_spline = use_spline,
    model_used = if (use_spline) paste0(spline_type, " spline") else "linear",
    cov_form = cov_form, cov_df = cov_df, cov_base_vars = cov_base_vars,
    covariate_handling = covariate_handling,
    k = k, spline_type = spline_type, preferred = preferred,
    rr = rr, rr_all = rr_all,
    spline_fell_back = spline_fell_back, refit_r2 = refit_r2
  )
}


# -----------------------------------------------------------------------------
# Internal: compute partial residuals for x.
#
# Partial residual for the x-component = (component of x in the linear
# predictor) + (overall residual). This isolates the x-relationship after
# adjusting for covariates, and is exactly what is overlaid by the fit line.
#
# Returns a list with:
#   px  : partial-residual y values aligned with x_proc (model scale)
#   fit_grid_x, fit_grid_y : the fitted x-effect curve on the model scale
# -----------------------------------------------------------------------------
.partial_components <- function(mod, n_grid = 200L) {
  fit    <- mod$fit
  adf    <- mod$adf
  x_proc <- mod$x_proc
  
  # Terms prediction isolates each term's contribution to the linear predictor.
  tp <- tryCatch(
    stats::predict(fit, type = "terms"),
    error = function(e) NULL
  )
  
  # Identify the x term column(s) in the terms matrix
  x_term_cols <- if (!is.null(tp)) {
    grep("^x$|rcs\\(x|ns\\(x|^x\\.", colnames(tp), value = TRUE)
  } else character(0)
  
  if (mod$covariate_handling == "partial" && !is.null(tp) &&
      length(x_term_cols) > 0L) {
    # x-component of the linear predictor (sum across spline basis term cols)
    x_component <- rowSums(tp[, x_term_cols, drop = FALSE])
    resid_all   <- stats::residuals(fit)
    px <- x_component + resid_all
  } else {
    # marginal / none: plot raw outcome (model scale)
    px <- mod$y_proc
  }
  
  # ---- Fitted x-effect curve on a grid ----
  x_lo <- quantile(x_proc, 0.01, na.rm = TRUE)
  x_hi <- quantile(x_proc, 0.99, na.rm = TRUE)
  xg   <- seq(x_lo, x_hi, length.out = n_grid)
  
  # Build newdata with covariates at their means / reference levels
  ref_covs <- if (ncol(mod$cov_df) > 0L) {
    ref <- lapply(mod$cov_df, function(col) {
      if (is.numeric(col)) mean(col, na.rm = TRUE)
      else if (is.factor(col)) factor(levels(col)[1L], levels = levels(col))
      else col[1L]
    })
    as.data.frame(ref)[rep(1L, n_grid), , drop = FALSE]
  } else data.frame(row.names = seq_len(n_grid))
  
  newdata <- cbind(data.frame(x = xg), ref_covs)
  
  if (mod$covariate_handling == "partial") {
    # Fitted x-effect only (terms), centred like the partial residuals
    tg <- tryCatch(
      stats::predict(fit, newdata = newdata, type = "terms"),
      error = function(e) NULL
    )
    if (!is.null(tg) && length(x_term_cols) > 0L) {
      yg <- rowSums(tg[, intersect(x_term_cols, colnames(tg)), drop = FALSE])
      # Add back the intercept so the curve sits in the partial-residual cloud
      yg <- yg + stats::coef(fit)[["(Intercept)"]]
    } else {
      yg <- stats::predict(fit, newdata = newdata)
    }
  } else {
    # marginal / none: full prediction at covariate means
    yg <- stats::predict(fit, newdata = newdata)
  }
  
  list(px = px, fit_grid_x = xg, fit_grid_y = yg, x_term_cols = x_term_cols)
}


# =============================================================================
# plot_regression_scatter()
#
# @param feature        feature/exposure column name (numeric)
# @param outcome        continuous outcome column name
# @param results        data.table from regression_analysis()
# @param dt             original data.table
# @param covariates     covariate spec used in the analysis (NULL / vector /
#                       one-sided formula) — must match for a consistent refit
# @param spline_type    "rcs" | "ns" — must match the analysis
# @param knots          integer — must match the analysis
# @param apply_rint_feature "auto"|"always"|"never" — must match the analysis
# @param skew_threshold numeric — must match the analysis
# @param normalise      logical — must match the analysis
# @param outcome_rint   logical — was RINT applied to this outcome?
# @param model_choice   "preferred" (default) | "spline" | "linear"
# @param covariate_handling
#                       "partial"  scatter shows partial residuals adjusted
#                                  for covariates; fit = x-effect curve (default)
#                       "marginal" raw scatter; fit at covariate means
#                       "none"     ignore covariates entirely in the refit
# @param scale          axis display scale:
#                       "model"    transformed scale the model was fitted on
#                                  (RINT/normalised x, RINT y) — straight line
#                                  for linear models (default)
#                       "original" raw data scale; transformations inverted,
#                                  fitted curve warps accordingly
# @param x_transform    optional EXTRA display transform for the x-axis applied
#                       on top of the chosen scale ("log","log10","sqrt", or a
#                       list(fn, inv, label)). NULL = none.
# @param n_grid         fit-line resolution (default 200)
# @param point_col      scatter point colour (default semi-transparent grey)
# @param point_cex      point size (default 0.6)
# @param point_pch      point symbol (default 16)
# @param fit_col        fit line colour (default "steelblue")
# @param fit_lwd        fit line width (default 2.5)
# @param ci             draw 95% CI band around the fit (default TRUE)
# @param ci_col         CI band fill (default fit_col at 15% alpha)
# @param show_linear    overlay the linear fit (dashed) when spline is shown
#                       (default TRUE)
# @param col_linear     colour of the linear overlay (default "firebrick")
# @param xlab,ylab,main labels (auto-generated if NULL)
# @param xlim,ylim      axis limits (auto if NULL)
# @param axis_lwd       axis line width (default 0.5)
# @param add_stats      annotate beta, p, R^2 from results (default TRUE)
# @param ...            passed to the initial plot() call
#
# @return invisibly: list(model, partials, fit_curve)
# =============================================================================
plot_regression_scatter <- function(feature,
                                    outcome,
                                    results,
                                    dt,
                                    covariates         = NULL,
                                    spline_type        = c("rcs","ns"),
                                    knots              = 4L,
                                    apply_rint_feature = c("auto","always","never"),
                                    skew_threshold     = 2,
                                    normalise          = TRUE,
                                    outcome_rint       = FALSE,
                                    model_choice       = c("preferred","spline","linear"),
                                    covariate_handling = c("partial","marginal","none"),
                                    scale              = c("model","original"),
                                    x_transform        = NULL,
                                    n_grid             = 200L,
                                    point_col          = NULL,
                                    point_cex          = 0.6,
                                    point_pch          = 16L,
                                    fit_col            = "steelblue",
                                    fit_lwd            = 2.5,
                                    ci                 = TRUE,
                                    ci_col             = NULL,
                                    show_linear        = TRUE,
                                    col_linear         = "firebrick",
                                    xlab               = NULL,
                                    ylab               = NULL,
                                    main               = NULL,
                                    xlim               = NULL,
                                    ylim               = NULL,
                                    axis_lwd           = 0.5,
                                    add_stats          = TRUE,
                                    diag               = FALSE,
                                    ...) {
  
  spline_type        <- match.arg(spline_type)
  apply_rint_feature <- match.arg(apply_rint_feature)
  model_choice       <- match.arg(model_choice)
  covariate_handling <- match.arg(covariate_handling)
  scale              <- match.arg(scale)
  
  if (is.null(point_col)) point_col <- adjustcolor("grey40", alpha.f = 0.35)
  if (is.null(ci_col))    ci_col    <- adjustcolor(fit_col, alpha.f = 0.15)
  
  xtf <- ..resolve_axis_transform(x_transform)
  
  # ---- Refit ----
  mod <- .refit_for_scatter(
    feature, outcome, results, dt,
    covariates = covariates, spline_type = spline_type, knots = knots,
    apply_rint_feature = apply_rint_feature, skew_threshold = skew_threshold,
    normalise = normalise, outcome_rint = outcome_rint,
    model_choice = model_choice, covariate_handling = covariate_handling,
    diag = diag
  )
  if (is.null(mod)) { warning("Model refit failed."); return(invisible(NULL)) }
  
  comp <- .partial_components(mod, n_grid = n_grid)
  
  # Linear overlay (refit with linear form, same covariate handling)
  lin_comp <- NULL
  if (show_linear && mod$use_spline) {
    mod_lin <- .refit_for_scatter(
      feature, outcome, results, dt,
      covariates = covariates, spline_type = spline_type, knots = knots,
      apply_rint_feature = apply_rint_feature, skew_threshold = skew_threshold,
      normalise = normalise, outcome_rint = outcome_rint,
      model_choice = "linear", covariate_handling = covariate_handling
    )
    if (!is.null(mod_lin)) lin_comp <- .partial_components(mod_lin, n_grid = n_grid)
  }
  
  # ---- Map model-scale x/y to the requested DISPLAY scale ----
  # scale = "model": display = model scale (optionally extra x_transform)
  # scale = "original": invert RINT/normalise for x, invert RINT for y
  to_display_x <- function(x_model) {
    if (scale == "original") {
      # invert normalisation
      xm <- if (mod$normalise && !is.na(mod$norm_sd) && mod$norm_sd > 0)
        x_model * mod$norm_sd + mod$norm_mean else x_model
      # invert RINT via empirical quantile mapping of the training data
      if (mod$do_rint) {
        ord  <- order(mod$x_raw)
        tfun <- stats::approxfun(rint(mod$x_raw)[ord], mod$x_raw[ord], rule = 2L)
        xm   <- tfun(xm)
      }
      x_disp <- xm
    } else {
      x_disp <- x_model
    }
    if (!xtf$is_identity) x_disp <- xtf$fn(x_disp)
    x_disp
  }
  to_display_y <- function(y_model) {
    if (scale == "original" && mod$outcome_rint) {
      ord  <- order(mod$y_raw)
      tfun <- stats::approxfun(rint(mod$y_raw)[ord], mod$y_raw[ord], rule = 2L)
      return(tfun(y_model))
    }
    y_model
  }
  
  # Scatter points (partial residuals or raw, on model scale) -> display
  sx <- to_display_x(mod$x_proc)
  sy <- to_display_y(comp$px)
  
  # Fit curve -> display
  fx <- to_display_x(comp$fit_grid_x)
  fy <- to_display_y(comp$fit_grid_y)
  
  # 95% CI band on the fitted x-effect (model scale), then map y to display
  band <- NULL
  if (ci) {
    ref_covs <- if (ncol(mod$cov_df) > 0L) {
      ref <- lapply(mod$cov_df, function(col) {
        if (is.numeric(col)) mean(col, na.rm = TRUE)
        else if (is.factor(col)) factor(levels(col)[1L], levels = levels(col))
        else col[1L]
      })
      as.data.frame(ref)[rep(1L, n_grid), , drop = FALSE]
    } else data.frame(row.names = seq_len(n_grid))
    nd <- cbind(data.frame(x = comp$fit_grid_x), ref_covs)
    pr <- tryCatch(stats::predict(mod$fit, newdata = nd, se.fit = TRUE),
                   error = function(e) NULL)
    if (!is.null(pr)) {
      lo <- pr$fit - 1.96 * pr$se.fit
      hi <- pr$fit + 1.96 * pr$se.fit
      # For partial scale the curve was centred via intercept; predict() gives
      # the full prediction, so shift band to match the drawn fy by aligning
      # midpoints.
      shift <- mean(fy - pr$fit, na.rm = TRUE)
      band <- list(x = fx,
                   lo = to_display_y(lo + shift),
                   hi = to_display_y(hi + shift))
    }
  }
  
  # ---- Labels ----
  if (is.null(xlab)) {
    base_x <- feature
    tags <- c()
    if (scale == "model") {
      if (mod$do_rint)   tags <- c(tags, "RINT")
      if (mod$normalise) tags <- c(tags, "z")
    }
    if (!xtf$is_identity) tags <- c(tags, xtf$label)
    xlab <- if (length(tags)) sprintf("%s (%s)", base_x, paste(tags, collapse=", ")) else base_x
  }
  if (is.null(ylab)) {
    ylab <- if (scale == "model" && mod$outcome_rint) sprintf("%s (RINT)", outcome) else outcome
    if (mod$covariate_handling == "partial" && length(mod$cov_base_vars) > 0L)
      ylab <- sprintf("%s | partial", ylab)
  }
  if (is.null(main)) main <- sprintf("%s ~ %s  [%s]", outcome, feature, mod$model_used)
  
  if (is.null(xlim)) xlim <- range(sx, fx, na.rm = TRUE, finite = TRUE)
  if (is.null(ylim)) {
    ylim <- range(sy, fy, if (!is.null(band)) c(band$lo, band$hi),
                  na.rm = TRUE, finite = TRUE)
  }
  
  # ---- Draw ----
  plot(NA, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab,
       main = main, las = 1L, axes = FALSE, ...)
  axis(1L, lwd = axis_lwd, lwd.ticks = axis_lwd)
  axis(2L, lwd = axis_lwd, lwd.ticks = axis_lwd, las = 1L)
  box(lwd = axis_lwd)
  
  # Scatter
  points(sx, sy, pch = point_pch, cex = point_cex, col = point_col)
  
  # CI band
  if (!is.null(band)) {
    ord <- order(band$x)
    polygon(c(band$x[ord], rev(band$x[ord])),
            c(band$lo[ord], rev(band$hi[ord])),
            col = ci_col, border = NA)
  }
  
  # Linear overlay (dashed)
  if (!is.null(lin_comp)) {
    lfx <- to_display_x(lin_comp$fit_grid_x)
    lfy <- to_display_y(lin_comp$fit_grid_y)
    ord <- order(lfx)
    lines(lfx[ord], lfy[ord], col = col_linear, lty = 2L, lwd = max(1, fit_lwd - 1))
  }
  
  # Main fit
  ord <- order(fx)
  lines(fx[ord], fy[ord], col = fit_col, lwd = fit_lwd)
  
  # ---- Stats annotation ----
  # beta / p_value live on TERM rows (linear, sp_1, ...), NOT the overall row.
  # Model-level p-values (p_overall_association, p_nonlinear, R2) live on every
  # row but are most reliably read from the overall row.
  if (add_stats) {
    rr     <- mod$rr        # overall row: model-level summaries (data.frame)
    rr_all <- mod$rr_all    # all term rows (data.frame)
    fmt_p  <- function(p) if (is.na(p)) "NA" else if (p < 1e-3) sprintf("%.1e", p) else sprintf("%.3f", p)
    
    # beta of the linear term (always reported, even when spline is preferred)
    lin_row <- rr_all[as.character(rr_all$term) == "linear", , drop = FALSE]
    beta    <- if (nrow(lin_row) > 0L) lin_row$beta[1L] else NA_real_
    p_lin   <- if (nrow(lin_row) > 0L) lin_row$p_value[1L] else NA_real_
    
    # Model-level stats from the overall row
    p_oa  <- rr$p_overall_association[1L]
    p_adj <- if ("p_overall_adj" %in% names(rr)) rr$p_overall_adj[1L] else NA_real_
    p_nl  <- if ("p_nonlinear" %in% names(rr)) rr$p_nonlinear[1L] else NA_real_
    # Use the R2 matching the model being drawn
    r2 <- if (mod$use_spline && "r2_nonlinear" %in% names(rr))
      rr$r2_nonlinear[1L] else rr$r2_linear[1L]
    
    txt <- c(
      sprintf("model: %s", mod$preferred),
      if (!is.na(beta))  sprintf("beta(lin) = %.3g", beta),
      if (!is.na(p_oa))  sprintf("p(overall) = %s", fmt_p(p_oa)),
      if (!is.na(p_adj)) sprintf("p.adj = %s", fmt_p(p_adj)),
      if (!is.na(r2))    sprintf("R2 = %.3f", r2),
      if (!is.na(p_nl))  sprintf("p(non-lin) = %s", fmt_p(p_nl))
    )
    if (isTRUE(mod$spline_fell_back))
      txt <- c(txt, "(spline refit FAILED)")
    legend("topright", legend = txt, bty = "n", cex = 0.74, text.col = "grey20")
  }
  
  # ---- Legend for fit lines ----
  leg_lab <- mod$model_used; leg_col <- fit_col; leg_lty <- 1L; leg_lwd <- fit_lwd
  if (!is.null(lin_comp)) {
    leg_lab <- c(leg_lab, "linear")
    leg_col <- c(leg_col, col_linear)
    leg_lty <- c(leg_lty, 2L)
    leg_lwd <- c(leg_lwd, max(1, fit_lwd - 1))
  }
  legend("topleft", legend = leg_lab, col = leg_col, lty = leg_lty,
         lwd = leg_lwd, bty = "n", cex = 0.8)
  
  invisible(list(model = mod, partials = comp,
                 fit_curve = data.frame(x = fx, y = fy)))
}