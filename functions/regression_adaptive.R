# =============================================================================
# regression_analysis.R
#
# Unified regression analysis across multiple features and multiple outcomes,
# with automatic outcome type detection and model selection.
#
# OUTCOME AUTODETECTION
#   Continuous (numeric, >2 unique values, not integer-count-like)
#     -> Linear regression  (stats::lm)
#        Optional RINT of outcome, F-tests for association/non-linearity
#   Binary (2 unique non-NA values, integer/logical/factor)
#     -> Logistic regression  (stats::glm, family = binomial)
#        LRT for association/non-linearity; OR and 95% CI reported
#   Count (non-negative integer, no upper bound, overdispersion tested)
#     -> Negative binomial regression  (MASS::glm.nb)
#        LRT for association/non-linearity; IRR and 95% CI reported
#   Categorical (factor/character with 3+ levels)
#     -> Multinomial logistic regression  (nnet::multinom)
#        One set of log-OR / OR per outcome level vs. reference
#        LRT for overall association via model comparison
#   Override via outcome_model = list(my_outcome = "linear") etc.
#
# EXPOSURE (feature) TYPES
#   Numeric  : RINT + z-score normalisation (optional), linear + spline models
#   Categorical : factor with user-defined reference level, omnibus test
#
# OUTPUT: LONG FORMAT — one row per (outcome, feature, term [, outcome_level])
#
#   outcome              outcome variable name
#   outcome_level        for multinomial: which level vs. reference (NA otherwise)
#   regression_model     "linear" | "logistic" | "negbinom" | "multinomial"
#   feature              feature/exposure variable name
#   feature_type         "numeric" | "categorical"
#   term                 "overall" | "linear" | "sp_1".."sp_N" | "<cat_level>"
#   term_type            "overall" | "linear" | "spline" | "categorical_level"
#   beta                 log-scale coefficient (log-OR, log-IRR, or raw beta)
#   se                   standard error of beta
#   lower95, upper95     95% CI on beta scale
#   effect               exponentiated effect (OR/IRR) for logistic/negbinom/
#                        multinomial; same as beta for linear
#   effect_lower95,
#   effect_upper95       exponentiated CI
#   p_value              Wald p for individual terms (NA for spline/overall)
#   p_overall_association  LRT/F-test: full model vs. null/covariates
#   p_overall_adj          BH-adjusted p_overall_association
#   p_overall_spline       LRT/F-test: spline vs. null (numeric only)
#   p_nonlinear            LRT/F-test: spline vs. linear (numeric only)
#   p_overall_spline_adj, p_nonlinear_adj
#   n_total, n_missing
#   skewness_raw           (numeric features only)
#   rint_applied_feature, outcome_rint, normalised, spline_knots
#   ref_level              reference level (categorical features or multinomial)
#   aic_linear, aic_nonlinear
#   r2_linear, r2_nonlinear   (linear models only; Nagelkerke pseudo-R2 for others)
#   preferred_model        "linear" | "non-linear (spline)" | etc.
#   error
#
# DEPENDENCIES
#   R packages: data.table, rms, splines, MASS, nnet
#
#   This script is SELF-CONTAINED. All shared helper functions
#   (parse_covariates, prune_covariate_formula, rint, z_normalise,
#    needs_rint, compute_skewness, make_sp_na_dt, x_form_str,
#    extract_x_coefs, is_categorical, prepare_factor) are defined below
#   in SECTION 0. No other scripts need to be sourced first.
# =============================================================================

# install.packages(c("data.table", "rms", "MASS", "nnet"))
library(data.table)
library(rms)
library(splines)
library(MASS)
library(nnet)


# =============================================================================
# SECTION 0: SHARED HELPER FUNCTIONS
#
# These were previously sourced from logistic_feature_analysis.R and
# linear_feature_analysis.R.  They are inlined here so this script runs
# standalone.  Definitions are identical to the canonical versions.
# =============================================================================

# -----------------------------------------------------------------------------
# Helper: parse covariate specification
#
# Accepts either:
#   (a) character vector  c("age", "sex")         -> "+ age + sex"
#   (b) one-sided formula ~ age + sex + age:sex   -> "+ age + sex + age:sex"
#       or                 ~ age * sex             -> "+ age * sex"
#
# Returns a list:
#   $formula_str  character string to paste after the feature term, e.g.
#                 "+ age + sex + age:sex"  (empty string if no covariates)
#   $base_vars    character vector of the *column names* actually needed from
#                 dt for complete-case filtering (just the variable names,
#                 not interaction terms like "age:sex")
# -----------------------------------------------------------------------------
parse_covariates <- function(covariates) {
  if (is.null(covariates)) {
    return(list(formula_str = "", base_vars = character(0)))
  }
  
  if (is.character(covariates)) {
    # Traditional vector: plain additive formula
    return(list(
      formula_str = paste("+", paste(covariates, collapse = " + ")),
      base_vars   = covariates
    ))
  }
  
  if (inherits(covariates, "formula")) {
    # One-sided formula ~ a + b + a:b  or  ~ a * b
    rhs <- covariates[[length(covariates)]]   # RHS of formula
    
    # formula_str: deparse the RHS and prefix with "+"
    rhs_str <- paste(deparse(rhs), collapse = "")
    formula_str <- paste("+", rhs_str)
    
    # base_vars: all *simple* variable names referenced in the formula
    # (excludes interaction notation like ":" — just the column names)
    all_vars <- all.vars(covariates)
    return(list(formula_str = formula_str, base_vars = all_vars))
  }
  
  stop("'covariates' must be NULL, a character vector, or a one-sided formula (e.g. ~ age * sex).")
}


# -----------------------------------------------------------------------------
# Helper: prune_covariate_formula
#
# Removes from a covariate formula string all terms that contain a given
# variable name (the interaction_by stratification variable).
#
# This is needed for per-stratum models: within a stratum, interaction_by
# is constant, so any term that involves it (main effect OR interaction with
# another covariate such as age:sex) is either undefined or perfectly
# collinear with the intercept.  Keeping such terms causes rank deficiency.
#
# Works on the formula_str representation (e.g. "+ age + sex + age:sex").
# Splits on "+", drops any token whose all.vars() contains strip_var,
# then reassembles.
#
# @param formula_str  character, e.g. "+ age + sex + age:sex + age:bmi"
# @param strip_var    character scalar, e.g. "sex"
# @return             pruned formula string, e.g. "+ age + age:bmi"
#                     or "" if nothing remains
# -----------------------------------------------------------------------------
prune_covariate_formula <- function(formula_str, strip_var) {
  if (nchar(trimws(formula_str)) == 0L || is.null(strip_var)) {
    return(formula_str)
  }
  
  # Split into individual terms (strip leading "+", trim whitespace)
  raw   <- strsplit(formula_str, "\\+")[[1L]]
  terms <- trimws(raw)
  terms <- terms[nchar(terms) > 0L]
  
  # Keep a term only if strip_var does NOT appear among its variable names.
  # We parse each term as a tiny formula to leverage all.vars().
  keep_term <- vapply(terms, function(tm) {
    vars <- tryCatch(
      all.vars(stats::as.formula(paste("~", tm))),
      error = function(e) character(0)
    )
    !strip_var %in% vars
  }, logical(1L))
  
  kept <- terms[keep_term]
  if (length(kept) == 0L) return("")
  paste("+", paste(kept, collapse = " + "))
}


# -----------------------------------------------------------------------------
# Helper: Rank Inverse Normal Transformation (Blom's formula)
# -----------------------------------------------------------------------------
rint <- function(x, c = 3/8) {
  r <- rank(x, na.last = "keep", ties.method = "average")
  qnorm((r - c) / (sum(!is.na(x)) - 2 * c + 1))
}


# -----------------------------------------------------------------------------
# Helper: Z-score normalisation
# -----------------------------------------------------------------------------
z_normalise <- function(x) {
  mu <- mean(x, na.rm = TRUE)
  s  <- sd(x,   na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mu) / s
}


# -----------------------------------------------------------------------------
# Helper: assess whether RINT is needed
# -----------------------------------------------------------------------------
needs_rint <- function(x, skew_threshold = 2) {
  x <- x[!is.na(x)];  n <- length(x)
  if (n < 10) return(FALSE)
  m <- mean(x);  s <- sd(x)
  if (s == 0) return(FALSE)
  abs((sum((x - m)^3) / n) / s^3) > skew_threshold
}


# -----------------------------------------------------------------------------
# Helper: compute skewness
# -----------------------------------------------------------------------------
compute_skewness <- function(x) {
  x <- x[!is.na(x)];  n <- length(x)
  if (n < 3) return(NA_real_)
  m <- mean(x);  s <- sd(x)
  if (s == 0) return(NA_real_)
  (sum((x - m)^3) / n) / s^3
}


# -----------------------------------------------------------------------------
# Helper: build spline-term NA columns (beta_sp_i / se_sp_i)
# -----------------------------------------------------------------------------
make_sp_na_dt <- function(k) {
  sp_betas <- setNames(as.list(rep(NA_real_, k - 1L)),
                       sprintf("beta_sp_%d", seq_len(k - 1L)))
  sp_ses   <- setNames(as.list(rep(NA_real_, k - 1L)),
                       sprintf("se_sp_%d",   seq_len(k - 1L)))
  as.data.table(c(sp_betas, sp_ses))
}


# =============================================================================
# Helper: build the functional-form expression string for x
# Returns either "x"  (linear)  or  "rms::rcs(x, k)" / "splines::ns(x, df)"
# =============================================================================
x_form_str <- function(use_spline, spline_type, k) {
  if (!use_spline) return("x")
  if (spline_type == "rcs") {
    sprintf("rms::rcs(x, %d)", k)
  } else {
    sprintf("splines::ns(x, df = %d)", k - 1L)
  }
}


# =============================================================================
# Helper: extract spline basis term rows from a glm coefficient table.
# Returns a submatrix of coef_mat whose rownames relate to "x" but are not
# the intercept or any covariate.
# =============================================================================
extract_x_coefs <- function(coef_mat, covariates) {
  rn   <- rownames(coef_mat)
  excl <- c("(Intercept)", if (!is.null(covariates)) covariates else character(0))
  is_x <- grepl("\\bx\\b|rcs\\(x|ns\\(x", rn) & !rn %in% excl
  coef_mat[is_x, , drop = FALSE]
}


# =============================================================================
# Helper: determine whether a feature should be treated as categorical
# =============================================================================
is_categorical <- function(x) {
  is.factor(x) || is.character(x) || is.logical(x)
}


# =============================================================================
# Helper: prepare a categorical feature as a factor with a given reference
#
# @param x          vector (factor / character / logical)
# @param ref_level  character scalar or NULL (auto = alphabetically first)
# @return           factor with ref_level as the first (reference) level
# =============================================================================
prepare_factor <- function(x, ref_level = NULL) {
  f <- if (is.factor(x)) x else factor(x)
  lvls <- levels(f)
  
  if (is.null(ref_level)) {
    # Default: alphabetically / numerically first level
    return(f)
  }
  
  ref_level <- as.character(ref_level)
  if (!ref_level %in% lvls)
    stop(sprintf(
      "Reference level '%s' not found in feature levels: %s",
      ref_level, paste(lvls, collapse = ", ")
    ))
  
  f <- stats::relevel(f, ref = ref_level)
  f
}


# =============================================================================
# SECTION 1: OUTCOME DETECTION
# =============================================================================

# -----------------------------------------------------------------------------
# Detect whether a numeric vector looks like a count outcome.
# Heuristics: all non-negative integers, no natural upper bound (i.e. max is
# not close to n), and max > 1 (excludes binary coded as 0/1).
# -----------------------------------------------------------------------------
.looks_like_count <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) return(FALSE)
  is_int   <- all(x == floor(x))
  all_nonneg <- all(x >= 0)
  max_gt1  <- max(x) > 1L
  # Reject if it looks binary (only 0 and 1)
  not_binary <- length(unique(x)) > 2L
  # Reject if it looks proportion-bounded (max <= 1 and non-integer variety)
  is_int && all_nonneg && max_gt1 && not_binary
}

# -----------------------------------------------------------------------------
# Detect outcome type for one column.
# Returns one of: "continuous", "binary", "count", "categorical"
# Can be overridden by the user via outcome_model.
# -----------------------------------------------------------------------------
detect_outcome_type <- function(x, name = "") {
  x_nona <- x[!is.na(x)]
  n_uniq  <- length(unique(x_nona))
  
  if (is.logical(x) || (is.integer(x) && n_uniq == 2L) ||
      (is.numeric(x) && n_uniq == 2L && all(x_nona %in% c(0, 1)))) {
    return("binary")
  }
  
  if (is.factor(x) || is.character(x)) {
    if (n_uniq == 2L) return("binary")
    return("categorical")
  }
  
  if (is.numeric(x) || is.integer(x)) {
    if (n_uniq == 2L) return("binary")
    if (.looks_like_count(x)) return("count")
    return("continuous")
  }
  
  stop(sprintf("Cannot determine outcome type for column '%s'.", name))
}

# -----------------------------------------------------------------------------
# Map outcome type to a human-readable model name
# -----------------------------------------------------------------------------
model_name_for_type <- function(type) {
  switch(type,
         continuous  = "linear",
         binary      = "logistic",
         count       = "negbinom",
         categorical = "multinomial",
         stop(sprintf("Unknown outcome type: '%s'", type))
  )
}


# =============================================================================
# SECTION 2: PSEUDO R-SQUARED (Nagelkerke) FOR GLM-FAMILY MODELS
# =============================================================================
nagelkerke_r2 <- function(fit, fit_null) {
  tryCatch({
    ll_full <- as.numeric(stats::logLik(fit))
    ll_null <- as.numeric(stats::logLik(fit_null))
    n       <- nobs(fit)
    r2_cs   <- 1 - exp(-(2/n) * (ll_full - ll_null))
    r2_max  <- 1 - exp((2/n) * ll_null)
    r2_cs / r2_max
  }, error = function(e) NA_real_)
}

r_squared_lm <- function(fit) {
  tryCatch({
    ss_res <- sum(stats::residuals(fit)^2)
    ss_tot <- sum((fit$model[[1L]] - mean(fit$model[[1L]]))^2)
    if (ss_tot == 0) NA_real_ else 1 - ss_res / ss_tot
  }, error = function(e) NA_real_)
}


# =============================================================================
# SECTION 3: LRT / F-TEST WRAPPER
# Returns p-value of model comparison (fit_null vs fit_full).
# Uses F-test for lm, LRT (Chi^2) for glm/multinom/glm.nb.
# =============================================================================
# -----------------------------------------------------------------------------
# Helper: refit a negbin model with theta fixed at a supplied value.
# Used by model_lrt() to hold theta constant across the null and full model
# so that the logLik difference is chi-squared distributed.
# -----------------------------------------------------------------------------
.refit_negbin_fixed_theta <- function(fit, theta) {
  # We cannot refit via formula + fit$model because covariate terms involving
  # functions (e.g. rcs(), ns(), I()) are stored in fit$model as pre-evaluated
  # matrix columns under compound names like "rcs(distance_to_footprint, 3)".
  # Re-parsing the formula against fit$model would require the original raw
  # column (e.g. distance_to_footprint) which is not present there.
  #
  # Solution: bypass formula re-parsing entirely.  Use the already-evaluated
  # model matrix X and response y directly with glm.fit(), which accepts
  # pre-computed matrices and avoids any column-lookup issue.
  tryCatch({
    X <- stats::model.matrix(fit)          # evaluated design matrix (n x p)
    y <- fit$model[[1L]]                   # response vector (first column)
    stats::glm.fit(
      x      = X,
      y      = y,
      family = MASS::negative.binomial(theta = theta)
    )
  }, error = function(e) NULL)
}


# Helper: extract rank from a glm.fit result (which lacks a $rank slot in
# some R versions — fall back to ncol of the qr decomposition).
.negbin_fixed_rank <- function(fit_obj) {
  if (!is.null(fit_obj$rank)) return(fit_obj$rank)
  if (!is.null(fit_obj$qr))   return(fit_obj$qr$rank)
  sum(!is.na(fit_obj$coefficients))
}


model_lrt <- function(fit_null, fit_full) {
  tryCatch({
    is_lm       <- inherits(fit_full, "lm") && !inherits(fit_full, "glm")
    is_negbin   <- inherits(fit_full, "negbin")
    is_multinom <- inherits(fit_full, "multinom")
    
    if (is_lm) {
      # Linear model: exact F-test
      ft <- stats::anova(fit_null, fit_full)
      if (nrow(ft) >= 2L) ft[2L, "Pr(>F)"] else NA_real_
      
    } else if (is_negbin) {
      # Negative binomial: LRT with theta fixed at the full-model MLE.
      #
      # Naive approach — comparing logLik of two glm.nb fits where each
      # independently re-estimates theta — is incorrect: the 2*delta-logLik
      # statistic is not chi-squared distributed when theta differs between
      # models, leading to inflated (or occasionally deflated) p-values.
      #
      # Correct approach (Lawless 1987; Venables & Ripley MASS 4th ed.):
      #   1. Fit the full model freely -> get theta_hat
      #   2. Refit BOTH models via glm.fit() with pre-evaluated model matrices
      #      and family = negative.binomial(theta = theta_hat) [theta fixed].
      #      Using glm.fit() with model.matrix() avoids re-parsing the formula
      #      against fit$model, which fails when covariates contain function
      #      calls (e.g. rcs(), ns()) stored as matrix columns.
      #   3. Compute LRT: 2*delta-logLik is chi-squared with df = difference
      #      in regression parameters only (theta no longer contributes to df).
      theta_hat <- fit_full$theta
      if (is.null(theta_hat) || !is.finite(theta_hat) || theta_hat <= 0)
        return(NA_real_)
      
      fit_full_fixed <- .refit_negbin_fixed_theta(fit_full, theta_hat)
      fit_null_fixed <- .refit_negbin_fixed_theta(fit_null, theta_hat)
      
      if (is.null(fit_full_fixed) || is.null(fit_null_fixed))
        return(NA_real_)
      
      # df: number of non-NA coefficients (= rank of the design matrix).
      # glm.fit() returns $rank; fall back to counting non-NA coefficients.
      df_full <- .negbin_fixed_rank(fit_full_fixed)
      df_null <- .negbin_fixed_rank(fit_null_fixed)
      df_diff <- df_full - df_null
      
      if (df_diff <= 0L) return(NA_real_)
      
      # Log-likelihood from glm.fit: -(deviance/2 + constant).
      # glm.fit stores the deviance; recover logLik as:
      #   logLik = -deviance/2  (up to a constant that cancels in the LRT)
      # This is exact for the negative binomial family evaluated at fixed theta.
      ll_full <- -fit_full_fixed$deviance / 2
      ll_null <- -fit_null_fixed$deviance / 2
      stats::pchisq(2 * (ll_full - ll_null), df = df_diff, lower.tail = FALSE)
      
    } else if (is_multinom) {
      # multinom: no anova() method; use logLik chi-squared test.
      # theta is not involved here so the test is straightforward.
      ll_full <- as.numeric(stats::logLik(fit_full))
      ll_null <- as.numeric(stats::logLik(fit_null))
      df_diff <- attr(stats::logLik(fit_full), "df") -
        attr(stats::logLik(fit_null), "df")
      if (df_diff <= 0L) return(NA_real_)
      stats::pchisq(2 * (ll_full - ll_null), df = df_diff, lower.tail = FALSE)
      
    } else {
      # Standard glm (binomial, Poisson, Gaussian): anova LRT
      lrt <- stats::anova(fit_null, fit_full, test = "LRT")
      if (nrow(lrt) >= 2L) lrt[2L, "Pr(>Chi)"] else NA_real_
    }
  }, error = function(e) NA_real_)
}


# =============================================================================
# SECTION 4: FIT WRAPPERS
# Each returns a fitted model object or NULL on failure.
# =============================================================================

.fit_model <- function(formula, data, outcome_type, dbg) {
  tryCatch({
    switch(outcome_type,
           continuous  = stats::lm(formula, data = data),
           binary      = stats::glm(formula, data = data,
                                    family = binomial(link = "logit")),
           count       = MASS::glm.nb(formula, data = data),
           categorical = {
             # suppress multinom iteration output
             suppressMessages(
               nnet::multinom(formula, data = data, trace = FALSE)
             )
           }
    )
  }, error = function(e) {
    dbg("fit ERROR [%s]: %s", outcome_type, conditionMessage(e))
    NULL
  })
}


# =============================================================================
# SECTION 5: COEFFICIENT EXTRACTION
# Returns a unified list:
#   $coef_tab   matrix with rows = terms, cols = Estimate/SE/lower95/upper95/p
#   $effect_tab matrix with same rows, cols = effect/eff_lower/eff_upper
#                (exponentiated for non-linear link; same as coef for linear)
#   $outcome_levels  character vector of outcome levels (multinomial only)
# For multinomial, rows are named "<outcome_level>:<term>".
# =============================================================================
.extract_coefs <- function(fit, outcome_type, ci_level = 0.95) {
  alpha <- 1 - ci_level
  z95   <- qnorm(1 - alpha / 2)
  
  exponentiate <- outcome_type %in% c("binary", "count", "categorical")
  
  if (outcome_type == "categorical") {
    # nnet::multinom: coef() returns a matrix (levels x terms) or vector
    cf  <- coef(fit)
    if (is.vector(cf)) {
      # Binary fallback (shouldn't happen since we use glm for binary)
      cf <- matrix(cf, nrow = 1L,
                   dimnames = list("1", names(cf)))
    }
    # Standard errors via vcov — multinom stores them in summary
    sm  <- summary(fit)
    se_mat <- sm$standard.errors
    if (is.vector(se_mat))
      se_mat <- matrix(se_mat, nrow = 1L,
                       dimnames = list(rownames(cf), names(se_mat)))
    
    out_levels <- rownames(cf)
    term_names <- colnames(cf)
    # Exclude intercept
    keep <- term_names != "(Intercept)"
    cf     <- cf[, keep, drop = FALSE]
    se_mat <- se_mat[, keep, drop = FALSE]
    
    # Build unified row names: "<level>:<term>"
    rn <- as.vector(outer(out_levels, colnames(cf),
                          function(a, b) paste0(a, ":", b)))
    beta_vec <- as.vector(cf)
    se_vec   <- as.vector(se_mat)
    p_vec    <- 2 * pnorm(-abs(beta_vec / se_vec))
    lo_vec   <- beta_vec - z95 * se_vec
    hi_vec   <- beta_vec + z95 * se_vec
    
    coef_tab <- matrix(c(beta_vec, se_vec, lo_vec, hi_vec, p_vec),
                       ncol = 5L,
                       dimnames = list(rn,
                                       c("Estimate","Std. Error","lower95","upper95","p")))
    
    eff_tab  <- coef_tab
    eff_tab[, c("Estimate","lower95","upper95")] <-
      exp(coef_tab[, c("Estimate","lower95","upper95")])
    
    return(list(coef_tab = coef_tab, effect_tab = eff_tab,
                outcome_levels = out_levels))
  }
  
  # linear / logistic / negbinom
  sc  <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
  if (is.null(sc)) return(NULL)
  
  # CI
  ci <- tryCatch(
    stats::confint(fit, level = ci_level),
    error = function(e) {
      # fallback: Wald CI
      b  <- sc[, "Estimate"]
      se <- sc[, if ("Std. Error" %in% colnames(sc)) "Std. Error" else 2L]
      cbind("2.5 %" = b - z95 * se, "97.5 %" = b + z95 * se)
    }
  )
  
  p_col <- if ("Pr(>|t|)" %in% colnames(sc)) "Pr(>|t|)" else "Pr(>|z|)"
  se_col <- if ("Std. Error" %in% colnames(sc)) "Std. Error" else 2L
  
  coef_tab <- matrix(
    c(sc[, "Estimate"], sc[, se_col],
      ci[, 1L], ci[, 2L], sc[, p_col]),
    ncol = 5L,
    dimnames = list(rownames(sc),
                    c("Estimate","Std. Error","lower95","upper95","p"))
  )
  
  eff_tab <- coef_tab
  if (exponentiate)
    eff_tab[, c("Estimate","lower95","upper95")] <-
    exp(coef_tab[, c("Estimate","lower95","upper95")])
  
  list(coef_tab = coef_tab, effect_tab = eff_tab, outcome_levels = NULL)
}


# =============================================================================
# SECTION 6: WIDE-TO-LONG CONVERSION
# Converts one wide per-feature result row into long-format term rows.
# =============================================================================
..wide_to_long_unified <- function(wide_dt, knots, regression_model) {
  k <- as.integer(knots)
  
  ia_cols  <- grep("^ia_", names(wide_dt), value = TRUE)
  # adj_cols excluded: adjustment runs post-stack in regression_analysis()
  
  meta_cols <- intersect(
    c("feature", "feature_type", "n_total", "n_missing",
      "skewness_raw", "rint_applied_feature", "outcome_rint",
      "normalised", "spline_knots", "ref_level",
      "aic_linear", "aic_nonlinear", "r2_linear", "r2_nonlinear",
      "p_overall_spline", "p_nonlinear", "nonlinear_df",
      "preferred_model", "error", ia_cols),
    names(wide_dt)
  )
  
  rows <- list()
  
  for (i in seq_len(nrow(wide_dt))) {
    r    <- wide_dt[i]
    ftyp <- r$feature_type
    meta <- r[, meta_cols, with = FALSE]
    
    # Unified p_overall_association
    p_oa <- if (ftyp == "numeric") {
      if ("p_linear" %in% names(r)) r$p_linear else NA_real_
    } else {
      if ("p_overall_cat" %in% names(r)) r$p_overall_cat else NA_real_
    }
    # p_overall_adj filled in by regression_analysis() post-stack
    meta[, p_overall_association := p_oa]
    meta[, p_overall_adj         := NA_real_]
    
    # Helper to build one term row
    make_row <- function(term, term_type, outcome_level = NA_character_,
                         beta = NA_real_, se = NA_real_,
                         lo = NA_real_, hi = NA_real_,
                         effect = NA_real_, eff_lo = NA_real_, eff_hi = NA_real_,
                         pv = NA_real_) {
      cbind(
        data.table(
          term              = term,
          term_type         = term_type,
          outcome_level     = outcome_level,
          regression_model  = regression_model,
          beta              = beta,
          se                = se,
          lower95           = lo,
          upper95           = hi,
          effect            = effect,
          effect_lower95    = eff_lo,
          effect_upper95    = eff_hi,
          p_value           = pv
        ),
        meta
      )
    }
    
    # --- OVERALL row ---
    rows <- c(rows, list(make_row("overall", "overall", pv = p_oa)))
    
    if (ftyp == "numeric") {
      
      # --- LINEAR term ---
      b  <- if ("beta_lin"     %in% names(r)) r$beta_lin     else NA_real_
      se <- if ("se_lin"       %in% names(r)) r$se_lin       else NA_real_
      lo <- if ("beta_lower95" %in% names(r)) r$beta_lower95 else NA_real_
      hi <- if ("beta_upper95" %in% names(r)) r$beta_upper95 else NA_real_
      pv <- if ("p_linear"     %in% names(r)) r$p_linear     else NA_real_
      
      eff    <- if ("effect_lin"       %in% names(r)) r$effect_lin       else b
      eff_lo <- if ("effect_lower95"   %in% names(r)) r$effect_lower95   else lo
      eff_hi <- if ("effect_upper95"   %in% names(r)) r$effect_upper95   else hi
      
      rows <- c(rows, list(make_row(
        "linear", "linear",
        beta = b, se = se, lo = lo, hi = hi,
        effect = eff, eff_lo = eff_lo, eff_hi = eff_hi,
        pv = pv
      )))
      
      # --- SPLINE terms ---
      for (si in seq_len(k - 1L)) {
        b_col <- sprintf("beta_sp_%d", si)
        s_col <- sprintf("se_sp_%d", si)
        rows <- c(rows, list(make_row(
          sprintf("sp_%d", si), "spline",
          beta = if (b_col %in% names(r)) r[[b_col]] else NA_real_,
          se   = if (s_col %in% names(r)) r[[s_col]] else NA_real_
        )))
      }
      
    } else {
      # --- CATEGORICAL level terms ---
      # Discover non-ref levels from cat_*_beta columns
      cat_beta_cols <- grep("^cat_.*_beta$", names(r), value = TRUE)
      cat_safelvls  <- sub("_beta$", "", sub("^cat_", "", cat_beta_cols))
      
      # Multinomial: columns are cat_<outlvl>_<safelvl>_beta
      # Regular categorical: columns are cat_<safelvl>_beta
      for (sl in cat_safelvls) {
        b_col  <- paste0("cat_", sl, "_beta")
        se_col <- paste0("cat_", sl, "_se")
        lo_col <- paste0("cat_", sl, "_lower95")
        hi_col <- paste0("cat_", sl, "_upper95")
        p_col  <- paste0("cat_", sl, "_p")
        ef_col <- paste0("cat_", sl, "_effect")
        elo_col <- paste0("cat_", sl, "_effect_lower95")
        ehi_col <- paste0("cat_", sl, "_effect_upper95")
        ol_col  <- paste0("cat_", sl, "_outcome_level")
        
        b  <- if (b_col  %in% names(r)) r[[b_col]]  else NA_real_
        se <- if (se_col %in% names(r)) r[[se_col]] else NA_real_
        lo <- if (lo_col %in% names(r)) r[[lo_col]] else NA_real_
        hi <- if (hi_col %in% names(r)) r[[hi_col]] else NA_real_
        pv <- if (p_col  %in% names(r)) r[[p_col]]  else NA_real_
        eff    <- if (ef_col  %in% names(r)) r[[ef_col]]  else NA_real_
        eff_lo <- if (elo_col %in% names(r)) r[[elo_col]] else NA_real_
        eff_hi <- if (ehi_col %in% names(r)) r[[ehi_col]] else NA_real_
        out_lv <- if (ol_col  %in% names(r)) r[[ol_col]]  else NA_character_
        
        rows <- c(rows, list(make_row(
          sl, "categorical_level",
          outcome_level = out_lv,
          beta = b, se = se, lo = lo, hi = hi,
          effect = eff, eff_lo = eff_lo, eff_hi = eff_hi,
          pv = pv
        )))
      }
    }
  }
  
  out <- data.table::rbindlist(rows, fill = TRUE)
  
  key_first <- c("feature", "feature_type", "term", "term_type",
                 "outcome_level", "regression_model",
                 "beta", "se", "lower95", "upper95",
                 "effect", "effect_lower95", "effect_upper95",
                 "p_value", "p_overall_association", "p_overall_adj")
  key_first <- intersect(key_first, names(out))
  data.table::setcolorder(out, c(key_first, setdiff(names(out), key_first)))
  out[]
}


# =============================================================================
# SECTION 7: CORE PER-FEATURE ANALYSIS (unified)
# =============================================================================
analyse_one_feature_unified <- function(feature_name,
                                        dt,
                                        outcome,
                                        y_proc,
                                        outcome_type,
                                        regression_model,
                                        ref_level_outcome  = NULL,
                                        covariates         = NULL,
                                        spline_type        = c("rcs","ns"),
                                        knots              = 4L,
                                        apply_rint_feature = c("auto","always","never"),
                                        skew_threshold     = 2,
                                        normalise          = TRUE,
                                        outcome_rint       = FALSE,
                                        ref_levels         = NULL,
                                        interaction_by     = NULL,
                                        nonlinear_threshold = 0.05,
                                        verbose            = FALSE) {
  
  spline_type        <- match.arg(spline_type)
  apply_rint_feature <- match.arg(apply_rint_feature)
  k                  <- as.integer(knots)
  
  dbg <- function(...) if (verbose) message("    [DBG] ", sprintf(...))
  exponentiate <- outcome_type %in% c("binary", "count", "categorical")
  
  ia_active <- !is.null(interaction_by)
  ia_levels <- if (ia_active) {
    sort(unique(as.character(
      dt[[interaction_by]][!is.na(dt[[interaction_by]])]
    )))
  } else NULL
  
  x_raw_full  <- dt[[feature_name]]
  cat_feature <- is_categorical(x_raw_full)
  dbg("feature type: %s  outcome type: %s", 
      if (cat_feature) "categorical" else "numeric", outcome_type)
  
  cov_parsed    <- parse_covariates(covariates)
  cov_form      <- cov_parsed$formula_str
  cov_base_vars <- cov_parsed$base_vars
  
  keep_cols <- c(feature_name, outcome, cov_base_vars)
  keep      <- complete.cases(dt[, keep_cols, with = FALSE])
  
  x_raw    <- x_raw_full[keep]
  y_use    <- y_proc[keep]
  n_total  <- sum(keep)
  n_missing <- sum(!keep)
  
  dbg("data: n=%d  missing=%d", n_total, n_missing)
  
  min_n <- if (outcome_type == "categorical") 20L else 10L
  if (n_total < min_n) {
    return(.na_row_unified(
      feature_name, k, cat_feature, normalise, outcome_rint,
      "insufficient data", "Too few observations",
      n_total = n_total, n_missing = n_missing,
      ia_levels = ia_levels, regression_model = regression_model
    ))
  }
  
  cov_df <- if (length(cov_base_vars) > 0L) {
    as.data.frame(dt[keep, cov_base_vars, with = FALSE])
  } else {
    data.frame(row.names = seq_len(n_total))
  }
  
  # ===========================================================================
  # PATH A: CATEGORICAL FEATURE
  # ===========================================================================
  if (cat_feature) {
    feat_ref  <- if (!is.null(ref_levels) && feature_name %in% names(ref_levels))
      ref_levels[[feature_name]] else NULL
    
    x_factor <- tryCatch(
      prepare_factor(x_raw, ref_level = feat_ref),
      error = function(e) { dbg("factor ERROR: %s", conditionMessage(e)); NULL }
    )
    if (is.null(x_factor)) {
      return(.na_row_unified(
        feature_name, k, TRUE, normalise, outcome_rint,
        "error", "Failed to prepare factor",
        n_total = n_total, n_missing = n_missing,
        ia_levels = ia_levels, regression_model = regression_model
      ))
    }
    
    lvls    <- levels(x_factor)
    ref_lvl <- lvls[1L]
    non_ref <- lvls[-1L]
    dbg("categorical: ref='%s'  levels=%s", ref_lvl, paste(lvls, collapse=", "))
    
    adf <- cbind(data.frame(outcome = y_use, x = x_factor), cov_df)
    
    f_null <- if (nchar(trimws(cov_form)) > 0L)
      stats::as.formula(paste("outcome ~", substring(cov_form, 3)))
    else stats::as.formula("outcome ~ 1")
    
    fit_null <- .fit_model(f_null, adf, outcome_type, dbg)
    f_cat    <- stats::as.formula(paste("outcome ~ x", cov_form))
    dbg("cat formula: %s", deparse(f_cat))
    fit_cat  <- .fit_model(f_cat, adf, outcome_type, dbg)
    
    if (is.null(fit_cat)) {
      return(.na_row_unified(
        feature_name, k, TRUE, normalise, outcome_rint,
        "error", "Categorical model failed",
        n_total = n_total, n_missing = n_missing,
        ia_levels = ia_levels, regression_model = regression_model
      ))
    }
    
    p_overall_cat <- if (!is.null(fit_null)) model_lrt(fit_null, fit_cat) else NA_real_
    dbg("cat overall p=%.4g", ifelse(is.na(p_overall_cat), -1, p_overall_cat))
    
    aic_cat <- tryCatch(stats::AIC(fit_cat), error = function(e) NA_real_)
    r2_cat  <- if (outcome_type == "continuous") r_squared_lm(fit_cat) else
      if (!is.null(fit_null)) nagelkerke_r2(fit_cat, fit_null) else NA_real_
    
    ce <- .extract_coefs(fit_cat, outcome_type)
    if (is.null(ce)) {
      return(.na_row_unified(
        feature_name, k, TRUE, normalise, outcome_rint,
        "error", "Coefficient extraction failed",
        n_total = n_total, n_missing = n_missing,
        ia_levels = ia_levels, regression_model = regression_model
      ))
    }
    
    cat_cols <- list()
    
    if (outcome_type == "categorical") {
      # multinomial: ce$coef_tab rows are "<outlvl>:<term>"
      # Terms matching "^x" are the x-related coefficients
      out_levels <- ce$outcome_levels
      for (ol in out_levels) {
        ol_safe <- make.names(ol)
        # x-related rows for this outcome level
        rn_pat <- paste0("^", ol, ":x")
        rn_match <- grep(rn_pat, rownames(ce$coef_tab), value = TRUE)
        for (rn in rn_match) {
          # extract the feature level from row name: "<outlvl>:x<featlvl>"
          feat_lvl      <- sub(paste0("^", ol, ":x"), "", rn)
          feat_lvl_safe <- make.names(feat_lvl)
          col_pfx <- paste0("cat_", ol_safe, "_", feat_lvl_safe)
          cat_cols[[paste0(col_pfx, "_beta")]]           <- ce$coef_tab[rn, "Estimate"]
          cat_cols[[paste0(col_pfx, "_se")]]             <- ce$coef_tab[rn, "Std. Error"]
          cat_cols[[paste0(col_pfx, "_lower95")]]        <- ce$coef_tab[rn, "lower95"]
          cat_cols[[paste0(col_pfx, "_upper95")]]        <- ce$coef_tab[rn, "upper95"]
          cat_cols[[paste0(col_pfx, "_p")]]              <- ce$coef_tab[rn, "p"]
          cat_cols[[paste0(col_pfx, "_effect")]]         <- ce$effect_tab[rn, "Estimate"]
          cat_cols[[paste0(col_pfx, "_effect_lower95")]] <- ce$effect_tab[rn, "lower95"]
          cat_cols[[paste0(col_pfx, "_effect_upper95")]] <- ce$effect_tab[rn, "upper95"]
          cat_cols[[paste0(col_pfx, "_outcome_level")]]  <- ol
        }
      }
    } else {
      # binary / continuous / count categorical feature
      ct <- ce$coef_tab
      et <- ce$effect_tab
      for (lvl in non_ref) {
        safe <- make.names(lvl)
        rn   <- paste0("x", lvl)
        if (rn %in% rownames(ct)) {
          cat_cols[[paste0("cat_", safe, "_beta")]]           <- ct[rn, "Estimate"]
          cat_cols[[paste0("cat_", safe, "_se")]]             <- ct[rn, "Std. Error"]
          cat_cols[[paste0("cat_", safe, "_lower95")]]        <- ct[rn, "lower95"]
          cat_cols[[paste0("cat_", safe, "_upper95")]]        <- ct[rn, "upper95"]
          cat_cols[[paste0("cat_", safe, "_p")]]              <- ct[rn, "p"]
          cat_cols[[paste0("cat_", safe, "_effect")]]         <- et[rn, "Estimate"]
          cat_cols[[paste0("cat_", safe, "_effect_lower95")]] <- et[rn, "lower95"]
          cat_cols[[paste0("cat_", safe, "_effect_upper95")]] <- et[rn, "upper95"]
          cat_cols[[paste0("cat_", safe, "_outcome_level")]]  <- NA_character_
        } else {
          for (sfx in c("_beta","_se","_lower95","_upper95","_p",
                        "_effect","_effect_lower95","_effect_upper95")) {
            cat_cols[[paste0("cat_", safe, sfx)]] <- NA_real_
          }
          cat_cols[[paste0("cat_", safe, "_outcome_level")]] <- NA_character_
        }
      }
    }
    
    ia_dt <- if (ia_active) {
      .run_interaction_unified(
        adf = adf, dt = dt, keep = keep,
        x_str = "x", cov_form = cov_form,
        ia_levels = ia_levels, interaction_by = interaction_by,
        term_names = "lin", use_spline = FALSE,
        outcome_type = outcome_type, spline_type = spline_type, k = k,
        nonlinear_threshold = nonlinear_threshold, dbg = dbg
      )
    } else NULL
    
    return(cbind(
      data.table(
        feature              = feature_name,
        feature_type         = "categorical",
        n_total              = as.integer(n_total),
        n_missing            = as.integer(n_missing),
        skewness_raw         = NA_real_,
        rint_applied_feature = FALSE,
        outcome_rint         = outcome_rint,
        normalised           = normalise,
        spline_knots         = NA_integer_,
        ref_level            = ref_lvl,
        beta_lin             = NA_real_,
        se_lin               = NA_real_,
        beta_lower95         = NA_real_,
        beta_upper95         = NA_real_,
        effect_lin           = NA_real_,
        effect_lower95       = NA_real_,
        effect_upper95       = NA_real_,
        p_linear             = NA_real_,
        p_overall_cat        = p_overall_cat,
        aic_linear           = aic_cat,
        aic_nonlinear        = NA_real_,
        r2_linear            = r2_cat,
        r2_nonlinear         = NA_real_,
        p_overall_spline     = NA_real_,
        p_nonlinear          = NA_real_,
        nonlinear_df         = NA_integer_,
        preferred_model      = "categorical",
        error                = NA_character_
      ),
      as.data.table(cat_cols),
      ia_dt
    ))
  }
  
  # ===========================================================================
  # PATH B: NUMERIC FEATURE
  # ===========================================================================
  skewness_raw <- compute_skewness(as.numeric(x_raw))
  do_rint <- switch(apply_rint_feature,
                    "always" = TRUE,  "never" = FALSE,
                    "auto"   = needs_rint(as.numeric(x_raw), skew_threshold)
  )
  x_proc <- if (do_rint) rint(as.numeric(x_raw)) else as.numeric(x_raw)
  if (normalise) x_proc <- z_normalise(x_proc)
  dbg("transform: RINT=%s  normalise=%s", do_rint, normalise)
  
  adf <- cbind(data.frame(outcome = y_use, x = x_proc), cov_df)
  
  f_null <- if (nchar(trimws(cov_form)) > 0L)
    stats::as.formula(paste("outcome ~", substring(cov_form, 3)))
  else stats::as.formula("outcome ~ 1")
  fit_null <- .fit_model(f_null, adf, outcome_type, dbg)
  
  # ---- Linear / additive model ----
  f_linear   <- stats::as.formula(paste("outcome ~ x", cov_form))
  dbg("linear formula: %s", deparse(f_linear))
  fit_linear <- .fit_model(f_linear, adf, outcome_type, dbg)
  
  if (is.null(fit_linear)) {
    return(.na_row_unified(
      feature_name, k, FALSE, normalise, outcome_rint,
      "error", "Linear model failed",
      n_total = n_total, n_missing = n_missing,
      skewness = skewness_raw, do_rint = do_rint,
      ia_levels = ia_levels, regression_model = regression_model
    ))
  }
  
  ce_lin  <- .extract_coefs(fit_linear, outcome_type)
  ct      <- ce_lin$coef_tab
  et      <- ce_lin$effect_tab
  
  x_rows <- grep("^x$|^x\\b", rownames(ct), value = TRUE)
  x_row  <- if (length(x_rows) > 0L) x_rows[1L] else "x"
  
  beta_lin  <- if (x_row %in% rownames(ct)) ct[x_row, "Estimate"]   else NA_real_
  se_lin    <- if (x_row %in% rownames(ct)) ct[x_row, "Std. Error"] else NA_real_
  lo_lin    <- if (x_row %in% rownames(ct)) ct[x_row, "lower95"]    else NA_real_
  hi_lin    <- if (x_row %in% rownames(ct)) ct[x_row, "upper95"]    else NA_real_
  p_lin     <- if (x_row %in% rownames(ct)) ct[x_row, "p"]          else NA_real_
  eff_lin   <- if (x_row %in% rownames(et)) et[x_row, "Estimate"]   else NA_real_
  eff_lo    <- if (x_row %in% rownames(et)) et[x_row, "lower95"]    else NA_real_
  eff_hi    <- if (x_row %in% rownames(et)) et[x_row, "upper95"]    else NA_real_
  
  aic_lin <- tryCatch(stats::AIC(fit_linear), error = function(e) NA_real_)
  r2_lin  <- if (outcome_type == "continuous") r_squared_lm(fit_linear) else
    if (!is.null(fit_null)) nagelkerke_r2(fit_linear, fit_null) else NA_real_
  dbg("linear OK: beta=%.4f  p=%.4g  AIC=%.1f", beta_lin, p_lin, aic_lin)
  
  # ---- Spline model ----
  fit_spline <- tryCatch({
    f_sp <- stats::as.formula(
      paste("outcome ~", x_form_str(TRUE, spline_type, k), cov_form)
    )
    dbg("spline formula: %s", deparse(f_sp))
    .fit_model(f_sp, adf, outcome_type, dbg)
  }, error = function(e) NULL)
  
  spline_failed <- is.null(fit_spline)
  sp_dt <- if (spline_failed) make_sp_na_dt(k) else NULL
  dbg("spline: %s", if (spline_failed) "FAILED" else "OK")
  
  if (!spline_failed) {
    ce_sp    <- .extract_coefs(fit_spline, outcome_type)
    sp_coefs <- extract_x_coefs(ce_sp$coef_tab, cov_base_vars)
    n_sp     <- nrow(sp_coefs)
    sp_betas <- setNames(
      lapply(seq_len(k-1L), function(i)
        if (i <= n_sp) sp_coefs[i,"Estimate"]   else NA_real_),
      sprintf("beta_sp_%d", seq_len(k-1L)))
    sp_ses <- setNames(
      lapply(seq_len(k-1L), function(i)
        if (i <= n_sp) sp_coefs[i,"Std. Error"] else NA_real_),
      sprintf("se_sp_%d", seq_len(k-1L)))
    sp_dt <- as.data.table(c(sp_betas, sp_ses))
  }
  
  p_overall_sp <- if (!spline_failed && !is.null(fit_null))
    model_lrt(fit_null, fit_spline) else NA_real_
  nonlinear_df <- as.integer(k - 2L)
  p_nonlin <- if (!spline_failed)
    model_lrt(fit_linear, fit_spline) else NA_real_
  aic_sp  <- if (!spline_failed)
    tryCatch(stats::AIC(fit_spline), error = function(e) NA_real_) else NA_real_
  r2_sp   <- if (!spline_failed && !is.null(fit_null))
    (if (outcome_type == "continuous") r_squared_lm(fit_spline) else
      nagelkerke_r2(fit_spline, fit_null)) else NA_real_
  
  if (!spline_failed)
    dbg("spline OK: p_nl=%.4g  AIC=%.1f", ifelse(is.na(p_nonlin),1,p_nonlin), aic_sp)
  
  # ---- Model preference ----
  main_error <- if (spline_failed) "Spline model failed to fit" else NA_character_
  preferred <- if (spline_failed) {
    "linear (spline failed)"
  } else {
    spline_overall_ok <- !is.na(p_overall_sp) && !is.na(p_lin) &&
      p_overall_sp <= p_lin
    pref <- if (!is.na(p_nonlin) && p_nonlin < 0.05 && spline_overall_ok) {
      "non-linear (spline)"
    } else if (!is.na(p_nonlin) && p_nonlin < 0.05 && !spline_overall_ok) {
      "linear (spline NL-test p<0.05 but overall spline fit weaker)"
    } else "linear"
    if (!is.na(p_nonlin) && p_nonlin >= 0.05 && p_nonlin < 0.10 &&
        !is.na(aic_sp) && !is.na(aic_lin) && aic_sp < (aic_lin - 2) &&
        spline_overall_ok)
      pref <- "non-linear (borderline, AIC-supported)"
    pref
  }
  dbg("preferred: '%s'", preferred)
  
  # ---- Interaction ----
  ia_dt <- if (ia_active) {
    use_sp_ia  <- grepl("non-linear", preferred) && !is.na(p_nonlin) &&
      p_nonlin < nonlinear_threshold
    x_str_ia   <- x_form_str(use_sp_ia, spline_type, k)
    term_names <- if (use_sp_ia) paste0("sp", seq_len(k-1L)) else "lin"
    .run_interaction_unified(
      adf = adf, dt = dt, keep = keep,
      x_str = x_str_ia, cov_form = cov_form,
      ia_levels = ia_levels, interaction_by = interaction_by,
      term_names = term_names, use_spline = use_sp_ia,
      outcome_type = outcome_type, spline_type = spline_type, k = k,
      nonlinear_threshold = nonlinear_threshold, dbg = dbg
    )
  } else NULL
  
  cbind(
    data.table(
      feature              = feature_name,
      feature_type         = "numeric",
      n_total              = as.integer(n_total),
      n_missing            = as.integer(n_missing),
      skewness_raw         = skewness_raw,
      rint_applied_feature = do_rint,
      outcome_rint         = outcome_rint,
      normalised           = normalise,
      spline_knots         = k,
      ref_level            = NA_character_,
      beta_lin             = beta_lin,
      se_lin               = se_lin,
      beta_lower95         = lo_lin,
      beta_upper95         = hi_lin,
      effect_lin           = eff_lin,
      effect_lower95       = eff_lo,
      effect_upper95       = eff_hi,
      p_linear             = p_lin,
      p_overall_cat        = NA_real_,
      aic_linear           = aic_lin,
      aic_nonlinear        = aic_sp,
      r2_linear            = r2_lin,
      r2_nonlinear         = r2_sp,
      p_overall_spline     = p_overall_sp,
      p_nonlinear          = p_nonlin,
      nonlinear_df         = nonlinear_df,
      preferred_model      = preferred,
      error                = main_error
    ),
    sp_dt,
    ia_dt
  )
}


# =============================================================================
# SECTION 8: INTERACTION HELPER (unified)
# =============================================================================
.run_interaction_unified <- function(adf, dt, keep, x_str, cov_form,
                                     ia_levels, interaction_by,
                                     term_names, use_spline,
                                     outcome_type, spline_type, k,
                                     nonlinear_threshold,
                                     dbg = function(...) NULL) {
  
  ia_vec  <- as.character(dt[[interaction_by]])[keep]
  ia_keep <- !is.na(ia_vec)
  adf_ia  <- adf[ia_keep, , drop = FALSE]
  adf_ia$strata <- factor(ia_vec[ia_keep])
  levels_present <- levels(adf_ia$strata)
  dbg("IA strata: %s", paste(levels_present, collapse=", "))
  
  ia_form_used <- if (use_spline) paste0(spline_type,"(x,k=",k,")") else "linear"
  
  make_na_ia <- function() {
    cols <- list(ia_form_used = ia_form_used, ia_p_interaction = NA_real_)
    for (lvl in ia_levels) {
      safe <- make.names(lvl)
      cols[[paste0("ia_", safe, "_n_total")]] <- NA_integer_
      for (tm in term_names) {
        cols[[paste0("ia_", safe, "_beta_", tm)]] <- NA_real_
        cols[[paste0("ia_", safe, "_se_",   tm)]] <- NA_real_
      }
      cols[[paste0("ia_", safe, "_error")]] <- NA_character_
    }
    as.data.table(cols)
  }
  
  if (length(levels_present) < 2L) return(make_na_ia())
  
  f_add <- stats::as.formula(paste("outcome ~", x_str, "+ strata", cov_form))
  f_int <- stats::as.formula(paste("outcome ~", x_str, "* strata", cov_form))
  dbg("IA add: %s", deparse(f_add))
  fit_add <- .fit_model(f_add, adf_ia, outcome_type, dbg)
  fit_int <- .fit_model(f_int, adf_ia, outcome_type, dbg)
  
  p_ia <- if (!is.null(fit_add) && !is.null(fit_int))
    model_lrt(fit_add, fit_int) else NA_real_
  dbg("IA p=%.4g", ifelse(is.na(p_ia), -1, p_ia))
  
  stratum_cols <- list(ia_form_used = ia_form_used, ia_p_interaction = p_ia)
  
  for (lvl in ia_levels) {
    safe  <- make.names(lvl)
    idx   <- adf_ia$strata == lvl
    n_lvl <- sum(idx)
    stratum_cols[[paste0("ia_", safe, "_n_total")]] <- as.integer(n_lvl)
    
    if (lvl %in% levels_present && n_lvl >= 10L) {
      sub_df       <- adf_ia[idx, , drop = FALSE]
      cov_form_sub <- prune_covariate_formula(cov_form, interaction_by)
      f_sub        <- stats::as.formula(paste("outcome ~", x_str, cov_form_sub))
      dbg("  stratum '%s' cov pruned '%s'->'%s'", lvl, cov_form, cov_form_sub)
      fit_sub      <- .fit_model(f_sub, sub_df, outcome_type, dbg)
      
      if (!is.null(fit_sub)) {
        ce_sub  <- .extract_coefs(fit_sub, outcome_type)
        sub_ct  <- if (!is.null(ce_sub)) ce_sub$coef_tab else NULL
        sub_cv  <- if (!is.null(cov_form_sub) && nchar(trimws(cov_form_sub)) > 0L)
          all.vars(stats::as.formula(paste("~", substring(cov_form_sub, 3))))
        else character(0)
        x_coefs <- if (!is.null(sub_ct)) extract_x_coefs(sub_ct, sub_cv) else NULL
        
        for (ti in seq_along(term_names)) {
          tm <- term_names[ti]
          if (!is.null(x_coefs) && ti <= nrow(x_coefs)) {
            stratum_cols[[paste0("ia_", safe, "_beta_", tm)]] <- x_coefs[ti,"Estimate"]
            stratum_cols[[paste0("ia_", safe, "_se_",   tm)]] <- x_coefs[ti,"Std. Error"]
          } else {
            stratum_cols[[paste0("ia_", safe, "_beta_", tm)]] <- NA_real_
            stratum_cols[[paste0("ia_", safe, "_se_",   tm)]] <- NA_real_
          }
        }
        stratum_cols[[paste0("ia_", safe, "_error")]] <- NA_character_
      } else {
        for (tm in term_names) {
          stratum_cols[[paste0("ia_", safe, "_beta_", tm)]] <- NA_real_
          stratum_cols[[paste0("ia_", safe, "_se_",   tm)]] <- NA_real_
        }
        stratum_cols[[paste0("ia_", safe, "_error")]] <- "stratum model failed"
      }
    } else {
      skip_reason <- if (!lvl %in% levels_present) "not present" else
        sprintf("n=%d < 10", n_lvl)
      for (tm in term_names) {
        stratum_cols[[paste0("ia_", safe, "_beta_", tm)]] <- NA_real_
        stratum_cols[[paste0("ia_", safe, "_se_",   tm)]] <- NA_real_
      }
      stratum_cols[[paste0("ia_", safe, "_error")]] <- skip_reason
    }
  }
  as.data.table(stratum_cols)
}


# =============================================================================
# SECTION 9: NA ROW BUILDER (unified)
# =============================================================================
.na_row_unified <- function(feature_name, k, cat_feature, normalise,
                            outcome_rint, preferred_model, error_msg,
                            n_total = NA_integer_, n_missing = NA_integer_,
                            skewness = NA_real_, do_rint = FALSE,
                            ia_levels = NULL, regression_model = NA_character_) {
  
  sp_dt  <- if (!cat_feature) make_sp_na_dt(k) else NULL
  ia_dt  <- if (!is.null(ia_levels))
    .run_interaction_unified_na(ia_levels, "lin") else NULL
  
  base <- data.table(
    feature              = feature_name,
    feature_type         = if (cat_feature) "categorical" else "numeric",
    n_total              = as.integer(n_total),
    n_missing            = as.integer(n_missing),
    skewness_raw         = as.numeric(skewness),
    rint_applied_feature = do_rint,
    outcome_rint         = outcome_rint,
    normalised           = normalise,
    spline_knots         = if (!cat_feature) as.integer(k) else NA_integer_,
    ref_level            = NA_character_,
    beta_lin             = NA_real_,
    se_lin               = NA_real_,
    beta_lower95         = NA_real_,
    beta_upper95         = NA_real_,
    effect_lin           = NA_real_,
    effect_lower95       = NA_real_,
    effect_upper95       = NA_real_,
    p_linear             = NA_real_,
    p_overall_cat        = NA_real_,
    aic_linear           = NA_real_,
    aic_nonlinear        = NA_real_,
    r2_linear            = NA_real_,
    r2_nonlinear         = NA_real_,
    p_overall_spline     = NA_real_,
    p_nonlinear          = NA_real_,
    nonlinear_df         = NA_integer_,
    preferred_model      = preferred_model,
    error                = error_msg
  )
  cbind(base, sp_dt, ia_dt)
}

.run_interaction_unified_na <- function(ia_levels, term_names) {
  cols <- list(ia_form_used = NA_character_, ia_p_interaction = NA_real_)
  for (lvl in ia_levels) {
    safe <- make.names(lvl)
    cols[[paste0("ia_", safe, "_n_total")]] <- NA_integer_
    for (tm in term_names) {
      cols[[paste0("ia_", safe, "_beta_", tm)]] <- NA_real_
      cols[[paste0("ia_", safe, "_se_",   tm)]] <- NA_real_
    }
    cols[[paste0("ia_", safe, "_error")]] <- NA_character_
  }
  as.data.table(cols)
}


# =============================================================================
# SECTION 10: PER-OUTCOME LOOP
# =============================================================================
..run_one_outcome_unified <- function(dt, outcome, outcome_type,
                                      regression_model, y_proc,
                                      ref_level_outcome,
                                      features, covariates,
                                      cov_bvars, cov_spec,
                                      outcome_rint, spline_type, knots,
                                      apply_rint_feature, skew_threshold,
                                      normalise, ref_levels, p_adjust,
                                      interaction_by, nonlinear_threshold,
                                      verbose) {
  
  if (verbose)
    message(sprintf(
      "  Outcome '%s' [%s | %s]%s | %d feature(s)",
      outcome, outcome_type, regression_model,
      if (outcome_rint) " [RINT]" else "",
      length(features)
    ))
  
  ia_lvls <- if (!is.null(interaction_by))
    sort(unique(as.character(dt[[interaction_by]][!is.na(dt[[interaction_by]])])))
  else NULL
  
  results <- vector("list", length(features))
  
  for (i in seq_along(features)) {
    feat <- features[i]
    if (verbose) message(sprintf("    [%d/%d] %s", i, length(features), feat))
    
    results[[i]] <- tryCatch(
      analyse_one_feature_unified(
        feature_name        = feat,
        dt                  = dt,
        outcome             = outcome,
        y_proc              = y_proc,
        outcome_type        = outcome_type,
        regression_model    = regression_model,
        ref_level_outcome   = ref_level_outcome,
        covariates          = covariates,
        spline_type         = spline_type,
        knots               = knots,
        apply_rint_feature  = apply_rint_feature,
        skew_threshold      = skew_threshold,
        normalise           = normalise,
        outcome_rint        = outcome_rint,
        ref_levels          = ref_levels,
        interaction_by      = interaction_by,
        nonlinear_threshold = nonlinear_threshold,
        verbose             = verbose
      ),
      error = function(e) {
        is_cat <- tryCatch(is_categorical(dt[[feat]]), error = function(e2) FALSE)
        .na_row_unified(feat, knots, is_cat, normalise, outcome_rint,
                        "error", conditionMessage(e),
                        ia_levels = ia_lvls,
                        regression_model = regression_model)
      }
    )
  }
  
  out <- data.table::rbindlist(results, fill = TRUE)
  
  # p-value adjustment is deferred: applied after all outcomes are
  # stacked in regression_analysis() so the full set of p-values is
  # corrected jointly across all (outcome x feature) combinations.
  
  
  # Convert to long format
  out <- ..wide_to_long_unified(out, knots = knots,
                                regression_model = regression_model)
  
  # Sort and tag
  data.table::setorder(out, p_overall_association, feature, term_type,
                       na.last = TRUE)
  out[, outcome := outcome]
  data.table::setcolorder(out, c("outcome", setdiff(names(out), "outcome")))
  out[]
}


# =============================================================================
# SECTION 11: MAIN EXPORTED FUNCTION
# =============================================================================
#
# @param dt                  data.table
# @param outcome             character scalar or vector of outcome column names
# @param features            character vector; default = all non-outcome columns
# @param covariates          NULL / character vector / one-sided formula
# @param outcome_type        optional override: named list or named character
#                            vector mapping outcome names to model type:
#                            "continuous","binary","count","categorical".
#                            E.g. list(my_var = "continuous")
#                            Unspecified outcomes are auto-detected.
# @param ref_level_outcome   named list: reference level for categorical or
#                            binary outcomes, e.g. list(my_outcome = "ctrl").
#                            For binary outcomes this sets which value = 0.
# @param outcome_rint        logical or named logical vector: apply RINT to
#                            continuous outcomes before fitting.  Only applied
#                            when regression_model == "linear".
# @param spline_type         "rcs" | "ns"
# @param knots               integer >= 3 (default 4)
# @param apply_rint_feature  "auto" | "always" | "never"
# @param skew_threshold      numeric (default 2)
# @param normalise           logical (default TRUE)
# @param ref_levels          named list: reference levels for categorical features
# @param p_adjust            p.adjust method (default "BH").
#                            Applied to p_overall_association across ALL
#                            (outcome x feature) combinations in the final
#                            stacked long table — i.e. the full set of
#                            omnibus tests is corrected jointly.
#                            Also applied to p_overall_spline, p_nonlinear,
#                            and ia_p_interaction across the same pool.
#                            Use "none" to suppress adjustment.
# @param interaction_by      character: stratification factor for IA testing
# @param nonlinear_threshold numeric (default 0.05)
# @param verbose             logical (default TRUE)
#
# @return data.table in long format with columns:
#   outcome, outcome_level, regression_model, feature, feature_type,
#   term, term_type, beta, se, lower95, upper95,
#   effect, effect_lower95, effect_upper95, p_value,
#   p_overall_association, p_overall_adj,
#   + model-level and IA columns
# =============================================================================
regression_analysis <- function(dt,
                                outcome,
                                features            = NULL,
                                covariates          = NULL,
                                outcome_type        = NULL,
                                ref_level_outcome   = NULL,
                                outcome_rint        = FALSE,
                                spline_type         = c("rcs","ns"),
                                knots               = 4L,
                                apply_rint_feature  = c("auto","always","never"),
                                skew_threshold      = 2,
                                normalise           = TRUE,
                                ref_levels          = NULL,
                                p_adjust            = "BH",
                                interaction_by      = NULL,
                                nonlinear_threshold = 0.05,
                                verbose             = TRUE) {
  
  spline_type        <- match.arg(spline_type)
  apply_rint_feature <- match.arg(apply_rint_feature)
  knots              <- as.integer(knots)
  
  if (knots < 3L) stop("'knots' must be >= 3.")
  if (!data.table::is.data.table(dt)) data.table::setDT(dt)
  
  outcomes    <- as.character(outcome)
  missing_oc  <- setdiff(outcomes, names(dt))
  if (length(missing_oc) > 0L)
    stop("Outcome column(s) not found: ", paste(missing_oc, collapse=", "))
  if (!is.null(interaction_by) && !interaction_by %in% names(dt))
    stop("'interaction_by' column not found.")
  
  # ---- Resolve outcome types ----
  oc_types <- vapply(outcomes, function(oc) {
    if (!is.null(outcome_type) && oc %in% names(outcome_type)) {
      as.character(outcome_type[[oc]])
    } else {
      detect_outcome_type(dt[[oc]], name = oc)
    }
  }, character(1L))
  
  oc_models <- vapply(oc_types, model_name_for_type, character(1L))
  
  # ---- Resolve outcome_rint per outcome (only applies to continuous) ----
  rint_vec <- if (is.null(names(outcome_rint))) {
    rep_len(as.logical(outcome_rint), length(outcomes))
  } else {
    vapply(outcomes, function(o)
      if (o %in% names(outcome_rint)) as.logical(outcome_rint[[o]]) else FALSE,
      logical(1L))
  }
  names(rint_vec) <- outcomes
  
  # ---- Validate covariates ----
  cov_spec  <- parse_covariates(covariates)
  cov_bvars <- cov_spec$base_vars
  missing_cv <- setdiff(cov_bvars, names(dt))
  if (length(missing_cv) > 0L)
    stop("Covariate column(s) not found: ", paste(missing_cv, collapse=", "))
  
  # ---- Resolve features ----
  if (is.null(features)) {
    exclude  <- c(outcomes, cov_bvars, interaction_by)
    features <- names(dt)[!names(dt) %in% exclude]
  }
  if (length(features) == 0L) stop("No features to analyse.")
  
  if (verbose) {
    message(sprintf("regression_analysis: %d outcome(s) x %d feature(s)",
                    length(outcomes), length(features)))
    for (oc in outcomes)
      message(sprintf("  %s -> %s (%s)%s", oc, oc_types[[oc]], oc_models[[oc]],
                      if (oc_types[[oc]] == "continuous" && rint_vec[[oc]])
                        " [RINT]" else ""))
    if (length(cov_bvars) > 0L)
      message(sprintf("  Covariates: %s", cov_spec$formula_str))
  }
  
  outcome_results <- vector("list", length(outcomes))
  
  for (oi in seq_along(outcomes)) {
    oc       <- outcomes[oi]
    oc_type  <- oc_types[[oc]]
    oc_model <- oc_models[[oc]]
    do_rint  <- rint_vec[[oc]] && oc_type == "continuous"
    
    ref_lv_oc <- if (!is.null(ref_level_outcome) && oc %in% names(ref_level_outcome))
      ref_level_outcome[[oc]] else NULL
    
    # Pre-process outcome
    y_raw  <- dt[[oc]]
    y_proc <- if (do_rint) {
      if (verbose) message(sprintf("    Applying RINT to '%s'", oc))
      rint(as.numeric(y_raw))
    } else if (oc_type == "binary") {
      # Ensure 0/1 numeric; apply ref level if supplied
      if (!is.null(ref_lv_oc)) {
        y_fac <- prepare_factor(y_raw, ref_level = ref_lv_oc)
        as.numeric(y_fac) - 1L
      } else {
        as.numeric(as.factor(y_raw)) - 1L
      }
    } else if (oc_type == "categorical") {
      # Factor with ref level
      if (!is.null(ref_lv_oc)) {
        prepare_factor(y_raw, ref_level = ref_lv_oc)
      } else {
        as.factor(y_raw)
      }
    } else {
      as.numeric(y_raw)
    }
    
    outcome_results[[oi]] <- tryCatch(
      ..run_one_outcome_unified(
        dt = dt, outcome = oc, outcome_type = oc_type,
        regression_model = oc_model,
        y_proc = y_proc, ref_level_outcome = ref_lv_oc,
        features = features, covariates = covariates,
        cov_bvars = cov_bvars, cov_spec = cov_spec,
        outcome_rint = do_rint, spline_type = spline_type,
        knots = knots, apply_rint_feature = apply_rint_feature,
        skew_threshold = skew_threshold, normalise = normalise,
        ref_levels = ref_levels, p_adjust = p_adjust,
        interaction_by = interaction_by,
        nonlinear_threshold = nonlinear_threshold,
        verbose = verbose
      ),
      error = function(e) {
        warning(sprintf("Outcome '%s' failed: %s", oc, conditionMessage(e)))
        NULL
      }
    )
  }
  
  out <- data.table::rbindlist(
    Filter(Negate(is.null), outcome_results), fill = TRUE
  )
  
  # Preserve outcome order
  out[, outcome := factor(outcome, levels = outcomes)]
  data.table::setorder(out, outcome, na.last = TRUE)
  out[, outcome := as.character(outcome)]
  
  # ---- Multiple testing correction on the full stacked long table ----
  # Adjustment runs across ALL (outcome x feature) combinations jointly,
  # on "overall" term rows only (one omnibus p-value per outcome x feature).
  # This correctly handles both many-features-per-outcome and
  # many-outcomes-per-feature use cases, and avoids the single-p-value
  # no-op that occurs when there is only one feature per outcome.
  if (p_adjust != "none" && "p_overall_association" %in% names(out)) {
    # Strategy: operate entirely on full-length column vectors using a
    # logical mask.  This avoids .I-based row indexing, which is unreliable
    # after sorting because .I returns positions within the *filtered* view,
    # not the full table — writing them back via out[.row, col := val] then
    # hits completely different rows.
    #
    # Pattern:
    #   1. Build a full-length NA vector
    #   2. Fill in adjusted values at the logical mask positions
    #   3. Assign the full vector to the column with set()
    # set() bypasses data.table's [i,j] scoping entirely and writes by
    # absolute column position, so no index misalignment is possible.
    
    is_overall  <- out[["term"]] == "overall"
    is_num_ov   <- is_overall & out[["feature_type"]] == "numeric"
    
    .adjust_col <- function(raw_col, mask) {
      # Returns a full-length numeric vector: adjusted where mask is TRUE,
      # NA elsewhere.  p.adjust receives only the non-NA masked values so
      # the pool size is exactly the number of tests being corrected.
      adj <- rep(NA_real_, length(raw_col))
      vals <- raw_col[mask]
      adj[mask] <- p.adjust(vals, method = p_adjust)
      adj
    }
    
    # p_overall_association -> p_overall_adj
    data.table::set(out, j = "p_overall_adj",
                    value = .adjust_col(out[["p_overall_association"]],
                                        is_overall))
    
    # p_overall_spline -> p_overall_spline_adj
    if ("p_overall_spline" %in% names(out))
      data.table::set(out, j = "p_overall_spline_adj",
                      value = .adjust_col(out[["p_overall_spline"]],
                                          is_num_ov))
    
    # p_nonlinear -> p_nonlinear_adj
    if ("p_nonlinear" %in% names(out))
      data.table::set(out, j = "p_nonlinear_adj",
                      value = .adjust_col(out[["p_nonlinear"]],
                                          is_num_ov))
    
    # ia_p_interaction -> ia_p_interaction_adj
    if (!is.null(interaction_by) && "ia_p_interaction" %in% names(out))
      data.table::set(out, j = "ia_p_interaction_adj",
                      value = .adjust_col(out[["ia_p_interaction"]],
                                          is_overall))
  }
  
  if (verbose) message("Done.")
  out[]
}


# =============================================================================
# EXAMPLE / DEMO
# =============================================================================
if (FALSE) {
  set.seed(42);  n <- 2000
  
  sim_dt <- data.table::data.table(
    age      = rnorm(n, 55, 12),
    bmi      = rnorm(n, 27, 5),
    crp      = exp(rnorm(n, 1, 1.2)),
    sex      = sample(c("M","F"), n, replace=TRUE),
    genotype = sample(c("AA","AB","BB"), n, replace=TRUE, prob=c(.5,.35,.15))
  )
  lp <- with(sim_dt, -3 + 0.03*age + 0.05*bmi + 0.2*log(crp+1))
  sim_dt[, outcome_bin  := rbinom(n, 1, plogis(lp))]
  sim_dt[, outcome_cont := 5 + 0.03*age + 0.04*bmi + rnorm(n,0,.6)]
  sim_dt[, outcome_cnt  := rpois(n, exp(0.5 + 0.02*age))]
  sim_dt[, outcome_cat  := factor(sample(c("low","mid","high"), n,
                                         replace=TRUE, prob=c(.3,.4,.3)))]
  
  # Auto-detect all four model types
  res <- regression_analysis(
    dt                 = sim_dt,
    outcome            = c("outcome_bin","outcome_cont",
                           "outcome_cnt","outcome_cat"),
    features           = c("age","bmi","crp","sex","genotype"),
    covariates         = NULL,
    outcome_rint       = c(outcome_cont = TRUE),
    apply_rint_feature = "auto",
    normalise          = TRUE,
    ref_levels         = list(sex="F", genotype="AA"),
    ref_level_outcome  = list(outcome_cat="low"),
    p_adjust           = "BH",
    verbose            = TRUE
  )
  
  # Model used per outcome
  print(unique(res[, .(outcome, regression_model)]))
  
  # Overall associations (model-level summary)
  print(res[term == "overall",
            .(outcome, regression_model, feature,
              p_overall_association, p_overall_adj,
              r2_linear, preferred_model)])
  
  # Linear/additive estimates with effects (OR/IRR for non-linear outcomes)
  print(res[term_type == "linear",
            .(outcome, regression_model, feature,
              beta, se, effect, effect_lower95, effect_upper95,
              p_value, p_overall_adj)])
  
  # Categorical feature contrasts (genotype vs AA)
  print(res[feature == "genotype",
            .(outcome, regression_model, term, outcome_level,
              beta, effect, effect_lower95, effect_upper95, p_value)])
  
  # Multinomial outcome: all level contrasts
  print(res[outcome == "outcome_cat" & term_type == "categorical_level",
            .(feature, term, outcome_level, effect, p_value)])
}