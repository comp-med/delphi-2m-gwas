run.cox <- function(dat, fol, inc, adj, lab, verbose = FALSE, hide.adj = FALSE,
                    centre.var = NULL, centre.method = c("fixed", "frailty", "strata"),
                    delphi.var = NULL) {
  
  centre.method <- match.arg(centre.method)
  
  # ── packages ────────────────────────────────────────────────────────────────
  if (centre.method == "frailty") {
    if (!requireNamespace("coxme", quietly = TRUE))
      stop("Package 'coxme' is required for centre.method = 'frailty'. Install with install.packages('coxme').")
    library(coxme)
  }
  
  registerDoMC(3)
  print(dim(dat))
  
  # ── helpers ──────────────────────────────────────────────────────────────────
  INT <- function(x) {
    qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
  }
  
  is.categorical <- function(x, max.levels = 10) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return(FALSE)
    is.factor(x) || (is.numeric(x) && all(x == floor(x)) && length(unique(x)) <= max.levels)
  }
  
  vcat <- function(...) if (verbose) cat("[VERBOSE]", ..., "\n")
  
  # ── rewrite adj string once, outside the trait loop ─────────────────────────
  # Handles the centre variable swap for all three methods.
  # Pattern is intentionally flexible: matches '+ centre' with or without
  # factor(), strata(), frailty(), or (1|centre) already present, so the
  # function is idempotent if called repeatedly.
  build.adj <- function(adj.str, sex.subset = FALSE) {
    
    if (!is.null(centre.var) && nchar(trimws(centre.var)) > 0) {
      
      # Strip any existing encoding of the centre variable so we can re-inject cleanly.
      # Matches: + factor(cv), + strata(cv), + frailty(cv,...), + (1|cv), + cv  (all optional spaces)
      strip.pattern <- paste0(
        "\\+\\s*(?:factor\\(|strata\\(|frailty\\([^)]*,?\\s*|\\(1\\s*\\|\\s*)?",
        centre.var, "\\)?(?:\\s*\\))?(?=\\s*\\+|$)"
      )
      clean.adj <- trimws(gsub(strip.pattern, "", adj.str, perl = TRUE))
      
      centre.term <- switch(
        centre.method,
        fixed   = paste0("+ factor(", centre.var, ")"),
        strata  = paste0("+ strata(", centre.var, ")"),
        frailty = paste0("+ (1 | ", centre.var, ")")
      )
      
      adj.str <- paste(clean.adj, centre.term)
    }
    
    # Drop sex term when fitting sex-stratified subset
    if (sex.subset) adj.str <- gsub("\\+\\s*sex", "", adj.str)
    
    trimws(adj.str)
  }
  
  # ── coxph / coxme dispatcher ─────────────────────────────────────────────────
  # coxme does not accept a Surv formula as a character string the same way
  # coxph does, so we parse explicitly and route by method.
  fit.model <- function(formula.str, data.sub) {
    f <- as.formula(formula.str)
    if (centre.method == "frailty") {
      coxme(f, data = data.sub)
    } else {
      coxph(f, data = data.sub, ties = "breslow")
    }
  }
  
  # ── convergence check: unified across coxph and coxme ───────────────────────
  check.convergence <- function(fit) {
    if (inherits(fit, "coxme")) {
      # coxme stores convergence info in fit$loglik; no iter.max field.
      # Treat as converged unless loglik contains NaN/NA.
      any(is.nan(unlist(fit$loglik))) || any(is.na(unlist(fit$loglik)))
    } else {
      iter.max <- if (length(fit$control$iter.max) > 0) fit$control$iter.max else Inf
      fit$iter >= iter.max
    }
  }
  
  # ── Schoenfeld residuals: coxme not supported by cox.zph ────────────────────
  compute.zph <- function(fit, trait) {
    if (inherits(fit, "coxme")) {
      vcat("  >> cox.zph not available for coxme -> NA")
      return(rep(NA_real_, 2))
    }
    tryCatch({
      res <- cox.zph(fit)$table[c(trait, "GLOBAL"), 3]
      vcat("  >> Schoenfeld residuals OK")
      res
    }, error = function(e) {
      vcat("  >> Schoenfeld residuals failed:", conditionMessage(e), "-> NA")
      rep(NA_real_, 2)
    })
  }
  
  # ── extract coefficients: unified across coxph and coxme ────────────────────
  # coxme does not have a summary()$coefficients matrix like coxph.
  # We reconstruct a comparable matrix: coef, exp(coef), se(coef), z, p.
  extract.coefs <- function(fit) {
    if (inherits(fit, "coxme")) {
      b    <- fixef(fit)                          # named numeric vector
      se   <- sqrt(diag(vcov(fit)))
      z    <- b / se
      p    <- 2 * pnorm(-abs(z))
      hr   <- exp(b)
      hr.lo <- exp(b - 1.96 * se)
      hr.hi <- exp(b + 1.96 * se)
      mat  <- cbind(
        coef          = b,
        `exp(coef)`   = hr,
        `se(coef)`    = se,
        `robust se`   = se,           # no robust SE in coxme; mirror se
        z             = z,
        `Pr(>|z|)`    = p,
        `lower .95`   = hr.lo,
        `upper .95`   = hr.hi
      )
      list(coefficients = mat, nevent = fit$n[2])
    } else {
      s <- summary(fit)
      list(coefficients = s$coefficients, nevent = s$nevent)
    }
  }
  
  # ── concordance helpers ─────────────────────────────────────────────────────
  # survival::concordance() is not defined for coxme, so frailty fits return NA.
  # cfit() fits a model from a formula string but swallows errors (-> NULL),
  # cstat() returns c(concordance, se) or c(NA, NA).
  cfit  <- function(formula.str, data.sub) tryCatch(fit.model(formula.str, data.sub),
                                                    error = function(e) NULL)
  cstat <- function(fit) {
    if (is.null(fit) || inherits(fit, "coxme")) return(c(NA_real_, NA_real_))
    cc <- tryCatch(concordance(fit), error = function(e) NULL)
    if (is.null(cc)) c(NA_real_, NA_real_) else c(cc$concordance, sqrt(cc$var))
  }
  
  # ── baseline cache (keyed by sex subset) ────────────────────────────────────
  # c.base (adj only) and c.delphi (adj + Delphi LP) do not depend on the
  # exposure, so they are fitted once per sex subset on all rows of that subset
  # and reused across the exposure loop (mirrors the cached c.base approach).
  base.cache <- list()
  
  # ── main trait loop ──────────────────────────────────────────────────────────
  res.cox <- lapply(1:nrow(lab), function(x) {
    
    trait <- lab$short_name[x]
    vcat("────────────────────────────────────────")
    vcat("Trait:", trait, "| centre.method:", centre.method)
    
    sex.subset <- lab$sex[x] != "Both"
    
    # Build adjusted adj string with correct centre encoding
    adj.mod <- build.adj(adj, sex.subset = sex.subset)
    
    if (lab$inv.transform[x]) {
      ff <- paste0("Surv(", fol, ",", inc, ") ~ INT(", trait, ") ", adj.mod)
      vcat("Formula (INT):", ff)
    } else {
      ff <- paste0("Surv(", fol, ",", inc, ") ~ ", trait, " ", adj.mod)
      vcat("Formula:", ff)
    }
    
    caught.warnings <- c()
    vcat("Fitting model [method =", centre.method, "] ...")
    
    ff <- tryCatch(
      withCallingHandlers(
        {
          data.sub <- if (sex.subset) dat[sex == lab$sex[x]] else dat
          fit      <- fit.model(ff, data.sub)
          
          if (check.convergence(fit)) {
            vcat("  >> non-convergence detected -> __no_converge__")
            stop("__no_converge__")
          }
          vcat("  >> convergence OK")
          fit
        },
        warning = function(w) {
          msg <- conditionMessage(w)
          vcat("  >> WARNING:", msg)
          if (grepl("Ran out of iterations", msg)) stop("__no_converge__")
          caught.warnings <<- c(caught.warnings, msg)
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) {
        msg <- conditionMessage(e)
        vcat("  >> ERROR:", msg)
        return(NA)
      }
    )
    
    if (length(caught.warnings) > 0) {
      vcat("Re-issuing", length(caught.warnings), "non-fatal warning(s)")
      for (w in caught.warnings) warning(w, call. = FALSE)
    }
    
    # ── results extraction ───────────────────────────────────────────────────
    # length(ff) > 1 is TRUE for fitted model objects; NA has length 1
    if (length(ff) > 1) {
      vcat("Model valid -> extracting results")
      
      data.sub <- if (sex.subset) dat[sex == lab$sex[x]] else dat
      col.vals <- data.sub[[trait]]
      
      n.intersect <- if (is.categorical(col.vals)) {
        levels.obs <- sort(unique(col.vals[!is.na(col.vals)]))
        levels.exp <- levels.obs[-1]
        vcat("  >> categorical levels:", paste(levels.obs, collapse = ","))
        sapply(levels.exp, function(lv) {
          data.sub[get(trait) == lv & get(inc) == 1, .N]
        })
      } else {
        vcat("  >> continuous -> n.intersect = NA")
        NA_integer_
      }
      
      ff.p        <- compute.zph(ff, trait)
      extracted   <- extract.coefs(ff)
      coef.mat    <- extracted$coefficients
      coef.names  <- row.names(coef.mat)
      vcat("  >> coefficient rows:", paste(coef.names, collapse = ", "))
      
      # ── incremental discrimination (Harrell's C) ──────────────────────────────
      # c.exposure       : marginal C of the exposure alone
      # c.model          : C of the full model (exposure + adj)        [= concordance(ff)]
      # c.base           : C of the adjustment set alone               (cached per sex)
      # c.delta          : c.model - c.base          (added value over adj)
      # --- only when a Delphi linear predictor is supplied via delphi.var: ---
      # c.delphi         : C of adj + Delphi LP                        (cached per sex)
      # c.model.delphi   : C of exposure + adj + Delphi LP
      # c.delta.delphi   : c.model.delphi - c.delphi  (added value BEYOND Delphi)
      # c.delphi.gain    : c.delphi - c.base          (what Delphi adds over adj)
      # Note: cached baselines use all rows of the sex subset; the full models use
      # the exposure's complete cases, so c.delta is the same mild-mismatch estimate
      # you already use. For an exactly matched delta, fit baselines on
      # dat[sex == .. & !is.na(get(trait))] instead of caching.
      adj.core  <- sub("^\\s*\\+\\s*", "", adj.mod)                 # strip leading '+'
      expo.term <- if (lab$inv.transform[x]) paste0("INT(", trait, ")") else trait
      
      base.key  <- if (sex.subset) lab$sex[x] else "all"
      if (is.null(base.cache[[base.key]])) {
        vcat("  >> caching baselines for subset:", base.key)
        c.cov <- cstat(cfit(paste0("Surv(", fol, ",", inc, ") ~ ", adj.core), data.sub))[1]
        if (!is.null(delphi.var)) {
          c.dlp <- cstat(cfit(paste0("Surv(", fol, ",", inc, ") ~ ", delphi.var, " + ", adj.core),
                              data.sub))[1]
        } else {
          c.dlp <- NA_real_
        }
        base.cache[[base.key]] <<- list(c.base = c.cov, c.delphi = c.dlp)
      }
      bc <- base.cache[[base.key]]
      
      c.exposure <- cstat(cfit(paste0("Surv(", fol, ",", inc, ") ~ ", expo.term), data.sub))[1]
      cm         <- cstat(ff)                                       # full model already fitted
      c.model    <- cm[1]; c.model.se <- cm[2]
      c.delta    <- c.model - bc$c.base
      
      if (!is.null(delphi.var)) {
        cmd <- cstat(cfit(paste0("Surv(", fol, ",", inc, ") ~ ", expo.term, " + ",
                                 delphi.var, " ", adj.mod), data.sub))
        c.model.delphi    <- cmd[1]; c.model.delphi.se <- cmd[2]
        c.delphi          <- bc$c.delphi
        c.delta.delphi    <- c.model.delphi - c.delphi
        c.delphi.gain     <- c.delphi - bc$c.base
      } else {
        c.model.delphi <- c.model.delphi.se <- NA_real_
        c.delphi <- c.delta.delphi <- c.delphi.gain <- NA_real_
      }
      vcat(sprintf("  >> C: base=%.3f model=%.3f delta=%.3f | delphi=%.3f model.delphi=%.3f delta.delphi=%.3f",
                   bc$c.base, c.model, c.delta, c.delphi, c.model.delphi, c.delta.delphi))
      
      n.intersect.vec <- if (length(n.intersect) == length(coef.names)) {
        n.intersect
      } else {
        rep_len(n.intersect, length(coef.names))
      }
      
      if (hide.adj) {
        keep.rows       <- grepl(paste0("^INT\\(", trait, "\\)|^", trait), coef.names)
        vcat("  >> hide.adj: keeping", sum(keep.rows), "of", length(coef.names), "rows")
        coef.names      <- coef.names[keep.rows]
        n.intersect.vec <- n.intersect.vec[keep.rows]
        coef.mat        <- coef.mat[keep.rows, , drop = FALSE]
      }
      
      warn.str <- if (length(caught.warnings) > 0) {
        paste(caught.warnings, collapse = "; ")
      } else {
        NA_character_
      }
      
      vcat("Returning result row(s)")
      return(data.table(
        short_name_new    = coef.names,
        coef.mat,
        c.exposure        = c.exposure,
        c.model           = c.model,
        c.model.se        = c.model.se,
        c.base            = bc$c.base,
        c.delta           = c.delta,
        c.delphi          = c.delphi,
        c.model.delphi    = c.model.delphi,
        c.model.delphi.se = c.model.delphi.se,
        c.delta.delphi    = c.delta.delphi,
        c.delphi.gain     = c.delphi.gain,
        nevent          = extracted$nevent,
        nall            = sum(!is.na(dat[, ..trait])),
        n.intersect     = n.intersect.vec,
        warnings        = warn.str,
        centre.method   = centre.method,        # audit column
        p.resid.prot    = ff.p[1],
        p.resid.overall = ff.p[2]
      ))
      
    } else {
      vcat("Model is NA -> skipping trait")
      return(NULL)
    }
  })
  
  
  res.cox <- rbindlist(res.cox, fill = TRUE)
  
  list.cols <- names(which(sapply(res.cox, is.list)))
  if (length(list.cols) > 0) {
    res.cox[, (list.cols) := lapply(.SD, function(col) {
      sapply(col, function(x) {
        if (is.null(x) || length(x) == 0) NA_character_
        else paste(x, collapse = "; ")
      })
    }), .SDcols = list.cols]
  }
  
  return(res.cox)
}