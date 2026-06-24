# =============================================================================
#  code/shared/_model_averaging.R — Satellite model averaging (BMA-style)
#
#  This module reduces *model-selection uncertainty* in the linear satellite
#  bridge by replacing a single AIC/OOS-selected specification with an
#  Akaike-weight average over the best specifications.
#
#  Two pieces:
#    select_satellite_fullsample()  full-sample (in-sample) AIC screen of all
#                                   candidate satellite designs. NO out-of-sample
#                                   evaluation — appropriate for short samples
#                                   (e.g. the EBA corporate default-rate panel).
#
#    make_bma_satellite()           builds a *mixture* Bayesian satellite that is
#                                   a drop-in replacement for the single-model
#                                   object returned by make_bayesian_satellite():
#                                   - model-selection uncertainty enters through
#                                     a mixture of posterior draws across the top
#                                     specifications (draw counts proportional to
#                                     Akaike/Schwarz weights);
#                                   - parameter uncertainty within each
#                                     specification enters through the conjugate
#                                     Jeffreys posterior (fit_bayesian_lm_jeffreys).
#                                   The mixture feeds the same GIRF -> Z -> PD
#                                   pipeline unchanged, so the full posterior of
#                                   Delta PD now also integrates over which
#                                   regressors enter the satellite.
#
#  Requires (sourced earlier): par_lapply_safe(), fit_bayesian_lm_jeffreys(),
#  make_bayesian_satellite(), parse_satellite_terms(), data.table.
# =============================================================================

# -----------------------------------------------------------------------------
# Full-sample AIC screen of all candidate satellite specifications (no OOS).
# Mirrors the screening stage of select_satellite_oos() but ranks purely by the
# in-sample AIC and returns the admissible column sets so they can be re-fit for
# model averaging.
# -----------------------------------------------------------------------------
select_satellite_fullsample <- function(model_df,
                                        max_vars = 6,
                                        p_threshold = NULL,
                                        preselect_by_aic_topn = Inf,
                                        verbose = TRUE) {
  model_df <- as.data.table(model_df)
  X_names <- setdiff(names(model_df), c("Date", "Y"))
  n <- nrow(model_df)

  y <- as.numeric(model_df$Y)
  Xfull <- as.matrix(model_df[, ..X_names])
  P <- length(X_names)

  combo_idx <- unlist(
    lapply(1:min(max_vars, P), function(k) combn(P, k, simplify = FALSE)),
    recursive = FALSE
  )

  total <- length(combo_idx)

  if (verbose) {
    message(sprintf(
      "Full-sample AIC screen: %d candidate models over %d predictors.",
      total, P
    ))
  }

  need_se <- !is.null(p_threshold)

  fast_eval <- function(cols) {
    Xs <- cbind(1, Xfull[, cols, drop = FALSE])
    fit <- .lm.fit(Xs, y)
    res <- fit$residuals
    rss <- sum(res * res)

    if (!is.finite(rss) || rss <= 0) return(NULL)

    npar <- ncol(Xs)

    if (need_se) {
      XtX_inv <- tryCatch(chol2inv(chol(crossprod(Xs))), error = function(e) NULL)
      if (is.null(XtX_inv)) return(NULL)

      sigma2 <- rss / (n - npar)
      se <- sqrt(pmax(diag(XtX_inv) * sigma2, 0))
      tval <- fit$coefficients / se
      pval <- 2 * pt(-abs(tval), df = n - npar)

      if (any(pval[-1] > p_threshold)) return(NULL)
    }

    list(cols = cols, aic = n * log(rss / n) + 2 * (npar + 1))
  }

  BATCH <- 20000L
  nb <- ceiling(total / BATCH)
  cand <- vector("list", total)

  for (b in seq_len(nb)) {
    lo <- (b - 1L) * BATCH + 1L
    hi <- min(b * BATCH, total)
    cand[lo:hi] <- par_lapply_safe(combo_idx[lo:hi], fast_eval)

    if (verbose) {
      message(sprintf("  screened %d / %d (%.0f%%)", hi, total, 100 * hi / total))
    }
  }

  cand <- Filter(Negate(is.null), cand)

  if (length(cand) == 0) {
    stop("No admissible satellite specification found (full-sample screen).")
  }

  ord <- order(vapply(cand, `[[`, numeric(1), "aic"))
  cand <- cand[ord]

  if (is.finite(preselect_by_aic_topn) && length(cand) > preselect_by_aic_topn) {
    cand <- cand[seq_len(preselect_by_aic_topn)]
  }

  best_vars <- X_names[cand[[1]]$cols]
  best_model <- lm(
    as.formula(paste("Y ~", paste(best_vars, collapse = " + "))),
    data = model_df
  )

  ranking <- rbindlist(lapply(cand, function(cc) {
    vars <- X_names[cc$cols]
    data.table(
      model = paste(vars, collapse = " + "),
      n_vars = length(vars),
      aic = cc$aic
    )
  }), use.names = TRUE, fill = TRUE)

  list(
    best = list(model = best_model, vars = best_vars),
    candidates = cand,
    X_names = X_names,
    ranking = ranking,
    n_tested = total,
    n_admissible = length(cand)
  )
}

# -----------------------------------------------------------------------------
# Largest-remainder allocation of n_draws across models with weights w.
# Guarantees the allocation sums exactly to n_draws.
# -----------------------------------------------------------------------------
.allocate_draws <- function(w, n_draws) {
  raw <- w * n_draws
  base <- floor(raw)
  rem <- n_draws - sum(base)

  if (rem > 0) {
    o <- order(raw - base, decreasing = TRUE)
    base[o[seq_len(rem)]] <- base[o[seq_len(rem)]] + 1
  }

  as.integer(base)
}

# -----------------------------------------------------------------------------
# Akaike-weight model averaging of the satellite bridge.
#
# Returns a list with the SAME fields as make_bayesian_satellite() so it can be
# passed unchanged to inject_girf_into_Z_bayes(),
# compute_factor_moment_draws_bayes() and compute_sobol_decomposition().
# Additional BMA bookkeeping is attached for the diagnostic exports.
# -----------------------------------------------------------------------------
make_bma_satellite <- function(candidates,
                               X_names,
                               model_df,
                               variables,
                               weight = c("aic", "bic", "equal"),
                               rule = c("cum_weight", "top_n", "delta"),
                               cum_threshold = 0.95,
                               top_n = 10L,
                               delta_max = 2,
                               n_draws = 10000L,
                               seed = 12345) {
  weight <- match.arg(weight)
  rule <- match.arg(rule)
  model_df <- as.data.frame(model_df)

  if (length(candidates) == 0L) {
    stop("make_bma_satellite(): no candidate specifications supplied.")
  }

  # Re-fit each admissible candidate as a proper lm so AIC/BIC/coef are exact.
  fits <- lapply(candidates, function(cc) {
    vars <- X_names[cc$cols]
    f <- as.formula(paste("Y ~", paste(vars, collapse = " + ")))
    lm(f, data = model_df)
  })

  ic <- if (weight == "bic") {
    vapply(fits, BIC, numeric(1))
  } else {
    vapply(fits, AIC, numeric(1))
  }

  ord <- order(ic)
  fits <- fits[ord]
  ic <- ic[ord]

  d_all <- ic - min(ic)
  w_all <- if (weight == "equal") {
    rep(1 / length(fits), length(fits))
  } else {
    ww <- exp(-0.5 * d_all)
    ww / sum(ww)
  }

  # Select the set of "top models" to average over.
  keep <- switch(
    rule,
    cum_weight = {
      cw <- cumsum(w_all)
      seq_len(max(1L, which(cw >= cum_threshold)[1]))
    },
    top_n = seq_len(min(as.integer(top_n), length(fits))),
    delta = which(d_all <= delta_max)
  )

  if (length(keep) == 0L) keep <- 1L

  fits <- fits[keep]
  ic <- ic[keep]
  d <- ic - min(ic)

  # Re-normalise weights over the retained set.
  w <- if (weight == "equal") {
    rep(1 / length(fits), length(fits))
  } else {
    ww <- exp(-0.5 * d)
    ww / sum(ww)
  }

  K <- length(fits)
  n_k <- .allocate_draws(w, n_draws)

  # ---------------------------------------------------------------------------
  # Per-model conjugate Jeffreys posteriors, then stack into the mixture.
  # Models allocated zero draws still contribute to the coefficient/inclusion
  # summaries (their weight is tiny) but not to the draw mixture.
  # ---------------------------------------------------------------------------
  beta0_chunks <- vector("list", K)
  sigma2_chunks <- vector("list", K)
  term_draw_chunks <- vector("list", K)  # named-column matrices, one per model
  term_lists <- vector("list", K)

  for (k in seq_len(K)) {
    bsat_k <- make_bayesian_satellite(
      sat_model = fits[[k]],
      model_df = model_df,
      variables = variables,
      n_draws = max(1L, n_k[k]),
      seed = seed + k
    )
    term_lists[[k]] <- bsat_k$terms

    if (n_k[k] >= 1L) {
      beta0_chunks[[k]] <- as.numeric(bsat_k$beta0_draws)
      sigma2_chunks[[k]] <- as.numeric(bsat_k$sigma2_draws)
      term_draw_chunks[[k]] <- bsat_k$beta_term_draws
    }
  }

  # Union of satellite terms across the retained models.
  union_terms <- unique(rbindlist(term_lists, use.names = TRUE, fill = TRUE),
                        by = "term")
  union_terms[, col := match(base, variables)]
  setorder(union_terms, lag, base)
  union_term_names <- union_terms$term

  # Assemble the mixture draw matrices (rows stacked model by model).
  beta0_draws <- unlist(beta0_chunks, use.names = FALSE)
  sigma2_draws <- unlist(sigma2_chunks, use.names = FALSE)
  total_draws <- length(beta0_draws)

  beta_term_draws <- matrix(
    0.0,
    nrow = total_draws,
    ncol = length(union_term_names),
    dimnames = list(NULL, union_term_names)
  )

  row0 <- 0L
  for (k in seq_len(K)) {
    if (is.null(term_draw_chunks[[k]])) next
    nk <- nrow(term_draw_chunks[[k]])
    cols_k <- colnames(term_draw_chunks[[k]])
    beta_term_draws[(row0 + 1L):(row0 + nk), cols_k] <- term_draw_chunks[[k]]
    row0 <- row0 + nk
  }

  # ---------------------------------------------------------------------------
  # Model-averaged point estimates and inclusion probabilities (over the full
  # retained weight vector w, including any zero-draw models).
  # ---------------------------------------------------------------------------
  avg_beta0 <- 0.0
  avg_terms <- setNames(numeric(length(union_term_names)), union_term_names)
  incl_prob <- setNames(numeric(length(union_term_names)), union_term_names)
  avg_sigma2 <- 0.0

  for (k in seq_len(K)) {
    cf <- coef(fits[[k]])
    avg_beta0 <- avg_beta0 + w[k] * unname(cf["(Intercept)"])
    avg_sigma2 <- avg_sigma2 + w[k] * (summary(fits[[k]])$sigma)^2
    tk <- term_lists[[k]]$term
    for (t in tk) {
      avg_terms[t] <- avg_terms[t] + w[k] * unname(cf[t])
      incl_prob[t] <- incl_prob[t] + w[k]
    }
  }

  beta_ols <- c(`(Intercept)` = avg_beta0, avg_terms)

  # Diagnostics: per-model table.
  models_dt <- rbindlist(lapply(seq_len(K), function(k) {
    cf <- coef(fits[[k]])
    data.table(
      rank = k,
      model = paste(setdiff(names(cf), "(Intercept)"), collapse = " + "),
      n_vars = length(cf) - 1L,
      ic = ic[k],
      delta_ic = d[k],
      weight = w[k],
      n_draws_alloc = n_k[k]
    )
  }), use.names = TRUE, fill = TRUE)

  inclusion_dt <- data.table(
    term = union_term_names,
    base = union_terms$base,
    lag = union_terms$lag,
    inclusion_prob = as.numeric(incl_prob[union_term_names]),
    avg_coefficient = as.numeric(avg_terms[union_term_names])
  )
  setorder(inclusion_dt, -inclusion_prob)

  list(
    # --- make_bayesian_satellite-compatible fields (feed the PD pipeline) ----
    terms = union_terms[, .(term, base, lag, col)],
    beta0_ols = avg_beta0,
    beta0_draws = beta0_draws,
    beta_term_draws = beta_term_draws,
    beta_ols = beta_ols,
    sigma2_draws = sigma2_draws,
    sigma2_ols = avg_sigma2,
    n_draws = total_draws,
    prior = paste0(
      "BMA mixture over top satellite specifications; within-model prior: ",
      "Jeffreys p(beta, sigma^2) proportional to 1/sigma^2; ",
      "weights: ", weight, "; selection rule: ", rule
    ),
    full = NULL,
    # --- BMA bookkeeping for diagnostic exports ------------------------------
    bma = list(
      weight_scheme = weight,
      rule = rule,
      cum_threshold = cum_threshold,
      top_n = top_n,
      delta_max = delta_max,
      n_models = K,
      models = models_dt,
      inclusion = inclusion_dt
    )
  )
}

# -----------------------------------------------------------------------------
# Publication-ready model-averaged regression table (CSV + LaTeX/booktabs).
# Built directly from the mixture posterior draws, so the reported coefficient,
# standard deviation and credible interval integrate BOTH parameter uncertainty
# (within-model Jeffreys posterior) AND model-selection uncertainty (mixture
# across specifications, including the zeros of excluded terms).
# -----------------------------------------------------------------------------
write_bma_regression_table <- function(bsat, out_dir, app_name, info_set_name,
                                       info_set_label = info_set_name,
                                       var_labels = NULL,
                                       n_obs = NA_integer_,
                                       app_label = app_name) {
  if (is.null(bsat$bma)) return(invisible(FALSE))

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  stem <- paste0(app_name, "_", info_set_name, "_")

  draws <- cbind(`(Intercept)` = bsat$beta0_draws, bsat$beta_term_draws)

  incl_map <- setNames(
    bsat$bma$inclusion$inclusion_prob,
    bsat$bma$inclusion$term
  )

  pretty_term <- function(term) {
    if (term == "(Intercept)") return("Intercept")
    g <- regmatches(term, regexec("^(.*)_lag([0-9]+)$", term))[[1]]
    if (length(g) == 3L) {
      base <- g[2]
      lab <- if (!is.null(var_labels) && base %in% names(var_labels)) {
        var_labels[[base]]
      } else {
        base
      }
      sprintf("%s (lag %s)", lab, g[3])
    } else {
      term
    }
  }

  rows <- lapply(colnames(draws), function(cn) {
    x <- draws[, cn]
    q <- stats::quantile(x, c(0.05, 0.16, 0.50, 0.84, 0.95), na.rm = TRUE)
    p_pos <- mean(x > 0, na.rm = TRUE)
    pdir <- max(p_pos, 1 - p_pos)
    stars <- if (pdir >= 0.99) {
      "***"
    } else if (pdir >= 0.975) {
      "**"
    } else if (pdir >= 0.95) {
      "*"
    } else {
      ""
    }
    incl <- if (cn == "(Intercept)") {
      1
    } else if (cn %in% names(incl_map)) {
      as.numeric(incl_map[[cn]])
    } else {
      NA_real_
    }
    data.table(
      term = cn,
      variable = pretty_term(cn),
      post_mean = mean(x, na.rm = TRUE),
      post_sd = stats::sd(x, na.rm = TRUE),
      q05 = q[[1]], q16 = q[[2]], median = q[[3]], q84 = q[[4]], q95 = q[[5]],
      prob_direction = pdir,
      inclusion_prob = incl,
      signif = stars
    )
  })

  tab <- rbindlist(rows, use.names = TRUE, fill = TRUE)
  intercept_row <- tab[term == "(Intercept)"]
  rest <- tab[term != "(Intercept)"]
  setorder(rest, -inclusion_prob, -prob_direction)
  tab <- rbind(intercept_row, rest)

  fwrite(tab, file.path(out_dir, paste0(stem, "averaged_regression.csv")))

  # ---- LaTeX (booktabs) -----------------------------------------------------
  fnum <- function(x, d = 4) {
    ifelse(is.na(x), "", formatC(x, format = "f", digits = d))
  }
  esc <- function(s) gsub("_", "\\\\_", s)

  con <- file(file.path(out_dir, paste0(stem, "averaged_regression.tex")),
              open = "wt")
  on.exit(close(con), add = TRUE)

  cat("\\begin{table}[!htbp]\n\\centering\n", file = con)
  cat("\\caption{", esc(app_label),
      " --- model-averaged satellite regression (",
      esc(info_set_label), " information set)}\n",
      sep = "", file = con)
  cat("\\label{tab:", app_name, "_", info_set_name, "_bma_satellite}\n",
      sep = "", file = con)
  cat("\\begin{tabular}{lrrrr}\n\\toprule\n", file = con)
  cat("Variable & Coef. & Post. SD & 90\\% CI & Incl. prob. \\\\\n", file = con)
  cat("\\midrule\n", file = con)

  for (i in seq_len(nrow(tab))) {
    ci <- sprintf("$[%s,\\ %s]$", fnum(tab$q05[i]), fnum(tab$q95[i]))
    coef_str <- paste0(fnum(tab$post_mean[i]), "$^{", tab$signif[i], "}$")
    coef_str <- gsub("\\$\\^\\{\\}\\$", "", coef_str)  # drop empty superscript
    incl_str <- if (is.na(tab$inclusion_prob[i])) {
      ""
    } else {
      formatC(tab$inclusion_prob[i], format = "f", digits = 2)
    }
    cat(esc(tab$variable[i]), " & ", coef_str, " & ", fnum(tab$post_sd[i]),
        " & ", ci, " & ", incl_str, " \\\\\n", sep = "", file = con)
  }

  cat("\\midrule\n", file = con)
  cat("Weighting & \\multicolumn{4}{r}{", bsat$bma$weight_scheme,
      " (", bsat$bma$rule, ")} \\\\\n", sep = "", file = con)
  cat("Models averaged & \\multicolumn{4}{r}{", bsat$bma$n_models,
      "} \\\\\n", sep = "", file = con)
  cat("Posterior draws & \\multicolumn{4}{r}{", bsat$n_draws,
      "} \\\\\n", sep = "", file = con)
  if (!is.na(n_obs)) {
    cat("Observations & \\multicolumn{4}{r}{", n_obs, "} \\\\\n",
        sep = "", file = con)
  }
  cat("\\bottomrule\n\\end{tabular}\n", file = con)
  cat("\\begin{tablenotes}\\small\n", file = con)
  cat("\\item Notes: Posterior mean, standard deviation and 90\\% credible ",
      "interval of each satellite coefficient under Akaike-weight model ",
      "averaging. Estimates integrate parameter uncertainty (Jeffreys ",
      "posterior within each model) and model-selection uncertainty (mixture ",
      "across the top specifications). Inclusion probability is the cumulative ",
      "Akaike weight of the retained models containing the term. Stars ",
      "$^{*}/^{**}/^{***}$ denote a posterior probability of the coefficient ",
      "sign exceeding $0.95/0.975/0.99$.\n", sep = "", file = con)
  cat("\\end{tablenotes}\n\\end{table}\n", file = con)

  invisible(tab)
}

# -----------------------------------------------------------------------------
# Persist BMA diagnostics to a clear, application-specific output folder.
# -----------------------------------------------------------------------------
write_bma_diagnostics <- function(bsat, out_dir, app_name, info_set_name,
                                  best_model = NULL) {
  if (is.null(bsat$bma)) {
    return(invisible(FALSE))
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  stem <- paste0(app_name, "_", info_set_name, "_")

  fwrite(
    bsat$bma$models,
    file.path(out_dir, paste0(stem, "top_models.csv"))
  )

  fwrite(
    bsat$bma$inclusion,
    file.path(out_dir, paste0(stem, "averaged_coefficients.csv"))
  )

  summ <- data.table(
    application = app_name,
    information_set = info_set_name,
    weight_scheme = bsat$bma$weight_scheme,
    selection_rule = bsat$bma$rule,
    cum_threshold = bsat$bma$cum_threshold,
    top_n = bsat$bma$top_n,
    delta_max = bsat$bma$delta_max,
    n_models_averaged = bsat$bma$n_models,
    n_posterior_draws = bsat$n_draws,
    averaged_intercept = bsat$beta0_ols,
    averaged_sigma2 = bsat$sigma2_ols,
    single_best_model = if (!is.null(best_model)) {
      paste(setdiff(names(coef(best_model)), "(Intercept)"), collapse = " + ")
    } else {
      NA_character_
    }
  )

  fwrite(summ, file.path(out_dir, paste0(stem, "summary.csv")))

  invisible(TRUE)
}
