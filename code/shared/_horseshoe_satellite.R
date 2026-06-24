# =============================================================================
#  code/shared/_horseshoe_satellite.R — Bayesian horseshoe satellite
#
#  Single Bayesian shrinkage regression of the systematic factor Z on the FULL
#  lagged satellite design, with a horseshoe prior on the slope coefficients.
#  This replaces selection + model averaging by one coherent estimator that
#  (i) handles the collinearity of lagged macro regressors, (ii) shrinks noise
#  hard while leaving strong signals nearly unbiased (so the PD-response
#  magnitude is preserved), and (iii) delivers a full posterior that is
#  propagated, unchanged, through the GIRF -> Z -> PD pipeline.
#
#  It is used as a ROBUSTNESS check on the model-averaged main specification:
#  the function returns a `bsat`-compatible object (same fields as
#  make_bayesian_satellite / make_bma_satellite), so it is a drop-in for
#  inject_girf_into_Z_bayes(), compute_factor_moment_draws_bayes() and
#  compute_pd_girf().
#
#  Sampler: Makalic & Schmidt (2016), "A Simple Sampler for the Horseshoe
#  Estimator", exact Gibbs via inverse-gamma scale mixtures. No dependency.
#
#  Requires (base R + data.table).
# =============================================================================

# Inverse-gamma with shape a and scale b: density proportional to
# x^{-a-1} exp(-b/x). Sampled as b / Gamma(shape = a, rate = 1).
.rinvgamma <- function(n, shape, scale) {
  scale / stats::rgamma(n, shape = shape, rate = 1)
}

.hs_parse_terms <- function(term_names, variables) {
  parse_one <- function(nm) {
    g <- regmatches(nm, regexec("^(.*)_lag([0-9]+)$", nm))[[1]]
    if (length(g) != 3L) stop("Cannot parse satellite term: ", nm)
    list(base = g[2], lag = as.integer(g[3]))
  }
  dt <- data.table(term = term_names, base = NA_character_, lag = NA_integer_)
  for (i in seq_len(nrow(dt))) {
    p <- parse_one(dt$term[i])
    dt$base[i] <- p$base
    dt$lag[i] <- p$lag
  }
  dt[, col := match(base, variables)]
  if (any(is.na(dt$col))) {
    stop("Horseshoe satellite terms not found in VAR variables: ",
         paste(dt$term[is.na(dt$col)], collapse = ", "))
  }
  dt
}

# -----------------------------------------------------------------------------
# Fit the horseshoe satellite and return a bsat-compatible object.
# -----------------------------------------------------------------------------
make_horseshoe_satellite <- function(model_df, variables,
                                     n_draws = 10000L,
                                     burnin = 1000L,
                                     thin = 1L,
                                     seed = 12345) {
  stopifnot(all(c("Date", "Y") %in% names(model_df)))
  model_df <- as.data.frame(model_df)
  model_df <- model_df[stats::complete.cases(
    model_df[, c("Y", setdiff(names(model_df), c("Date", "Y")))]
  ), , drop = FALSE]

  X_names <- setdiff(names(model_df), c("Date", "Y"))
  y <- as.numeric(model_df$Y)
  X <- as.matrix(model_df[, X_names, drop = FALSE])

  # Drop zero-variance columns; standardize the rest.
  sds <- apply(X, 2, stats::sd, na.rm = TRUE)
  keep <- is.finite(sds) & sds > 1e-12
  X <- X[, keep, drop = FALSE]
  X_names <- X_names[keep]
  xbar <- colMeans(X)
  xsd <- apply(X, 2, stats::sd)
  Xs <- scale(X, center = xbar, scale = xsd)
  ybar <- mean(y)
  y_c <- y - ybar

  n <- nrow(Xs)
  p <- ncol(Xs)
  if (p < 1L) stop("Horseshoe satellite: empty design after cleaning.")

  XtX <- crossprod(Xs)
  Xty <- crossprod(Xs, y_c)

  set.seed(seed)

  # Initial values.
  beta <- rep(0, p)
  sigma2 <- stats::var(y_c)
  lambda2 <- rep(1, p)
  nu <- rep(1, p)
  tau2 <- 1
  xi <- 1

  total <- as.integer(burnin + n_draws * thin)
  keep_idx <- burnin + seq_len(n_draws) * thin

  beta_std_draws <- matrix(NA_real_, nrow = n_draws, ncol = p)
  sigma2_draws <- numeric(n_draws)
  lambda2_draws <- matrix(NA_real_, nrow = n_draws, ncol = p)
  tau2_draws <- numeric(n_draws)
  store <- 0L

  for (it in seq_len(total)) {
    # --- beta | rest  (Rue 2001 sampler) ---
    A <- XtX + diag(1 / (tau2 * lambda2), p, p)
    R <- chol(A)                         # R'R = A
    w <- forwardsolve(t(R), Xty)
    mu <- backsolve(R, w)                # mu = A^{-1} X'y
    beta <- as.numeric(mu + sqrt(sigma2) * backsolve(R, stats::rnorm(p)))

    # --- sigma2 | rest ---
    resid <- y_c - Xs %*% beta
    rss <- sum(resid^2)
    quad <- sum(beta^2 / (tau2 * lambda2))
    sigma2 <- .rinvgamma(1, shape = (n + p) / 2, scale = (rss + quad) / 2)

    # --- local scales lambda2_j | rest and auxiliaries nu_j ---
    lambda2 <- .rinvgamma(p, shape = 1,
                          scale = 1 / nu + beta^2 / (2 * tau2 * sigma2))
    nu <- .rinvgamma(p, shape = 1, scale = 1 + 1 / lambda2)

    # --- global scale tau2 | rest and auxiliary xi ---
    tau2 <- .rinvgamma(1, shape = (p + 1) / 2,
                       scale = 1 / xi + sum(beta^2 / lambda2) / (2 * sigma2))
    xi <- .rinvgamma(1, shape = 1, scale = 1 + 1 / tau2)

    if (it %in% keep_idx) {
      store <- store + 1L
      beta_std_draws[store, ] <- beta
      sigma2_draws[store] <- sigma2
      lambda2_draws[store, ] <- lambda2
      tau2_draws[store] <- tau2
    }
  }

  # Back-transform standardized slopes to the raw scale used by psi_Z.
  beta_term_draws <- sweep(beta_std_draws, 2, xsd, `/`)
  colnames(beta_term_draws) <- X_names
  beta0_draws <- as.numeric(ybar - beta_term_draws %*% xbar)

  terms_df <- .hs_parse_terms(X_names, variables)

  # Posterior summaries and an (approximate) horseshoe signal weight
  # 1 - kappa_j, with kappa_j = 1 / (1 + n * tau2 * lambda2_j).
  kappa <- 1 / (1 + n * sweep(lambda2_draws, 1, tau2_draws, `*`))
  signal_weight <- 1 - colMeans(kappa)

  q <- function(x, pr) stats::quantile(x, pr, na.rm = TRUE)
  diag_dt <- data.table(
    term = X_names,
    base = terms_df$base,
    lag = terms_df$lag,
    post_mean = colMeans(beta_term_draws),
    post_sd = apply(beta_term_draws, 2, stats::sd),
    q05 = apply(beta_term_draws, 2, q, pr = 0.05),
    median = apply(beta_term_draws, 2, q, pr = 0.50),
    q95 = apply(beta_term_draws, 2, q, pr = 0.95),
    prob_direction = apply(beta_term_draws, 2,
                           function(b) max(mean(b > 0), mean(b < 0))),
    signal_weight = signal_weight
  )
  setorder(diag_dt, -signal_weight)

  beta_ols <- c(`(Intercept)` = mean(beta0_draws), colMeans(beta_term_draws))

  list(
    terms = terms_df[, .(term, base, lag, col)],
    beta0_ols = mean(beta0_draws),
    beta0_draws = beta0_draws,
    beta_term_draws = beta_term_draws,
    beta_ols = beta_ols,
    sigma2_draws = sigma2_draws,
    sigma2_ols = mean(sigma2_draws),
    n_draws = n_draws,
    prior = "Horseshoe prior on slopes (Makalic-Schmidt 2016 Gibbs); flat intercept; Jeffreys sigma^2",
    full = NULL,
    hs = list(
      diagnostics = diag_dt,
      tau2_post_mean = mean(tau2_draws),
      n_obs = n,
      n_terms = p,
      burnin = burnin,
      thin = thin
    )
  )
}

# -----------------------------------------------------------------------------
# Publication-ready horseshoe regression table (CSV + LaTeX/booktabs).
# -----------------------------------------------------------------------------
write_horseshoe_table <- function(bsat, out_dir, app_name, info_set_name,
                                  info_set_label = info_set_name,
                                  var_labels = NULL,
                                  n_obs = NA_integer_,
                                  app_label = app_name) {
  if (is.null(bsat$hs)) return(invisible(FALSE))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  stem <- paste0(app_name, "_", info_set_name, "_")

  tab <- copy(bsat$hs$diagnostics)

  pretty_term <- function(base, lag) {
    lab <- if (!is.null(var_labels) && base %in% names(var_labels)) {
      var_labels[[base]]
    } else {
      base
    }
    sprintf("%s (lag %d)", lab, lag)
  }
  tab[, variable := mapply(pretty_term, base, lag)]
  tab[, signif := fifelse(prob_direction >= 0.99, "***",
                  fifelse(prob_direction >= 0.975, "**",
                  fifelse(prob_direction >= 0.95, "*", "")))]

  fwrite(tab, file.path(out_dir, paste0(stem, "horseshoe_regression.csv")))

  fnum <- function(x, d = 4) ifelse(is.na(x), "", formatC(x, format = "f", digits = d))
  esc <- function(s) gsub("_", "\\\\_", s)

  con <- file(file.path(out_dir, paste0(stem, "horseshoe_regression.tex")),
              open = "wt")
  on.exit(close(con), add = TRUE)

  cat("\\begin{table}[!htbp]\n\\centering\n", file = con)
  cat("\\caption{", esc(app_label),
      " --- horseshoe satellite regression (",
      esc(info_set_label), " information set)}\n", sep = "", file = con)
  cat("\\label{tab:", app_name, "_", info_set_name, "_horseshoe}\n",
      sep = "", file = con)
  cat("\\begin{tabular}{lrrrr}\n\\toprule\n", file = con)
  cat("Variable & Coef. & Post. SD & 90\\% CI & Signal wt. \\\\\n", file = con)
  cat("\\midrule\n", file = con)
  cat("Intercept & ", fnum(bsat$beta0_ols), " & & & \\\\\n", sep = "", file = con)
  for (i in seq_len(nrow(tab))) {
    ci <- sprintf("$[%s,\\ %s]$", fnum(tab$q05[i]), fnum(tab$q95[i]))
    coef_str <- if (nzchar(tab$signif[i])) {
      paste0(fnum(tab$post_mean[i]), "$^{", tab$signif[i], "}$")
    } else {
      fnum(tab$post_mean[i])
    }
    cat(esc(tab$variable[i]), " & ", coef_str, " & ", fnum(tab$post_sd[i]),
        " & ", ci, " & ", fnum(tab$signal_weight[i], 2), " \\\\\n",
        sep = "", file = con)
  }
  cat("\\midrule\n", file = con)
  cat("Prior & \\multicolumn{4}{r}{horseshoe (shrinkage)} \\\\\n", file = con)
  cat("Posterior draws & \\multicolumn{4}{r}{", bsat$n_draws, "} \\\\\n",
      sep = "", file = con)
  if (!is.na(n_obs)) {
    cat("Observations & \\multicolumn{4}{r}{", n_obs, "} \\\\\n",
        sep = "", file = con)
  }
  cat("\\bottomrule\n\\end{tabular}\n", file = con)
  cat("\\begin{tablenotes}\\small\n", file = con)
  cat("\\item Notes: Posterior mean, standard deviation and 90\\% credible ",
      "interval of each satellite slope under a horseshoe shrinkage prior, on ",
      "the full lagged design. The signal weight $1-\\kappa_j$ (with ",
      "$\\kappa_j$ the local shrinkage factor) is an inclusion-like measure ",
      "in $[0,1]$, near one for retained signals and near zero for shrunk ",
      "noise. Stars $^{*}/^{**}/^{***}$ denote a posterior probability of the ",
      "coefficient sign exceeding $0.95/0.975/0.99$.\n", sep = "", file = con)
  cat("\\end{tablenotes}\n\\end{table}\n", file = con)

  invisible(tab)
}

# -----------------------------------------------------------------------------
# Engine driver: fit the horseshoe satellite, write its table, and propagate it
# through the GIRF -> Z -> PD pipeline (both shocks). Pure side-effect.
# -----------------------------------------------------------------------------
run_horseshoe_robustness <- function(out_dir, sat_df, variables, n_draws, burnin,
                                     seed, kern, impulse_ix, p_lags, H,
                                     shock_scale_2001, p, rho,
                                     app_name, app_label, info_set_name,
                                     info_set_label, var_labels, n_obs) {
  tryCatch({
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    bsat_hs <- make_horseshoe_satellite(
      model_df = sat_df, variables = variables,
      n_draws = n_draws, burnin = burnin, seed = seed
    )
    saveRDS(bsat_hs, file.path(out_dir, "horseshoe_satellite.rds"))

    write_horseshoe_table(
      bsat = bsat_hs, out_dir = out_dir, app_name = app_name,
      info_set_name = info_set_name, info_set_label = info_set_label,
      var_labels = var_labels, n_obs = n_obs, app_label = app_label
    )

    mom_hs <- compute_factor_moment_draws_bayes(kern, bsat_hs, impulse_ix,
                                                p_lags, H)

    tag <- function(dt, lab) {
      dt <- as.data.table(copy(dt))
      dt[, `:=`(information_set = info_set_name,
                info_set_label = info_set_label, shock = lab)]
      dt[]
    }

    psiZ_hs_1sd <- inject_girf_into_Z_bayes(kern$Psi_draws, kern$Sigma_kept,
                                            bsat_hs, impulse_ix, shock_scale = 1)
    pd_hs_1sd <- compute_pd_girf(psiZ_hs_1sd$draws, mom_hs$mu_draws,
                                 mom_hs$s2_draws, mom_hs$s2_delta_draws,
                                 p = p, rho = rho)

    psiZ_hs_2001 <- inject_girf_into_Z_bayes(kern$Psi_draws, kern$Sigma_kept,
                                             bsat_hs, impulse_ix,
                                             shock_scale = shock_scale_2001)
    pd_hs_2001 <- compute_pd_girf(psiZ_hs_2001$draws, mom_hs$mu_draws,
                                  mom_hs$s2_draws, mom_hs$s2_delta_draws,
                                  p = p, rho = rho)

    fwrite(tag(pd_hs_1sd$bands, "one_sd"),
           file.path(out_dir, "pd_girf_bands_1sd.csv"))
    fwrite(tag(pd_hs_2001$bands, "2001Q3"),
           file.path(out_dir, "pd_girf_bands_2001Q3.csv"))
    fwrite(tag(psiZ_hs_1sd$bands, "one_sd"),
           file.path(out_dir, "psiZ_girf_bands_1sd.csv"))

    save_png_plot(
      plot_single_bands(pd_hs_1sd$bands, ylab_txt = expression(Delta ~ PD)),
      out_dir, "pd_girf_1sd", width = 6.5, height = 4.2
    )
    save_png_plot(
      plot_single_bands(pd_hs_2001$bands, ylab_txt = expression(Delta ~ PD)),
      out_dir, "pd_girf_2001Q3", width = 6.5, height = 4.2
    )
  }, error = function(e) {
    cat("    Horseshoe robustness skipped for ", info_set_name, ": ",
        conditionMessage(e), "\n", sep = "")
  })
  invisible(TRUE)
}
