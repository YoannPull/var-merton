# =============================================================================
#  code/_bayesian_satellite_utils.R — Bayesian linear satellite helpers
#
#  The main pipeline uses the OLS-selected satellite specification only as a
#  design-selection device. Conditional on that fixed design, inference is based
#  on a conjugate Bayesian linear regression with the standard non-informative
#  prior p(beta, sigma^2) ∝ 1/sigma^2. This gives posterior beta draws centred
#  on the OLS coefficients while propagating beta and residual-variance
#  uncertainty to Z and PD responses.
# =============================================================================

rinvchisq <- function(n, df, scale) {
  if (!is.finite(df) || df <= 0) stop("df must be positive in rinvchisq().")
  if (!is.finite(scale) || scale <= 0) stop("scale must be positive in rinvchisq().")
  df * scale / stats::rchisq(n, df = df)
}

parse_lagged_satellite_terms <- function(term_names, allowed_bases, gpr_name, vars_in_var) {
  parse_one <- function(name) {
    m <- regexec("^(.*)_lag([0-9]+)$", name)
    g <- regmatches(name, m)[[1]]
    if (length(g) != 3L) return(list(base = NA_character_, lag = NA_integer_))
    list(base = g[2], lag = as.integer(g[3]))
  }

  terms_df <- data.frame(
    term = term_names,
    base = NA_character_,
    lag = NA_integer_,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(terms_df))) {
    p <- parse_one(terms_df$term[i])
    terms_df$base[i] <- p$base
    terms_df$lag[i] <- p$lag
  }

  keep <- terms_df$base %in% vars_in_var & terms_df$base != gpr_name & terms_df$base %in% allowed_bases
  if (any(!keep)) {
    warning("Ignored satellite terms outside the VAR information set: ", paste(terms_df$term[!keep], collapse = ", "))
    terms_df <- terms_df[keep, , drop = FALSE]
  }
  if (nrow(terms_df) == 0L) stop("No valid lagged satellite term remains after VAR-alignment checks.")

  terms_df$col <- match(terms_df$base, vars_in_var)
  terms_df
}

fit_bayesian_lm_jeffreys <- function(model_df, formula_obj, n_draws, seed = 123) {
  stopifnot(is.data.frame(model_df), inherits(formula_obj, "formula"))
  if (!all(c("Date", "Y") %in% names(model_df))) stop("model_df must contain Date and Y columns.")
  if (!is.finite(n_draws) || n_draws < 1) stop("n_draws must be a positive integer.")

  mf <- stats::model.frame(formula_obj, data = model_df, na.action = stats::na.omit)
  y <- stats::model.response(mf)
  X <- stats::model.matrix(formula_obj, data = mf)

  if (anyNA(y) || anyNA(X)) stop("Bayesian satellite design contains NA after model.frame construction.")
  if (nrow(X) != length(y)) stop("Incompatible X/y dimensions in Bayesian satellite.")
  if (qr(X)$rank < ncol(X)) stop("Bayesian satellite design is rank-deficient.")

  n <- nrow(X)
  p <- ncol(X)
  df <- n - p
  if (df <= 0L) stop("Not enough observations for the selected satellite design: n <= p.")

  XtX_inv <- solve(crossprod(X))
  beta_hat <- as.numeric(XtX_inv %*% crossprod(X, y))
  names(beta_hat) <- colnames(X)

  resid <- as.numeric(y - X %*% beta_hat)
  rss <- sum(resid^2)
  sigma2_hat_unbiased <- rss / df

  set.seed(seed)
  sigma2_draws <- rinvchisq(n_draws, df = df, scale = sigma2_hat_unbiased)

  R <- chol(XtX_inv)
  z <- matrix(stats::rnorm(n_draws * p), nrow = n_draws, ncol = p)
  beta_draws <- sweep(z %*% R, 1, sqrt(sigma2_draws), `*`)
  beta_draws <- sweep(beta_draws, 2, beta_hat, `+`)
  colnames(beta_draws) <- colnames(X)

  fitted_mean <- as.numeric(X %*% beta_hat)

  list(
    prior = "Jeffreys non-informative prior: p(beta, sigma^2) proportional to 1/sigma^2",
    formula = formula_obj,
    n_obs = n,
    n_coef = p,
    df_residual = df,
    y_name = all.vars(formula_obj)[1],
    coef_names = colnames(X),
    beta_ols = beta_hat,
    sigma2_ols_unbiased = sigma2_hat_unbiased,
    rss = rss,
    XtX_inv = XtX_inv,
    beta_draws = beta_draws,
    sigma2_draws = sigma2_draws,
    fitted_mean = fitted_mean,
    residuals_ols = resid,
    dates = as.Date(model_df$Date[as.integer(rownames(mf))]),
    model_df = model_df
  )
}

summarise_draw_matrix <- function(draws, ols = NULL, probs = c(0.05, 0.16, 0.50, 0.84, 0.95)) {
  out <- data.frame(
    term = colnames(draws),
    mean = colMeans(draws, na.rm = TRUE),
    sd = apply(draws, 2, stats::sd, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  qs <- t(apply(draws, 2, stats::quantile, probs = probs, na.rm = TRUE))
  colnames(qs) <- c("q05", "q16", "median", "q84", "q95")
  out <- cbind(out, as.data.frame(qs, row.names = NULL))
  if (!is.null(ols)) {
    out$ols <- as.numeric(ols[out$term])
    out$posterior_mean_minus_ols <- out$mean - out$ols
  }
  out
}

summarise_vector_draws <- function(x, name = "sigma2", probs = c(0.05, 0.16, 0.50, 0.84, 0.95)) {
  qs <- stats::quantile(x, probs = probs, na.rm = TRUE)
  data.frame(
    quantity = name,
    mean = mean(x, na.rm = TRUE),
    sd = stats::sd(x, na.rm = TRUE),
    q05 = qs[[1]], q16 = qs[[2]], median = qs[[3]], q84 = qs[[4]], q95 = qs[[5]],
    stringsAsFactors = FALSE
  )
}

load_bayesian_satellite <- function(path, allowed_bases, gpr_name, vars_in_var) {
  if (!file.exists(path)) stop("Missing Bayesian satellite posterior: ", path)
  obj <- readRDS(path)
  required <- c("beta_draws", "sigma2_draws", "beta_ols", "formula", "model_df")
  missing <- setdiff(required, names(obj))
  if (length(missing) > 0L) stop("Invalid Bayesian satellite object; missing: ", paste(missing, collapse = ", "))

  beta_draws <- obj$beta_draws
  intercept_col <- match("(Intercept)", colnames(beta_draws))
  if (is.na(intercept_col)) stop("Bayesian satellite draws must contain an intercept column.")

  term_names <- setdiff(colnames(beta_draws), "(Intercept)")
  terms_df <- parse_lagged_satellite_terms(term_names, allowed_bases, gpr_name, vars_in_var)

  keep_terms <- c("(Intercept)", terms_df$term)
  obj$beta_draws <- beta_draws[, keep_terms, drop = FALSE]
  obj$beta0_draws <- obj$beta_draws[, "(Intercept)"]
  obj$beta_term_draws <- obj$beta_draws[, terms_df$term, drop = FALSE]
  obj$terms <- terms_df
  obj
}

pair_satellite_draw <- function(m, n_sat) {
  ((m - 1L) %% n_sat) + 1L
}

# =============================================================================
#  Sobol variance decomposition of the PD response (pick-freeze design)
# -----------------------------------------------------------------------------
#  Exact variance-based decomposition of Var(Delta PD_h) into the contributions
#  of the two independent posterior blocks:
#     V  = BVAR draw (Psi, Sigma, B)
#     S  = Bayesian satellite draw (beta0, beta, sigma_eta2)
#  Using the law of total variance / functional ANOVA:
#     Var = V_V + V_S + V_VS
#  with first-order effects V_V = Var_V[E_S g], V_S = Var_S[E_V g] estimated by
#  the Saltelli (2010) pick-freeze estimator, the interaction V_VS by residual,
#  and total-effect indices by the Jansen (1999) estimator. No "fix-at-mean" and
#  no Gaussian-band approximation: this is the assumption-free decomposition.
#
#  Requires (available at call time): build_selection_from_terms(), build_G_hq(),
#  forecast_baseline_path().
# =============================================================================

# Per-draw conditional-variance kernels of Z (without folding in beta), so that
# different satellite draws can be applied to the same BVAR draw cheaply.
compute_s2_kernels_one_draw <- function(Psi_m, Sigma_m, terms_df, variables,
                                        impulse_idx = 1, jitter = 1e-10) {
  H <- dim(Psi_m)[1] - 1
  Sigma_m <- 0.5 * (Sigma_m + t(Sigma_m))
  k <- ncol(Sigma_m)
  sjj <- Sigma_m[impulse_idx, impulse_idx]
  if (!is.finite(sjj) || sjj <= 0) {
    Sigma_m <- Sigma_m + diag(jitter, k)
    sjj <- Sigma_m[impulse_idx, impulse_idx]
  }
  Sigma_cond <- Sigma_m - tcrossprod(Sigma_m[, impulse_idx], Sigma_m[, impulse_idx]) / sjj
  Sigma_cond <- 0.5 * (Sigma_cond + t(Sigma_cond))

  sel <- build_selection_from_terms(terms_df, variables)
  S_list <- sel$S_list
  Lmax <- sel$Lmax
  m <- nrow(S_list[[1]])

  Ss <- array(0.0, dim = c(m, m, H + 1))
  Ssc <- array(0.0, dim = c(m, m, H + 1))
  for (hh in 0:H) {
    A <- matrix(0.0, m, m)
    Ac <- matrix(0.0, m, m)
    for (q in 0:hh) {
      G <- build_G_hq(hh, q, S_list, Lmax, Psi_m)
      A <- A + G %*% Sigma_m %*% t(G)
      if (q == 0) Ac <- Ac + G %*% Sigma_cond %*% t(G) else Ac <- Ac + G %*% Sigma_m %*% t(G)
    }
    Ss[, , hh + 1] <- A
    Ssc[, , hh + 1] <- Ac
  }
  list(Ss = Ss, Ssc = Ssc)
}

compute_sobol_decomposition <- function(Psi_draws, Sigma_kept, B_kept, variables, DT_var,
                                        terms_df, beta0_draws, beta_term_draws, sigma2_draws,
                                        impulse_idx, p_lags, H, p_val, rho_val,
                                        shock_scale = 1, seed = 1L, jitter = 1e-10) {
  M <- dim(Psi_draws)[4]
  n_sat <- nrow(beta_term_draws)
  k <- length(variables)
  m_terms <- nrow(terms_df)
  qpi <- qnorm(p_val)

  N <- min(M %/% 2L, n_sat %/% 2L)
  if (N < 2L) stop("compute_sobol_decomposition(): not enough draws for a pick-freeze design.")

  set.seed(seed)
  Vsel <- sample.int(M, 2L * N)
  Ssel <- sample.int(n_sat, 2L * N)
  posA <- 1:N
  posB <- (N + 1L):(2L * N)

  # Baseline history (shared across draws).
  DT_var <- as.data.table(DT_var)
  Tn <- nrow(DT_var)
  Lmax <- max(terms_df$lag)
  Y_all <- as.matrix(DT_var[, ..variables])
  Y_hist <- Y_all[(Tn - p_lags + 1):Tn, , drop = FALSE]
  hist_tail <- Y_all[(Tn - Lmax):Tn, , drop = FALSE]
  cols <- terms_df$col
  lags <- terms_df$lag

  # Per-V precomputation (impulse path, baseline path, s2 kernels), indexed 1..2N.
  girfY_list <- vector("list", 2L * N)
  Yfore_list <- vector("list", 2L * N)
  Ss_list <- vector("list", 2L * N)
  Ssc_list <- vector("list", 2L * N)

  for (pos in seq_len(2L * N)) {
    v <- Vsel[pos]
    Sigma_v <- 0.5 * (Sigma_kept[, , v] + t(Sigma_kept[, , v]))
    sjj <- Sigma_v[impulse_idx, impulse_idx]
    if (!is.finite(sjj) || sjj <= 0) {
      Sigma_v <- Sigma_v + diag(jitter, k)
      sjj <- Sigma_v[impulse_idx, impulse_idx]
    }
    delta_unit <- Sigma_v[, impulse_idx] / sqrt(sjj)
    Psi_v <- Psi_draws[, , , v, drop = FALSE][, , , 1]
    girfY <- matrix(0.0, H + 1, k)
    for (hh in 0:H) girfY[hh + 1, ] <- Psi_v[hh + 1, , ] %*% delta_unit
    girfY_list[[pos]] <- girfY
    Yfore_list[[pos]] <- forecast_baseline_path(B_kept[, , v], Y_hist, p_lags, H)
    ker <- compute_s2_kernels_one_draw(Psi_v, Sigma_kept[, , v], terms_df, variables, impulse_idx, jitter)
    Ss_list[[pos]] <- ker$Ss
    Ssc_list[[pos]] <- ker$Ssc
  }

  Bmat <- beta_term_draws[, terms_df$term, drop = FALSE]

  # Evaluate Delta PD over horizons for paired positions (V-position, S-index).
  eval_g <- function(vpos, sidx) {
    out <- matrix(NA_real_, H + 1, length(vpos))
    for (j in seq_along(vpos)) {
      pos <- vpos[j]
      s <- sidx[j]
      beta_s <- as.numeric(Bmat[s, ])
      b0 <- beta0_draws[s]
      se2 <- sigma2_draws[s]
      girfY <- girfY_list[[pos]]
      Yfore <- Yfore_list[[pos]]
      Ss <- Ss_list[[pos]]
      Ssc <- Ssc_list[[pos]]

      psiZ <- numeric(H + 1)
      mu <- rep(b0, H + 1)
      for (r in seq_len(m_terms)) {
        lag_r <- lags[r]
        col_r <- cols[r]
        contrib <- if (lag_r == 0) girfY[, col_r] else c(rep(0, lag_r), girfY[1:(H + 1 - lag_r), col_r])
        psiZ <- psiZ + beta_s[r] * contrib
        for (hh in 0:H) {
          val <- if (hh >= lag_r) Yfore[hh - lag_r + 1, col_r] else hist_tail[(Lmax + 1) - (lag_r - hh), col_r]
          mu[hh + 1] <- mu[hh + 1] + beta_s[r] * val
        }
      }
      psiZ <- shock_scale * psiZ

      s2 <- numeric(H + 1)
      s2d <- numeric(H + 1)
      for (hh in 0:H) {
        s2[hh + 1] <- as.numeric(t(beta_s) %*% Ss[, , hh + 1] %*% beta_s) + se2
        s2d[hh + 1] <- as.numeric(t(beta_s) %*% Ssc[, , hh + 1] %*% beta_s) + se2
      }
      s2 <- pmax(s2, 0)
      s2d <- pmax(s2d, 0)
      s_h <- sqrt((1 - rho_val) + rho_val * s2)
      s_hd <- sqrt((1 - rho_val) + rho_val * s2d)
      pd_base <- pnorm((qpi - sqrt(rho_val) * mu) / s_h)
      pd_shock <- pnorm((qpi - sqrt(rho_val) * (mu + psiZ)) / s_hd)
      out[, j] <- pd_shock - pd_base
    }
    out
  }

  yA <- eval_g(posA, Ssel[posA])                 # (V_A, S_A)
  yB <- eval_g(posB, Ssel[posB])                 # (V_B, S_B)
  yCV <- eval_g(posB, Ssel[posA])                # A_B^{(V)} = (V_B, S_A)
  yCS <- eval_g(posA, Ssel[posB])                # A_B^{(S)} = (V_A, S_B)

  var_total <- var_bvar <- var_satellite <- var_interaction <- numeric(H + 1)
  te_bvar <- te_sat <- numeric(H + 1)
  for (h in seq_len(H + 1)) {
    a <- yA[h, ]; b <- yB[h, ]; cv <- yCV[h, ]; cs <- yCS[h, ]
    tot <- stats::var(c(a, b))
    VV <- mean(b * (cv - a))          # Saltelli 2010 first-order, factor V
    VS <- mean(b * (cs - a))          # first-order, factor S
    VTV <- mean((a - cv)^2) / 2       # Jansen total-effect, V
    VTS <- mean((a - cs)^2) / 2       # Jansen total-effect, S
    var_total[h] <- tot
    var_bvar[h] <- VV
    var_satellite[h] <- VS
    var_interaction[h] <- tot - VV - VS
    te_bvar[h] <- VTV
    te_sat[h] <- VTS
  }

  out <- data.table(
    horizon = 0:H,
    var_total = var_total,
    var_bvar = var_bvar,
    var_satellite = var_satellite,
    var_interaction = var_interaction
  )
  out[, share_bvar := fifelse(var_total > 0, 100 * var_bvar / var_total, NA_real_)]
  out[, share_satellite := fifelse(var_total > 0, 100 * var_satellite / var_total, NA_real_)]
  out[, share_interaction := fifelse(var_total > 0, 100 * var_interaction / var_total, NA_real_)]
  out[, total_effect_bvar := fifelse(var_total > 0, 100 * te_bvar / var_total, NA_real_)]
  out[, total_effect_satellite := fifelse(var_total > 0, 100 * te_sat / var_total, NA_real_)]
  out[, n_pick_freeze := N]
  out[, method := "Sobol pick-freeze (Saltelli 2010 first-order; Jansen total-effect); V=BVAR, S=satellite, independent blocks"]
  out[]
}
