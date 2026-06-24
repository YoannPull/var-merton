# =============================================================================
#  code/shared/_helpers_var_merton.R
#
#  Shared engine for the robustness scripts. This file defines, in one place,
#  every reusable building block: the Bayesian VAR sampler, the GIRF machinery
#  (Pesaran-Shin / Koop-Pesaran-Potter), the Merton-Vasicek mapping, the
#  satellite selection routines, the Z-factor reconstruction, the plotting
#  helpers and the LaTeX writers.
#
#  It is meant to be sourced AFTER code/00_setup.R, e.g.
#
#      source("code/00_setup.R")
#      source("code/shared/_helpers_var_merton.R")
#
#  No script-specific paths or parameters live here: everything that depends on
#  a given robustness exercise (output directories, shock targets, etc.) is
#  passed as an argument so the same functions serve every script.
#
#  Methodological note (identical to the main pipeline):
#  - GIRFs are computed explicitly as Pesaran-Shin / Koop-Pesaran-Potter GIRFs:
#        delta_j = Sigma[, j] / sqrt(Sigma[j, j])
#    For a GPR shock with GPR ordered first this is numerically equal to the
#    first Cholesky impulse response of a Gaussian VAR, but the shock definition
#    does not rely on Cholesky.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(mniw)
  library(parallel)
  library(ggplot2)
  library(scales)
  library(matrixStats)
})

# -----------------------------------------------------------------------------
# Parallelism control
# -----------------------------------------------------------------------------

if (is.null(getOption("robustness.core_fraction"))) {
  options(robustness.core_fraction = 0.6)
}

Sys.setenv(
  VECLIB_MAXIMUM_THREADS = "1",
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1"
)

if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  try(RhpcBLASctl::blas_set_num_threads(1), silent = TRUE)
  try(RhpcBLASctl::omp_set_num_threads(1), silent = TRUE)
}

robustness_cores <- function(fraction = getOption("robustness.core_fraction", 0.6)) {
  nc <- parallel::detectCores(logical = TRUE)
  if (!is.finite(nc) || nc < 1L) nc <- 1L
  max(1L, floor(nc * fraction))
}

# -----------------------------------------------------------------------------
# 0. Small utilities
# -----------------------------------------------------------------------------

safe_fread <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
  fread(path)
}

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "", sprintf(paste0("%.", digits, "f"), x))
}

make_quarter <- function(date) {
  paste0(
    format(date, "%Y"),
    "Q",
    ((as.integer(format(date, "%m")) - 1L) %/% 3L) + 1L
  )
}

par_lapply_safe <- function(X, FUN, mc.cores = robustness_cores()) {
  if (.Platform$OS.type == "windows" || mc.cores <= 1L) {
    lapply(X, FUN)
  } else {
    parallel::mclapply(X, FUN, mc.cores = mc.cores)
  }
}

# -----------------------------------------------------------------------------
# 1. Z-factor reconstruction
# -----------------------------------------------------------------------------

f_Z_estimation <- function(vect_TD, rho_fixed = NULL, eps = 1e-8) {
  vect_TD <- as.numeric(vect_TD)
  vect_TD <- pmin(pmax(vect_TD, eps), 1 - eps)
  
  p_ttc <- mean(vect_TD, na.rm = TRUE)
  
  compute_Z <- function(rho) {
    (qnorm(p_ttc) - qnorm(vect_TD) * sqrt(1 - rho)) / sqrt(rho)
  }
  
  if (!is.null(rho_fixed)) {
    if (!is.finite(rho_fixed) || rho_fixed <= 0 || rho_fixed >= 1) {
      stop("rho_fixed must belong to (0, 1).")
    }
    Z_estim <- compute_Z(rho_fixed)
    return(list(
      Z = Z_estim,
      rho = rho_fixed,
      p_ttc = p_ttc,
      rho_type = "fixed",
      var_Z_raw = var(Z_estim, na.rm = TRUE)
    ))
  }
  
  var_Z_rho <- function(rho) {
    if (rho <= 0 || rho >= 1) return(Inf)
    var(compute_Z(rho), na.rm = TRUE) - 1
  }
  
  result <- tryCatch(
    uniroot(var_Z_rho, lower = 1e-6, upper = 1 - 1e-6),
    error = function(e) NULL
  )
  
  if (is.null(result)) {
    return(list(
      Z = rep(NA_real_, length(vect_TD)),
      rho = NA_real_,
      p_ttc = p_ttc,
      rho_type = "estimated_failed",
      var_Z_raw = NA_real_
    ))
  }
  
  rho_opt <- result$root
  Z_estim <- compute_Z(rho_opt)
  
  list(
    Z = Z_estim,
    rho = rho_opt,
    p_ttc = p_ttc,
    rho_type = "estimated",
    var_Z_raw = var(Z_estim, na.rm = TRUE)
  )
}

# -----------------------------------------------------------------------------
# 2. Bayesian VAR and companion-matrix diagnostics
# -----------------------------------------------------------------------------
# simulate_bvar_niw(), build_companion(), max_root(), is_stable(),
# extract_c_A_list(), compute_ma_coefficients() and forecast_baseline_path()
# are defined once in code/_bvar_niw_utils.R and loaded via code/00_setup.R
# (this engine is always sourced AFTER 00_setup.R). Do not redefine them here.
# -----------------------------------------------------------------------------

estimate_bvar_kernel <- function(DT_raw, spec_name, vars_extra,
                                 p_lags, nrep, H, seed, M_target,
                                 kernel_dir, impulse_name = "log_GPRD") {
  i_var <- c(impulse_name, vars_extra)
  
  missing <- setdiff(i_var, names(DT_raw))
  if (length(missing) > 0) {
    stop("Missing variables for ", spec_name, ": ", paste(missing, collapse = ", "))
  }
  
  DT <- copy(DT_raw[, ..i_var])
  DT <- na.omit(DT)
  
  Y <- as.matrix(DT)
  colnames(Y) <- i_var
  k <- ncol(Y)
  
  res <- simulate_bvar_niw(Y, p = p_lags, nrep = nrep, seed = seed)
  
  roots_all <- vapply(
    seq_len(dim(res$B)[3]),
    function(i) max_root(res$B[, , i], p_lags, k),
    numeric(1)
  )
  
  stable <- roots_all < 1
  idx_kept <- which(stable)
  
  if (length(idx_kept) >= M_target) {
    set.seed(seed)
    idx_kept <- sample(idx_kept, M_target)
  }
  
  if (length(idx_kept) == 0) {
    stop("No stable draw for ", spec_name)
  }
  
  B_kept <- res$B[, , idx_kept, drop = FALSE]
  Sigma_kept <- res$S[, , idx_kept, drop = FALSE]
  
  psi_list <- par_lapply_safe(seq_along(idx_kept), function(jj) {
    compute_ma_coefficients(B_kept[, , jj], p = p_lags, H = H)
  })
  
  Psi_draws <- array(
    NA_real_,
    dim = c(H + 1, k, k, length(idx_kept)),
    dimnames = list(horizon = 0:H, variable = i_var, shock_col = i_var, draw = NULL)
  )
  
  for (jj in seq_along(idx_kept)) {
    Psi_draws[, , , jj] <- psi_list[[jj]]
  }
  
  kernel <- list(
    spec_name = spec_name,
    variables = i_var,
    DT = DT,
    Psi_draws = Psi_draws,
    B_kept = B_kept,
    Sigma_kept = Sigma_kept,
    p = p_lags,
    h = H,
    roots_all = roots_all,
    stable = stable,
    idx_kept = idx_kept
  )
  
  if (!is.null(kernel_dir)) {
    dir.create(kernel_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(kernel, file = file.path(kernel_dir, paste0(spec_name, "_var_kernel.rds")))
  }
  
  stable_roots <- roots_all[stable]
  
  table_row <- data.table(
    spec = spec_name,
    k = k,
    k_minus_1 = length(vars_extra),
    p = p_lags,
    variables = paste(vars_extra, collapse = ", "),
    n_draws = nrep,
    n_stable = sum(stable),
    stable_share = mean(stable),
    lambda_max_median_all = median(roots_all, na.rm = TRUE),
    lambda_max_p95_all = as.numeric(quantile(roots_all, 0.95, na.rm = TRUE)),
    lambda_max_median_stable = median(stable_roots, na.rm = TRUE),
    lambda_max_p95_stable = as.numeric(quantile(stable_roots, 0.95, na.rm = TRUE))
  )
  
  list(kernel = kernel, table_row = table_row)
}

get_or_estimate_kernel <- function(spec_name, vars_extra, DT_raw,
                                   p_lags, nrep, H, seed, M_target,
                                   kernel_dir, impulse_name = "log_GPRD") {
  kpath <- file.path(kernel_dir, paste0(spec_name, "_var_kernel.rds"))
  
  if (file.exists(kpath)) {
    message("Reusing existing kernel: ", kpath)
    return(readRDS(kpath))
  }
  
  message("Kernel not found, estimating ", spec_name, " ...")
  
  estimate_bvar_kernel(
    DT_raw = DT_raw,
    spec_name = spec_name,
    vars_extra = vars_extra,
    p_lags = p_lags,
    nrep = nrep,
    H = H,
    seed = seed,
    M_target = M_target,
    kernel_dir = kernel_dir,
    impulse_name = impulse_name
  )$kernel
}

# -----------------------------------------------------------------------------
# 3. GIRFs for the macro block
# -----------------------------------------------------------------------------

compute_var_girf <- function(Psi_draws, Sigma_kept, impulse_idx = 1, shock_scale = 1,
                             qs = c(0.05, 0.16, 0.50, 0.84, 0.95)) {
  H <- dim(Psi_draws)[1] - 1
  k <- dim(Psi_draws)[2]
  M <- dim(Psi_draws)[4]
  variables <- dimnames(Psi_draws)$variable
  
  draws <- array(
    NA_real_,
    dim = c(H + 1, k, M),
    dimnames = list(horizon = 0:H, variable = variables, draw = NULL)
  )
  
  for (m in seq_len(M)) {
    Sigma_m <- Sigma_kept[, , m]
    Sigma_m <- 0.5 * (Sigma_m + t(Sigma_m))
    
    sigma_jj <- Sigma_m[impulse_idx, impulse_idx]
    if (!is.finite(sigma_jj) || sigma_jj <= 0) next
    
    delta <- shock_scale * Sigma_m[, impulse_idx] / sqrt(sigma_jj)
    Psi_m <- Psi_draws[, , , m, drop = FALSE][, , , 1]
    
    for (hh in 0:H) {
      draws[hh + 1, , m] <- Psi_m[hh + 1, , ] %*% delta
    }
  }
  
  out <- vector("list", k)
  for (j in seq_len(k)) {
    mat <- draws[, j, , drop = FALSE]
    qh <- t(apply(mat[, 1, ], 1, quantile, probs = qs, na.rm = TRUE))
    out[[j]] <- data.table(
      horizon = 0:H,
      lower90 = qh[, 1],
      lower68 = qh[, 2],
      median = qh[, 3],
      upper68 = qh[, 4],
      upper90 = qh[, 5],
      variable = variables[j]
    )
  }
  
  list(draws = draws, bands = rbindlist(out))
}

compute_structural_e1_median <- function(kernel, dates_raw, target_quarter = "2001Q3",
                                         impulse_idx = 1) {
  variables <- kernel$variables
  p <- kernel$p
  B_bar <- apply(kernel$B_kept, c(1, 2), median)
  Sigma_bar <- apply(kernel$Sigma_kept, c(1, 2), median)
  Sigma_bar <- 0.5 * (Sigma_bar + t(Sigma_bar))
  
  Y_mat <- as.matrix(kernel$DT[, ..variables])
  
  Y_t <- Y_mat[(p + 1):nrow(Y_mat), , drop = FALSE]
  X_t <- do.call(
    cbind,
    lapply(1:p, function(l) Y_mat[(p + 1 - l):(nrow(Y_mat) - l), , drop = FALSE])
  )
  X_t <- cbind(1, X_t)
  
  U <- Y_t - X_t %*% B_bar
  
  sigma_gg <- Sigma_bar[impulse_idx, impulse_idx]
  e1 <- as.numeric(U[, impulse_idx] / sqrt(sigma_gg))
  
  dates_eff <- if (!is.null(dates_raw)) {
    as.Date(dates_raw[(p + 1):length(dates_raw)])
  } else {
    seq_len(length(e1))
  }
  
  out <- data.table(date = dates_eff, quarter = make_quarter(dates_eff), e1 = e1)
  target <- out[quarter == target_quarter]
  
  if (nrow(target) == 0) {
    warning("Target quarter ", target_quarter, " not found. Using largest 2001 |e1|.")
    target <- out[format(date, "%Y") == "2001"][which.max(abs(e1))]
  }
  
  list(e1_dt = out, shock_scale = as.numeric(target$e1[1]), target = target[1])
}

# -----------------------------------------------------------------------------
# 4. Satellite: data preparation and model selection
# -----------------------------------------------------------------------------

prepare_lagged_satellite_data <- function(z_input, DT_raw, vars_sat, nb_lags = 4) {
  z_dt <- if (is.character(z_input)) safe_fread(z_input) else as.data.table(copy(z_input))
  
  if (!"Date" %in% names(z_dt)) stop("z input must contain Date.")
  if (!"Z" %in% names(z_dt)) stop("z input must contain Z.")
  
  z_dt[, Date := as.Date(Date)]
  z_dt[, Date := as.Date(format(Date, "%Y-%m-01"))]
  
  macro <- copy(DT_raw)
  macro[, Date := as.Date(Date)]
  macro[, Date := as.Date(format(Date + 62, "%Y-%m-01"))]
  
  lag_indices <- 0:nb_lags
  
  lagged_list <- lapply(vars_sat, function(v) {
    sapply(lag_indices, function(k) data.table::shift(macro[[v]], n = k, type = "lag"))
  })
  
  X <- as.data.frame(do.call(cbind, lagged_list))
  names(X) <- paste0(rep(vars_sat, each = length(lag_indices)), "_lag", lag_indices)
  X$Date <- macro$Date
  
  df <- merge(
    z_dt[, .(Date, Y = Z)],
    as.data.table(X),
    by = "Date",
    all = FALSE
  )
  
  na.omit(df)
}

select_satellite_ols <- function(model_df, max_vars = 4, p_threshold = 0.10) {
  X_names <- setdiff(names(model_df), c("Date", "Y"))
  
  combs <- unlist(
    lapply(1:min(max_vars, length(X_names)), function(k) {
      combn(X_names, k, simplify = FALSE)
    }),
    recursive = FALSE
  )
  
  evaluate <- function(vars) {
    f <- as.formula(paste("Y ~", paste(vars, collapse = " + ")))
    fit <- try(lm(f, data = model_df), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    
    pvals <- try(summary(fit)$coefficients[-1, "Pr(>|t|)"], silent = TRUE)
    if (inherits(pvals, "try-error") || any(is.na(pvals))) return(NULL)
    if (any(pvals > p_threshold)) return(NULL)
    
    list(model = fit, vars = vars, aic = AIC(fit), bic = BIC(fit))
  }
  
  results <- lapply(combs, evaluate)
  valid <- Filter(Negate(is.null), results)
  
  if (length(valid) == 0) {
    stop("No valid satellite specification found.")
  }
  
  best <- valid[[which.min(sapply(valid, `[[`, "aic"))]]
  list(best = best, all_valid = valid, n_tested = length(combs))
}

rolling_oos_error <- function(df, form, initial, h = 1, step = 1) {
  n <- nrow(df)
  preds <- rep(NA_real_, n)
  
  for (t in seq(initial, n - h, by = step)) {
    fit <- try(lm(form, data = df[1:t, ]), silent = TRUE)
    if (inherits(fit, "try-error")) next
    preds[t + h] <- predict(fit, newdata = df[t + h, , drop = FALSE])
  }
  
  ok <- !is.na(preds)
  y <- df$Y
  
  list(
    rmse = sqrt(mean((y[ok] - preds[ok])^2)),
    mae  = mean(abs(y[ok] - preds[ok])),
    r2   = 1 - sum((y[ok] - preds[ok])^2) / sum((y[ok] - mean(y[ok]))^2),
    preds = preds
  )
}

select_satellite_oos <- function(model_df, max_vars = 6, p_threshold = NULL,
                                 initial = NULL, h = 1,
                                 preselect_by_aic_topn = 500L,
                                 verbose = TRUE) {
  model_df <- as.data.table(model_df)
  X_names <- setdiff(names(model_df), c("Date", "Y"))
  n <- nrow(model_df)
  
  if (is.null(initial)) {
    initial <- max(24L, ceiling(n / 2))
  }
  
  y <- as.numeric(model_df$Y)
  Xfull <- as.matrix(model_df[, ..X_names])
  P <- length(X_names)
  
  combo_idx <- unlist(
    lapply(1:min(max_vars, P), function(k) combn(P, k, simplify = FALSE)),
    recursive = FALSE
  )
  
  total <- length(combo_idx)
  
  if (verbose) {
    message(sprintf("OOS selection: %d candidate models over %d predictors.", total, P))
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
    stop("No admissible satellite specification found.")
  }
  
  if (is.finite(preselect_by_aic_topn) && length(cand) > preselect_by_aic_topn) {
    ord <- order(vapply(cand, `[[`, numeric(1), "aic"))
    cand <- cand[ord[seq_len(preselect_by_aic_topn)]]
  }
  
  if (verbose) {
    message(sprintf("Scoring %d models by rolling OOS RMSE ...", length(cand)))
  }
  
  oos_fast <- function(cols) {
    preds <- rep(NA_real_, n)
    
    for (t in seq(initial, n - h, by = 1)) {
      fit <- .lm.fit(cbind(1, Xfull[1:t, cols, drop = FALSE]), y[1:t])
      preds[t + h] <- sum(c(1, Xfull[t + h, cols]) * fit$coefficients)
    }
    
    ok <- !is.na(preds)
    
    list(
      rmse = sqrt(mean((y[ok] - preds[ok])^2)),
      mae  = mean(abs(y[ok] - preds[ok])),
      r2   = 1 - sum((y[ok] - preds[ok])^2) / sum((y[ok] - mean(y[ok]))^2)
    )
  }
  
  scored <- par_lapply_safe(cand, function(cc) {
    o <- oos_fast(cc$cols)
    vars <- X_names[cc$cols]
    fit <- lm(as.formula(paste("Y ~", paste(vars, collapse = " + "))), data = model_df)
    
    list(vars = vars, aic = AIC(fit), rmse = o$rmse, mae = o$mae, r2_oos = o$r2)
  })
  
  tab <- rbindlist(lapply(scored, function(s) {
    data.table(
      model = paste(s$vars, collapse = " + "),
      n_vars = length(s$vars),
      aic = s$aic,
      oos_rmse = s$rmse,
      oos_mae = s$mae,
      oos_r2 = s$r2_oos
    )
  }), use.names = TRUE, fill = TRUE)
  
  setorder(tab, oos_rmse, aic)
  
  rank_idx <- order(
    vapply(scored, `[[`, numeric(1), "rmse"),
    vapply(scored, `[[`, numeric(1), "aic")
  )
  
  best_vars <- scored[[rank_idx[1]]]$vars
  
  best_model <- lm(
    as.formula(paste("Y ~", paste(best_vars, collapse = " + "))),
    data = model_df
  )
  
  best_oos <- rolling_oos_error(model_df, formula(best_model), initial = initial, h = h)
  
  list(
    best = list(model = best_model, vars = best_vars, oos = best_oos),
    ranking = tab,
    n_tested = total,
    n_admissible = length(cand),
    initial = initial
  )
}

# -----------------------------------------------------------------------------
# 5. Propagating GIRFs through satellite and PD
# -----------------------------------------------------------------------------

parse_satellite_terms <- function(model, vars_in_var) {
  betas <- coef(model)
  beta0 <- unname(betas["(Intercept)"])
  betas <- betas[names(betas) != "(Intercept)"]
  
  parse_one <- function(term) {
    m <- regexec("^(.*)_lag([0-9]+)$", term)
    g <- regmatches(term, m)[[1]]
    if (length(g) != 3) stop("Cannot parse satellite term: ", term)
    list(base = g[2], lag = as.integer(g[3]))
  }
  
  out <- data.table(
    term = names(betas),
    base = NA_character_,
    lag = NA_integer_,
    beta = as.numeric(betas)
  )
  
  for (i in seq_len(nrow(out))) {
    z <- parse_one(out$term[i])
    out$base[i] <- z$base
    out$lag[i] <- z$lag
  }
  
  out[, col := match(base, vars_in_var)]
  
  if (any(is.na(out$col))) {
    stop(
      "Satellite terms not found in VAR variables: ",
      paste(out[is.na(col)]$term, collapse = ", ")
    )
  }
  
  list(beta0 = beta0, terms = out)
}

build_selection_from_terms <- function(terms_df, variables) {
  bases <- terms_df$base
  lags <- terms_df$lag
  
  m <- length(bases)
  k <- length(variables)
  Lmax <- max(lags)
  
  S_list <- vector("list", Lmax + 1)
  
  for (ell in 0:Lmax) {
    S <- matrix(0, nrow = m, ncol = k)
    rows <- which(lags == ell)
    
    if (length(rows) > 0) {
      for (r in rows) {
        S[r, match(bases[r], variables)] <- 1
      }
    }
    
    S_list[[ell + 1]] <- S
  }
  
  list(S_list = S_list, Lmax = Lmax)
}

build_G_hq <- function(hh, q, S_list, Lmax, Psi_h_all) {
  m <- nrow(S_list[[1]])
  k <- ncol(S_list[[1]])
  G <- matrix(0.0, nrow = m, ncol = k)
  
  for (ell in 0:Lmax) {
    idx <- (hh - ell) - q
    if (idx >= 0) {
      G <- G + S_list[[ell + 1]] %*% Psi_h_all[idx + 1, , ]
    }
  }
  
  G
}

inject_girf_into_Z <- function(Psi_draws, Sigma_kept, terms_df,
                               impulse_idx = 1, shock_scale = 1) {
  H <- dim(Psi_draws)[1] - 1
  k <- dim(Psi_draws)[2]
  M <- dim(Psi_draws)[4]
  
  psiZ <- matrix(0.0, nrow = H + 1, ncol = M)
  
  for (m in seq_len(M)) {
    Sigma_m <- Sigma_kept[, , m]
    Sigma_m <- 0.5 * (Sigma_m + t(Sigma_m))
    
    sigma_jj <- Sigma_m[impulse_idx, impulse_idx]
    if (!is.finite(sigma_jj) || sigma_jj <= 0) next
    
    delta <- shock_scale * Sigma_m[, impulse_idx] / sqrt(sigma_jj)
    Psi_m <- Psi_draws[, , , m, drop = FALSE][, , , 1]
    
    girfY <- matrix(0.0, nrow = H + 1, ncol = k)
    
    for (hh in 0:H) {
      girfY[hh + 1, ] <- Psi_m[hh + 1, , ] %*% delta
    }
    
    acc <- numeric(H + 1)
    
    for (r in seq_len(nrow(terms_df))) {
      lag_r <- terms_df$lag[r]
      col_r <- terms_df$col[r]
      beta_r <- terms_df$beta[r]
      
      contrib <- if (lag_r == 0) {
        girfY[, col_r]
      } else {
        c(rep(0, lag_r), girfY[1:(H + 1 - lag_r), col_r])
      }
      
      acc <- acc + beta_r * contrib
    }
    
    psiZ[, m] <- acc
  }
  
  qs <- c(0.05, 0.16, 0.50, 0.84, 0.95)
  qmat <- t(apply(psiZ, 1, quantile, probs = qs, na.rm = TRUE))
  
  bands <- data.table(
    horizon = 0:H,
    lower90 = qmat[, 1],
    lower68 = qmat[, 2],
    median = qmat[, 3],
    upper68 = qmat[, 4],
    upper90 = qmat[, 5]
  )
  
  list(draws = psiZ, bands = bands)
}

compute_mu_baseline_draws <- function(DT_var, variables, terms_df, beta0, B_kept, p, H) {
  Tn <- nrow(DT_var)
  M <- dim(B_kept)[3]
  Lmax <- max(terms_df$lag)
  
  Y_all <- as.matrix(DT_var[, ..variables])
  Y_hist <- Y_all[(Tn - p + 1):Tn, , drop = FALSE]
  hist_tail <- Y_all[(Tn - Lmax):Tn, , drop = FALSE]
  
  mu_draws <- matrix(NA_real_, nrow = H + 1, ncol = M)
  
  for (m in seq_len(M)) {
    Y_fore <- forecast_baseline_path(B_kept[, , m], Y_hist, p, H)
    mu_h <- rep(beta0, H + 1)
    
    for (r in seq_len(nrow(terms_df))) {
      base_j <- terms_df$col[r]
      lag_r <- terms_df$lag[r]
      beta_r <- terms_df$beta[r]
      
      for (hh in 0:H) {
        val <- if (hh >= lag_r) {
          Y_fore[hh - lag_r + 1, base_j]
        } else {
          j <- lag_r - hh
          hist_tail[(Lmax + 1) - j, base_j]
        }
        
        mu_h[hh + 1] <- mu_h[hh + 1] + beta_r * val
      }
    }
    
    mu_draws[, m] <- mu_h
  }
  
  mu_draws
}

compute_s2_one_draw <- function(Psi_m, Sigma_m, terms_df, beta_vec,
                                sigma_eta2, impulse_idx = 1) {
  H <- dim(Psi_m)[1] - 1
  variables <- dimnames(Psi_m)$variable
  
  Sigma_m <- 0.5 * (Sigma_m + t(Sigma_m))
  sigma_jj <- Sigma_m[impulse_idx, impulse_idx]
  
  Sigma_cond <- Sigma_m - tcrossprod(Sigma_m[, impulse_idx], Sigma_m[, impulse_idx]) / sigma_jj
  Sigma_cond <- 0.5 * (Sigma_cond + t(Sigma_cond))
  
  sel <- build_selection_from_terms(terms_df, variables)
  S_list <- sel$S_list
  Lmax <- sel$Lmax
  
  m <- nrow(S_list[[1]])
  s2 <- numeric(H + 1)
  s2_delta <- numeric(H + 1)
  
  for (hh in 0:H) {
    Sigma_s <- matrix(0.0, m, m)
    Sigma_sd <- matrix(0.0, m, m)
    
    for (q in 0:hh) {
      G <- build_G_hq(hh, q, S_list, Lmax, Psi_m)
      Sigma_s <- Sigma_s + G %*% Sigma_m %*% t(G)
      
      if (q == 0) {
        Sigma_sd <- Sigma_sd + G %*% Sigma_cond %*% t(G)
      } else {
        Sigma_sd <- Sigma_sd + G %*% Sigma_m %*% t(G)
      }
    }
    
    s2[hh + 1] <- as.numeric(t(beta_vec) %*% Sigma_s %*% beta_vec) + sigma_eta2
    s2_delta[hh + 1] <- as.numeric(t(beta_vec) %*% Sigma_sd %*% beta_vec) + sigma_eta2
  }
  
  list(s2 = s2, s2_delta = s2_delta)
}

compute_factor_moment_draws <- function(kernel, terms_df, beta0, beta_vec,
                                        sigma_eta2, impulse_idx, p_lags, H) {
  variables <- kernel$variables
  M <- dim(kernel$Psi_draws)[4]
  
  mu_draws <- compute_mu_baseline_draws(
    DT_var = as.data.table(kernel$DT),
    variables = variables,
    terms_df = terms_df,
    beta0 = beta0,
    B_kept = kernel$B_kept,
    p = p_lags,
    H = H
  )
  
  s2_draws <- matrix(NA_real_, nrow = H + 1, ncol = M)
  s2_delta_draws <- matrix(NA_real_, nrow = H + 1, ncol = M)
  
  for (m in seq_len(M)) {
    Psi_m <- kernel$Psi_draws[, , , m, drop = FALSE][, , , 1]
    Sigma_m <- kernel$Sigma_kept[, , m]
    
    s2_m <- compute_s2_one_draw(
      Psi_m = Psi_m,
      Sigma_m = Sigma_m,
      terms_df = terms_df,
      beta_vec = beta_vec,
      sigma_eta2 = sigma_eta2,
      impulse_idx = impulse_idx
    )
    
    s2_draws[, m] <- pmax(s2_m$s2, 0)
    s2_delta_draws[, m] <- pmax(s2_m$s2_delta, 0)
  }
  
  list(
    mu_draws = mu_draws,
    s2_draws = s2_draws,
    s2_delta_draws = s2_delta_draws
  )
}

compute_pd_girf <- function(psiZ_draws, mu_draws, s2_draws, s2_delta_draws, p, rho) {
  H <- nrow(psiZ_draws) - 1
  M <- ncol(psiZ_draws)
  qpi <- qnorm(p)
  
  out <- matrix(NA_real_, nrow = H + 1, ncol = M)
  
  for (m in seq_len(M)) {
    mu_base <- mu_draws[, m]
    mu_shock <- mu_base + psiZ_draws[, m]
    
    s_base <- sqrt((1 - rho) + rho * pmax(s2_draws[, m], 0))
    s_shock <- sqrt((1 - rho) + rho * pmax(s2_delta_draws[, m], 0))
    
    pd_base <- pnorm((qpi - sqrt(rho) * mu_base) / s_base)
    pd_shock <- pnorm((qpi - sqrt(rho) * mu_shock) / s_shock)
    
    out[, m] <- pd_shock - pd_base
  }
  
  qs <- c(0.05, 0.16, 0.50, 0.84, 0.95)
  qmat <- t(apply(out, 1, quantile, probs = qs, na.rm = TRUE))
  
  bands <- data.table(
    horizon = 0:H,
    lower90 = qmat[, 1],
    lower68 = qmat[, 2],
    median = qmat[, 3],
    upper68 = qmat[, 4],
    upper90 = qmat[, 5],
    p = p,
    rho = rho
  )
  
  list(draws = out, bands = bands)
}

# -----------------------------------------------------------------------------
# 6. Plotting
# -----------------------------------------------------------------------------

theme_robustness <- function() {
  theme_minimal(base_size = 12) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      axis.title = element_text(color = "#1A1A1A"),
      axis.text = element_text(color = "grey20"),
      strip.text = element_text(face = "bold", color = "#1A1A1A", size = 11),
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(linetype = "dotted", linewidth = 0.25),
      panel.grid.major.y = element_line(linetype = "dotted", linewidth = 0.25)
    )
}

plot_bands_facets <- function(bands, ylab_txt = "Response", var_labels = NULL) {
  dt <- as.data.table(copy(bands))

  if ("variable" %in% names(dt)) {
    dt[, variable_label := as.character(variable)]

    if (!is.null(var_labels)) {
      dt[, variable_label := unname(var_labels[variable_label])]
      dt[is.na(variable_label), variable_label := as.character(variable)]
    }

    # Order facets by the VAR variable ordering (impulse variable first),
    # instead of ggplot's default alphabetical ordering.
    dt[, variable_label := factor(variable_label, levels = unique(variable_label))]
  }
  
  ggplot(dt, aes(x = horizon, y = median)) +
    geom_ribbon(aes(ymin = lower90, ymax = upper90), fill = "#F3D6E3", alpha = 0.65) +
    geom_ribbon(aes(ymin = lower68, ymax = upper68), fill = "#AB4A7D", alpha = 0.42) +
    geom_line(color = "#6F1732", linewidth = 1.10) +
    geom_hline(yintercept = 0, color = "grey45", linewidth = 0.35, linetype = "dashed") +
    facet_wrap(~ variable_label, scales = "free_y") +
    scale_x_continuous(
      breaks = function(x) seq(max(0, ceiling(x[1])), floor(x[2]), by = 1),
      minor_breaks = NULL
    ) +
    scale_y_continuous(labels = label_number(big.mark = " ")) +
    labs(x = "Horizon (quarters)", y = ylab_txt) +
    theme_robustness()
}

plot_single_bands <- function(bands, ylab_txt = "Response") {
  dt <- as.data.table(copy(bands))
  
  ggplot(dt, aes(x = horizon, y = median)) +
    geom_ribbon(aes(ymin = lower90, ymax = upper90), fill = "#F3D6E3", alpha = 0.65) +
    geom_ribbon(aes(ymin = lower68, ymax = upper68), fill = "#AB4A7D", alpha = 0.42) +
    geom_line(color = "#6F1732", linewidth = 1.10) +
    geom_hline(yintercept = 0, color = "grey45", linewidth = 0.35, linetype = "dashed") +
    scale_x_continuous(
      breaks = function(x) seq(max(0, ceiling(x[1])), floor(x[2]), by = 1),
      minor_breaks = NULL
    ) +
    scale_y_continuous(labels = label_number(big.mark = " ")) +
    labs(x = "Horizon (quarters)", y = ylab_txt) +
    theme_robustness()
}

plot_timeseries <- function(dt, xcol = "Date", ycol = "value",
                            ylab_txt = "Value", scale_y = 1) {
  d <- data.table(
    xx = as.Date(dt[[xcol]]),
    yy = as.numeric(dt[[ycol]]) * scale_y
  )
  
  ggplot(d, aes(x = xx, y = yy)) +
    geom_line(color = "#6F1732", linewidth = 0.9) +
    labs(x = "Date", y = ylab_txt) +
    theme_robustness()
}

# PNG only (PDF generation intentionally removed). Name kept for compatibility.
ggsave_both <- function(plot, dir, stem, width, height, dpi = 300) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(dir, paste0(stem, ".png")), plot, width = width, height = height, dpi = dpi)
  # Vector PDF twin (preferred for the paper). grDevices::pdf is used rather
  # than cairo_pdf because the latter is not available in every R build.
  tryCatch(
    ggsave(file.path(dir, paste0(stem, ".pdf")), plot,
           width = width, height = height, device = grDevices::pdf),
    error = function(e) message("PDF export failed for ", stem, ": ",
                                conditionMessage(e))
  )
}

# -----------------------------------------------------------------------------
# 7. LaTeX writers
# -----------------------------------------------------------------------------

stars_factory <- function(breaks = c(0.01, 0.05, 0.10),
                          symbols = c("$^{***}$", "$^{**}$", "$^{*}$")) {
  function(p) {
    ifelse(
      p < breaks[1], symbols[1],
      ifelse(
        p < breaks[2], symbols[2],
        ifelse(p < breaks[3], symbols[3], "")
      )
    )
  }
}

write_stability_table_tex <- function(dt, path, caption, label) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  
  cat("\\begin{table}[!htbp]\n\\centering\n", file = con)
  cat("\\caption{", caption, "}\n", sep = "", file = con)
  cat("\\label{", label, "}\n", sep = "", file = con)
  cat("\\begin{tabular}{lrrrrrrr}\n", file = con)
  cat("\\toprule\n", file = con)
  cat("Spec & $k$ & $p$ & Draws & Stable share & Stable draws & $|\\lambda|_{\\max}$ median & $|\\lambda|_{\\max}$ p95 \\\\\n", file = con)
  cat("\\midrule\n", file = con)
  
  for (i in seq_len(nrow(dt))) {
    cat(
      dt$spec[i], " & ", dt$k[i], " & ", dt$p[i], " & ",
      format(dt$n_draws[i], big.mark = " "), " & ",
      fmt_num(dt$stable_share[i], 4), " & ",
      format(dt$n_stable[i], big.mark = " "), " & ",
      fmt_num(dt$lambda_max_median_stable[i], 4), " & ",
      fmt_num(dt$lambda_max_p95_stable[i], 4), " \\\\\n",
      sep = "", file = con
    )
  }
  
  cat("\\bottomrule\n\\end{tabular}\n", file = con)
  cat("\\begin{tablenotes}\n\\small\n", file = con)
  cat(
    "\\item Notes: $k$ is the number of variables in the VAR (including the GPR index) and $p$ the lag order. ",
    "Stable share is the fraction of posterior draws whose companion-matrix spectral radius is strictly below one. ",
    "The last two columns report the posterior median and 95th percentile of the largest-modulus eigenvalue over stable draws.\n",
    sep = "", file = con
  )
  cat("\\end{tablenotes}\n\\end{table}\n", file = con)
}

write_satellite_table_tex <- function(model, path, caption, label,
                                      extra_rows = list(),
                                      notes = NULL,
                                      stars = stars_factory()) {
  sm <- summary(model)
  co <- as.data.frame(sm$coefficients)
  co$term <- rownames(co)
  rownames(co) <- NULL
  
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  
  cat("\\begin{table}[!htbp]\n\\centering\n", file = con)
  cat("\\caption{", caption, "}\n", sep = "", file = con)
  cat("\\label{", label, "}\n", sep = "", file = con)
  cat("\\begin{tabular}{lrr}\n\\toprule\n", file = con)
  cat("Variable & Estimate & Std. Error \\\\\n\\midrule\n", file = con)
  
  for (i in seq_len(nrow(co))) {
    cat(
      co$term[i], " & ",
      sprintf("%.6f", co$Estimate[i]), stars(co$`Pr(>|t|)`[i]), " & ",
      sprintf("%.6f", co$`Std. Error`[i]), " \\\\\n",
      sep = "", file = con
    )
  }
  
  cat("\\midrule\n", file = con)
  cat("Observations & \\multicolumn{2}{r}{", nobs(model), "} \\\\\n", sep = "", file = con)
  cat("Residual std. error & \\multicolumn{2}{r}{", sprintf("%.4f", sm$sigma), "} \\\\\n", sep = "", file = con)
  cat("$R^2$ / adj. $R^2$ & \\multicolumn{2}{r}{",
      sprintf("%.4f / %.4f", sm$r.squared, sm$adj.r.squared), "} \\\\\n", sep = "", file = con)
  cat("AIC & \\multicolumn{2}{r}{", sprintf("%.3f", AIC(model)), "} \\\\\n", sep = "", file = con)
  
  for (nm in names(extra_rows)) {
    cat(nm, " & \\multicolumn{2}{r}{", extra_rows[[nm]], "} \\\\\n", sep = "", file = con)
  }
  
  cat("\\bottomrule\n\\end{tabular}\n", file = con)
  cat("\\begin{tablenotes}\n\\small\n", file = con)
  
  if (is.null(notes)) {
    notes <- "Notes: The dependent variable is the reconstructed and standardized systematic factor $Z$. Significance: $^{*}p<0.10$, $^{**}p<0.05$, $^{***}p<0.01$."
  }
  
  cat("\\item ", notes, "\n", sep = "", file = con)
  cat("\\end{tablenotes}\n\\end{table}\n", file = con)
}

write_coeftest_table_tex <- function(ct, path, caption, label,
                                     notes = NULL,
                                     stars = stars_factory(
                                       breaks = c(0.001, 0.01, 0.05),
                                       symbols = c("$^{***}$", "$^{**}$", "$^{*}$")
                                     )) {
  m <- as.matrix(ct)
  est <- m[, "Estimate"]
  se  <- m[, "Std. Error"]
  pv  <- m[, ncol(m)]
  terms <- rownames(m)
  
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  
  cat("\\begin{table}[!htbp]\n\\centering\n", file = con)
  cat("\\caption{", caption, "}\n", sep = "", file = con)
  cat("\\label{", label, "}\n", sep = "", file = con)
  cat("\\begin{tabular}{lrr}\n\\toprule\n", file = con)
  cat("Variable & Estimate & Std. Error \\\\\n\\midrule\n", file = con)
  
  for (i in seq_along(terms)) {
    cat(
      terms[i], " & ",
      sprintf("%.6f", est[i]), stars(pv[i]), " & ",
      sprintf("%.6f", se[i]), " \\\\\n",
      sep = "", file = con
    )
  }
  
  cat("\\bottomrule\n\\end{tabular}\n", file = con)
  cat("\\begin{tablenotes}\n\\small\n", file = con)
  
  if (is.null(notes)) {
    notes <- "Notes: OLS point estimates with robust standard errors. Significance: $^{*}p<0.05$, $^{**}p<0.01$, $^{***}p<0.001$."
  }
  
  cat("\\item ", notes, "\n", sep = "", file = con)
  cat("\\end{tablenotes}\n\\end{table}\n", file = con)
}

write_regularized_table_tex <- function(coef_dt, footer_rows, path, caption, label,
                                        notes = NULL) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  
  cat("\\begin{table}[!htbp]\n\\centering\n", file = con)
  cat("\\caption{", caption, "}\n", sep = "", file = con)
  cat("\\label{", label, "}\n", sep = "", file = con)
  cat("\\begin{tabular}{lr}\n\\toprule\n", file = con)
  cat("Variable & Estimate \\\\\n\\midrule\n", file = con)
  
  for (i in seq_len(nrow(coef_dt))) {
    cat(
      coef_dt$variable[i], " & ",
      sprintf("%.5f", coef_dt$coef[i]), " \\\\\n",
      sep = "", file = con
    )
  }
  
  cat("\\midrule\n", file = con)
  
  for (nm in names(footer_rows)) {
    cat(nm, " & ", footer_rows[[nm]], " \\\\\n", sep = "", file = con)
  }
  
  cat("\\bottomrule\n\\end{tabular}\n", file = con)
  cat("\\begin{tablenotes}\n\\small\n", file = con)
  
  if (is.null(notes)) {
    notes <- "Notes: Coefficients from the regularized regression, ordered by absolute value. The optimal penalty is selected by cross-validation (mean squared error)."
  }
  
  cat("\\item ", notes, "\n", sep = "", file = con)
  cat("\\end{tablenotes}\n\\end{table}\n", file = con)
}

write_variable_definitions_tex <- function(def_dt, path, caption, label) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  
  cat("\\begin{table}[!htbp]\n\\centering\n", file = con)
  cat("\\caption{", caption, "}\n", sep = "", file = con)
  cat("\\label{", label, "}\n", sep = "", file = con)
  cat("\\begin{tabular}{lll}\n\\toprule\n", file = con)
  cat("Code & Definition & Transform / Unit \\\\\n", file = con)
  
  groups <- unique(def_dt$group)
  
  for (g in groups) {
    cat("\\midrule\n\\multicolumn{3}{l}{\\textit{", g, "}} \\\\\n", sep = "", file = con)
    sub <- def_dt[group == g]
    
    for (i in seq_len(nrow(sub))) {
      cat(
        "\\texttt{", sub$code[i], "} & ",
        sub$definition[i], " & ",
        sub$transform[i],
        " \\\\\n",
        sep = "", file = con
      )
    }
  }
  
  cat("\\bottomrule\n\\end{tabular}\n", file = con)
  cat("\\begin{tablenotes}\n\\small\n", file = con)
  cat(
    "\\item Notes: ``Real'' variables are deflated by CPI (CPIAUCSL); per-capita measures divide by the civilian noninstitutional population aged 16+ (CNP16OV). Higher-frequency series are averaged to calendar quarters. ``pp'' = percentage points.\n",
    file = con
  )
  cat("\\end{tablenotes}\n\\end{table}\n", file = con)
}

# =============================================================================
# 7. Bayesian satellite + BVAR/regression variance decomposition
# -----------------------------------------------------------------------------
# These helpers mirror, for the DRALACBN application, the EBA pipeline in
# code/04_girf.R. The OLS satellite selected by select_satellite_oos() is used
# ONLY to fix the design; conditional on that design, inference uses the
# conjugate Bayesian linear regression (Jeffreys prior) from
# code/_bayesian_satellite_utils.R. Coefficient and residual-variance
# uncertainty are then propagated, jointly with the BVAR posterior draws, into
# psi_Z and PD. The decomposition splits the psi_Z and PD response variance into
# a BVAR-only part (satellite fixed at its posterior mean), a satellite/
# regression-only part (BVAR fixed at its posterior-mean kernel) and an
# interaction/non-linearity residual.
# =============================================================================

# Build a Bayesian satellite object on the design selected for the OLS model.
make_bayesian_satellite <- function(sat_model, model_df, variables,
                                    n_draws = 10000L, seed = 12345) {
  f <- stats::formula(sat_model)
  sat <- fit_bayesian_lm_jeffreys(
    model_df = as.data.frame(model_df),
    formula_obj = f,
    n_draws = n_draws,
    seed = seed
  )

  parsed <- parse_satellite_terms(sat_model, variables)
  terms_df <- parsed$terms

  missing_cols <- setdiff(terms_df$term, colnames(sat$beta_draws))
  if (length(missing_cols) > 0L) {
    stop("Bayesian satellite draws are missing columns: ", paste(missing_cols, collapse = ", "))
  }

  list(
    terms = terms_df,
    beta0_ols = parsed$beta0,
    beta0_draws = sat$beta_draws[, "(Intercept)"],
    beta_term_draws = sat$beta_draws[, terms_df$term, drop = FALSE],
    beta_ols = sat$beta_ols,
    sigma2_draws = sat$sigma2_draws,
    sigma2_ols = sat$sigma2_ols_unbiased,
    n_draws = nrow(sat$beta_draws),
    prior = sat$prior,
    full = sat
  )
}


# psi_Z(h) with satellite coefficients varying by posterior draw (Bayesian).
inject_girf_into_Z_bayes <- function(Psi_draws, Sigma_kept, bsat,
                                     impulse_idx = 1, shock_scale = 1) {
  H <- dim(Psi_draws)[1] - 1
  k <- dim(Psi_draws)[2]
  M <- dim(Psi_draws)[4]

  terms_df <- bsat$terms
  beta_term_draws <- bsat$beta_term_draws
  n_sat <- nrow(beta_term_draws)

  psiZ <- matrix(0.0, nrow = H + 1, ncol = M)

  for (m in seq_len(M)) {
    sm <- pair_satellite_draw(m, n_sat)
    beta_m <- as.numeric(beta_term_draws[sm, terms_df$term, drop = TRUE])

    Sigma_m <- Sigma_kept[, , m]
    Sigma_m <- 0.5 * (Sigma_m + t(Sigma_m))
    sigma_jj <- Sigma_m[impulse_idx, impulse_idx]
    if (!is.finite(sigma_jj) || sigma_jj <= 0) next

    delta <- shock_scale * Sigma_m[, impulse_idx] / sqrt(sigma_jj)
    Psi_m <- Psi_draws[, , , m, drop = FALSE][, , , 1]

    girfY <- matrix(0.0, nrow = H + 1, ncol = k)
    for (hh in 0:H) girfY[hh + 1, ] <- Psi_m[hh + 1, , ] %*% delta

    acc <- numeric(H + 1)
    for (r in seq_len(nrow(terms_df))) {
      lag_r <- terms_df$lag[r]
      col_r <- terms_df$col[r]
      contrib <- if (lag_r == 0) girfY[, col_r] else c(rep(0, lag_r), girfY[1:(H + 1 - lag_r), col_r])
      acc <- acc + beta_m[r] * contrib
    }
    psiZ[, m] <- acc
  }

  qs <- c(0.05, 0.16, 0.50, 0.84, 0.95)
  qmat <- t(apply(psiZ, 1, quantile, probs = qs, na.rm = TRUE))
  bands <- data.table(
    horizon = 0:H,
    lower90 = qmat[, 1], lower68 = qmat[, 2], median = qmat[, 3],
    upper68 = qmat[, 4], upper90 = qmat[, 5]
  )

  list(draws = psiZ, bands = bands)
}

# Baseline mean mu(h) and conditional variances s2(h), satellite varying by draw.
compute_factor_moment_draws_bayes <- function(kernel, bsat, impulse_idx, p_lags, H) {
  variables <- kernel$variables
  M <- dim(kernel$Psi_draws)[4]

  terms_df <- bsat$terms
  beta_term_draws <- bsat$beta_term_draws
  beta0_draws <- bsat$beta0_draws
  sigma2_draws <- bsat$sigma2_draws
  n_sat <- nrow(beta_term_draws)

  DT_var <- as.data.table(kernel$DT)
  Tn <- nrow(DT_var)
  Lmax <- max(terms_df$lag)
  Y_all <- as.matrix(DT_var[, ..variables])
  Y_hist <- Y_all[(Tn - p_lags + 1):Tn, , drop = FALSE]
  hist_tail <- Y_all[(Tn - Lmax):Tn, , drop = FALSE]

  mu_draws <- matrix(NA_real_, nrow = H + 1, ncol = M)
  s2_draws <- matrix(NA_real_, nrow = H + 1, ncol = M)
  s2_delta_draws <- matrix(NA_real_, nrow = H + 1, ncol = M)

  for (m in seq_len(M)) {
    sm <- pair_satellite_draw(m, n_sat)
    beta_m <- as.numeric(beta_term_draws[sm, terms_df$term, drop = TRUE])
    beta0_m <- beta0_draws[sm]
    sigma2_m <- sigma2_draws[sm]

    Y_fore <- forecast_baseline_path(kernel$B_kept[, , m], Y_hist, p_lags, H)
    mu_h <- rep(beta0_m, H + 1)
    for (r in seq_len(nrow(terms_df))) {
      base_j <- terms_df$col[r]
      lag_r <- terms_df$lag[r]
      for (hh in 0:H) {
        val <- if (hh >= lag_r) Y_fore[hh - lag_r + 1, base_j] else hist_tail[(Lmax + 1) - (lag_r - hh), base_j]
        mu_h[hh + 1] <- mu_h[hh + 1] + beta_m[r] * val
      }
    }
    mu_draws[, m] <- mu_h

    Psi_m <- kernel$Psi_draws[, , , m, drop = FALSE][, , , 1]
    Sigma_m <- kernel$Sigma_kept[, , m]
    s2m <- compute_s2_one_draw(
      Psi_m = Psi_m, Sigma_m = Sigma_m, terms_df = terms_df,
      beta_vec = beta_m, sigma_eta2 = sigma2_m, impulse_idx = impulse_idx
    )
    s2_draws[, m] <- pmax(s2m$s2, 0)
    s2_delta_draws[, m] <- pmax(s2m$s2_delta, 0)
  }

  list(mu_draws = mu_draws, s2_draws = s2_draws, s2_delta_draws = s2_delta_draws)
}


# -----------------------------------------------------------------------------
# Plotting: colour-blind-safe variance-decomposition share plot.
# Colours stay in the burgundy/rose house style; the share plot adds redundant
# line type + point shape so it reads under colour-vision deficiency / greyscale.
# -----------------------------------------------------------------------------

decomp_source_colors_robust <- c(
  "BVAR" = "#6F1732",
  "Satellite (regression)" = "#AB4A7D",
  "Interaction / non-linearity" = "#6B6B6B"
)
decomp_source_linetypes_robust <- c(
  "BVAR" = "solid",
  "Satellite (regression)" = "longdash",
  "Interaction / non-linearity" = "dotted"
)
decomp_source_shapes_robust <- c(
  "BVAR" = 16,
  "Satellite (regression)" = 17,
  "Interaction / non-linearity" = 15
)

plot_variance_decomposition_shares_robust <- function(summary_df) {
  dt <- as.data.table(copy(summary_df))
  long <- melt(
    dt, id.vars = "horizon",
    measure.vars = c("share_bvar", "share_satellite", "share_interaction"),
    variable.name = "source", value.name = "share"
  )
  long[, source := factor(
    source,
    levels = c("share_bvar", "share_satellite", "share_interaction"),
    labels = c("BVAR", "Satellite (regression)", "Interaction / non-linearity")
  )]

  ggplot(long, aes(x = horizon, y = share, color = source, linetype = source, shape = source)) +
    geom_hline(yintercept = 0, color = "grey55", linewidth = 0.30, linetype = "dashed") +
    geom_line(linewidth = 1.05) +
    geom_point(size = 1.8) +
    scale_color_manual(values = decomp_source_colors_robust) +
    scale_linetype_manual(values = decomp_source_linetypes_robust) +
    scale_shape_manual(values = decomp_source_shapes_robust) +
    scale_x_continuous(
      breaks = function(x) seq(max(0, ceiling(x[1])), floor(x[2]), by = 1),
      minor_breaks = NULL
    ) +
    labs(x = "Horizon (quarters)", y = "Share of total variance (%)") +
    theme_robustness()
}

invisible(TRUE)