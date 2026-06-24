source("code/00_setup.R")


# =============================================================================
# Script Z-factor / Merton satellite — USA
#
# Baseline:
#   - rho estimated by imposing Var(Z) = 1.
#
# Robustness:
#   - fixed rho values, including a Basel-style corporate asset correlation.
#
# Important:
#   - NO centering of Z.
#   - NO centering of macro variables for the main model.
#   - A centered-X version of the selected model is estimated ONLY for display.
#
# Main outputs:
#   - output/models/best_model.rds
#   - output/models/best_model_centerX.rds
#   - output/models/by_rho/<model_tag>/best_model.rds
#   - output/models/by_rho/<model_tag>/best_model_centerX.rds
#   - output/zfactor/by_rho/<model_tag>/model_summary.txt
#   - output/zfactor/by_rho/<model_tag>/model_summary_centerX.txt
#   - output/zfactor/merton_params_by_rho_usa.csv
#   - output/zfactor/model_comparison_by_rho.csv
#   - output/zfactor/coefficients_by_rho.csv
# =============================================================================

# =============================================================================
# Packages
# =============================================================================

packages <- c(
  "dplyr",
  "lubridate",
  "parallel",
  "data.table",
  "glmnet",
  "lmtest",
  "sandwich",
  "MASS"
)

missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0L) {
  stop(
    "Missing required packages: ", paste(missing_packages, collapse = ", "),
    "\nRun renv::restore() from the repository root before running this script.",
    call. = FALSE
  )
}

invisible(lapply(packages, library, character.only = TRUE))

# =============================================================================
# Global parameters
# =============================================================================

P_VALUE_THRESHOLD      <- if (exists("P_VALUE_THRESHOLD")) P_VALUE_THRESHOLD else 0.05
NB_LAGS                <- if (exists("NB_LAGS")) NB_LAGS else 4
MAX_VARIABLES_IN_MODEL <- if (exists("MAX_VARIABLES_IN_MODEL")) MAX_VARIABLES_IN_MODEL else 4

constraint_sign <- FALSE

RUN_REGULARIZATION_FOR_BASELINE  <- TRUE
RUN_REGULARIZATION_FOR_FIXED_RHO <- FALSE

RUN_RIDGE  <- TRUE
RUN_LASSO  <- TRUE
USE_LAMBDA <- "lambda.1se"

CORE_FRACTION <- 0.75

set.seed(if (exists("SEED_ZFACTOR")) SEED_ZFACTOR else 12345)

# =============================================================================
# Output directories
# =============================================================================

output_dir <- if (exists("DIR_ZFACTOR")) DIR_ZFACTOR else "output/zfactor"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

models_dir <- if (exists("DIR_MODELS")) DIR_MODELS else "output/models"
dir.create(models_dir, showWarnings = FALSE, recursive = TRUE)

models_by_rho_dir <- file.path(models_dir, "by_rho")
dir.create(models_by_rho_dir, showWarnings = FALSE, recursive = TRUE)

z_by_rho_dir <- file.path(output_dir, "by_rho")
dir.create(z_by_rho_dir, showWarnings = FALSE, recursive = TRUE)

selection_by_rho_dir <- file.path(output_dir, "model_selection_by_rho")
dir.create(selection_by_rho_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# Helper functions
# =============================================================================

get_num_cores <- function(fraction = CORE_FRACTION) {
  max(1L, floor(parallel::detectCores() * fraction))
}

rho_tag_from_value <- function(rho) {
  paste0("rho_", gsub("\\.", "p", sprintf("%.4f", rho)))
}

fmt_dur <- function(sec) {
  sec <- max(0, as.numeric(sec))
  h <- floor(sec / 3600)
  m <- floor((sec - 3600 * h) / 60)
  s <- round(sec - 3600 * h - 60 * m)
  sprintf("%02d:%02d:%02d", h, m, s)
}

drop_zero_variance <- function(X, tol = 1e-12) {
  sds <- apply(X, 2, sd, na.rm = TRUE)
  keep <- is.finite(sds) & sds > tol
  
  list(
    X = X[, keep, drop = FALSE],
    dropped = colnames(X)[!keep],
    kept = colnames(X)[keep]
  )
}

vif_from_matrix <- function(X) {
  p <- ncol(X)
  out <- numeric(p)
  
  for (j in seq_len(p)) {
    y <- X[, j]
    Z <- X[, -j, drop = FALSE]
    mod <- lm(y ~ Z)
    r2 <- summary(mod)$r.squared
    
    out[j] <- if (is.finite(r2) && r2 < 1) {
      1 / (1 - r2)
    } else {
      Inf
    }
  }
  
  setNames(out, colnames(X))
}

top_abs_cor_pairs <- function(C, top_k = 40, thresh = 0.95) {
  cn <- colnames(C)
  p <- ncol(C)
  
  if (p < 2) {
    return(list(top = data.frame(), above_thresh = data.frame()))
  }
  
  pairs <- which(upper.tri(C), arr.ind = TRUE)
  
  res <- data.frame(
    i = cn[pairs[, 1]],
    j = cn[pairs[, 2]],
    corr = C[pairs],
    stringsAsFactors = FALSE
  )
  
  res$abs_corr <- abs(res$corr)
  res <- res[order(-res$abs_corr), ]
  
  list(
    top = head(res, top_k),
    above_thresh = subset(res, abs_corr >= thresh)
  )
}

condition_index <- function(X) {
  Xs <- scale(X)
  s <- svd(Xs, nu = 0, nv = 0)$d
  
  list(
    kappa = max(s) / min(s),
    singvals = s
  )
}

make_blocked_foldid <- function(n, K = 5) {
  K <- max(3, min(K, n))
  rep(rep(seq_len(K), each = ceiling(n / K)), length.out = n)
}

basel_corporate_rho <- function(pd) {
  pd <- pmin(pmax(pd, 1e-8), 1 - 1e-8)
  w <- (1 - exp(-50 * pd)) / (1 - exp(-50))
  0.12 * w + 0.24 * (1 - w)
}

# =============================================================================
# Z-factor construction
# =============================================================================

f_Z_estimation <- function(vect_TD, rho_fixed = NULL, eps = 1e-8) {
  
  vect_TD <- as.numeric(vect_TD)
  vect_TD <- pmin(pmax(vect_TD, eps), 1 - eps)
  
  p_ttc <- mean(vect_TD, na.rm = TRUE)
  
  compute_Z <- function(rho) {
    (
      qnorm(p_ttc) -
        qnorm(vect_TD) * sqrt(1 - rho)
    ) / sqrt(rho)
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
      mean_Z_raw = mean(Z_estim, na.rm = TRUE),
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
      mean_Z_raw = NA_real_,
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
    mean_Z_raw = mean(Z_estim, na.rm = TRUE),
    var_Z_raw = var(Z_estim, na.rm = TRUE)
  )
}

# =============================================================================
# Step 1. Load default data and define rho specifications
# =============================================================================

cat("--- Step 1: Computing the systemic factor Z and rho values for the USA ---\n")

data_risk_all <- read.csv("data/raw/data_risk_EBA.csv", header = TRUE, sep = ";")

data_risk_usa <- data_risk_all %>%
  filter(Country == "United States")

dates_risk <- as.Date(data_risk_usa$Date, format = "%d/%m/%Y")

Z_result_baseline <- f_Z_estimation(data_risk_usa$Corpo_DR_WA)

rho_usa   <- Z_result_baseline$rho
p_ttc_usa <- Z_result_baseline$p_ttc
rho_basel <- basel_corporate_rho(p_ttc_usa)

cat(sprintf("Estimated baseline rho = %.8f\n", rho_usa))
cat(sprintf("Estimated TTC PD       = %.8f\n", p_ttc_usa))
cat(sprintf("Basel corporate rho    = %.8f\n", rho_basel))

rho_specs <- data.table(
  model_tag = c(
    "rho_estimated",
    "rho_0p1000",
    "rho_basel_corporate"
  ),
  rho_value = c(
    rho_usa,
    0.10,
    rho_basel
  ),
  rho_source = c(
    "estimated_baseline",
    "fixed_sensitivity",
    "regulatory_basel_corporate"
  )
)

rho_specs[, p_ttc := p_ttc_usa]
rho_specs[, country := "United States"]
rho_specs[, date_min := min(dates_risk, na.rm = TRUE)]
rho_specs[, date_max := max(dates_risk, na.rm = TRUE)]
rho_specs[, n_obs := sum(!is.na(data_risk_usa$Corpo_DR_WA))]
rho_specs[, mean_default_rate := mean(data_risk_usa$Corpo_DR_WA, na.rm = TRUE)]

fwrite(
  rho_specs,
  file.path(output_dir, "merton_params_by_rho_usa.csv")
)

fwrite(
  rho_specs[model_tag == "rho_estimated"],
  file.path(output_dir, "merton_params_usa.csv")
)

cat("Merton parameters saved.\n\n")

# =============================================================================
# Step 2. Prepare macro predictors and lagged design
# =============================================================================

cat("--- Step 2: Preparing macro predictors ---\n")

data_macro_usa <- read.csv(if (exists("PATH_DATA_VAR")) PATH_DATA_VAR else "data/processed/data_var_for_model.csv", header = TRUE)
data_macro_usa$Date <- as.Date(data_macro_usa$Date)

data_macro_usa$Date <- data_macro_usa$Date %m+% months(2)

vars_to_lag <- if (exists("VAR_BASELINE")) {
  setdiff(VAR_BASELINE, if (exists("IMPULSE_NAME")) IMPULSE_NAME else "log_GPRD")
} else {
  c("vix", "log_sp500_real", "log_oil_real", "log_hours_pc", "log_gdp_pc")
}

missing_macro <- setdiff(vars_to_lag, names(data_macro_usa))

if (length(missing_macro) > 0) {
  stop("Missing macro variables: ", paste(missing_macro, collapse = ", "))
}

lag_indices <- 0:NB_LAGS

lagged_vars_list <- lapply(vars_to_lag, function(var_name) {
  sapply(lag_indices, function(k) {
    dplyr::lag(data_macro_usa[[var_name]], n = k)
  })
})

data_lags <- as.data.frame(do.call(cbind, lagged_vars_list))

names(data_lags) <- paste0(
  rep(vars_to_lag, each = length(lag_indices)),
  "_lag",
  lag_indices
)

data_lags$Date <- data_macro_usa$Date
data_lags$Date <- floor_date(data_lags$Date, unit = "month")

cat("Lagged predictors created successfully.\n\n")

# =============================================================================
# Pre-build model combinations once
# =============================================================================

X_names_global <- setdiff(names(data_lags), "Date")

combs_global <- unlist(
  lapply(1:min(MAX_VARIABLES_IN_MODEL, length(X_names_global)), function(k) {
    combn(X_names_global, k, simplify = FALSE)
  }),
  recursive = FALSE
)

total_models_global <- length(combs_global)

cat(sprintf(
  "Number of candidate model combinations per rho: %d\n\n",
  total_models_global
))

# =============================================================================
# Optional sign constraints
# =============================================================================

if (constraint_sign) {
  expected_signs <- list(
    vix = -1,
    log_sp500_real = 1,
    log_oil_real = -1,
    log_hours_pc = 1,
    log_gdp_pc = 1
  )
} else {
  expected_signs <- NULL
}

# =============================================================================
# Collinearity diagnostics
# =============================================================================

run_collinearity_diagnostics <- function(model_df, output_dir) {
  
  cat("--- Collinearity diagnostics on the design matrix ---\n")
  
  X_names <- names(model_df)[!names(model_df) %in% c("Date", "Y")]
  X_raw <- as.matrix(model_df[, X_names, drop = FALSE])
  
  nzv <- drop_zero_variance(X_raw)
  X2 <- nzv$X
  
  if (length(nzv$dropped)) {
    cat("Zero-variance columns removed:",
        paste(nzv$dropped, collapse = ", "), "\n")
  }
  
  Xs <- scale(X2)
  C <- cor(Xs)
  tp <- top_abs_cor_pairs(C, top_k = 40, thresh = 0.8)
  vif <- vif_from_matrix(Xs)
  ci <- condition_index(X2)
  
  cat("Top 10 absolute correlations:\n")
  print(head(tp$top[, c("i", "j", "corr")], 10))
  
  vif_dt <- data.frame(
    variable = names(vif),
    VIF = as.numeric(vif),
    stringsAsFactors = FALSE
  )
  
  vif_dt <- vif_dt[order(-vif_dt$VIF), ]
  
  cat(sprintf("Condition index (SVD kappa) = %.2f\n", ci$kappa))
  
  colld <- file.path(output_dir, "collinearity")
  dir.create(colld, showWarnings = FALSE, recursive = TRUE)
  
  write.csv(vif_dt, file.path(colld, "vif.csv"), row.names = FALSE)
  write.csv(tp$top, file.path(colld, "top_correlations.csv"), row.names = FALSE)
  write.csv(tp$above_thresh, file.path(colld, "corr_above_0p95.csv"), row.names = FALSE)
  write.csv(
    data.frame(singval = ci$singvals),
    file.path(colld, "singular_values.csv"),
    row.names = FALSE
  )
  
  cat("Collinearity diagnostics exported.\n\n")
  
  invisible(list(nzv = nzv, vif = vif_dt, top_cor = tp, condition = ci))
}

# =============================================================================
# OLS model selection
# =============================================================================

select_best_ols_model <- function(
    model_df,
    model_tag,
    rho_value,
    rho_source,
    combs = combs_global,
    total_models = total_models_global
) {
  
  cat("\n============================================================\n")
  cat("OLS selection for ", model_tag, " | rho = ", sprintf("%.8f", rho_value), "\n", sep = "")
  cat("============================================================\n")
  
  cat(sprintf("Testing %d model combinations...\n", total_models))
  
  num_cores <- get_num_cores()
  cat(sprintf("Using %d fork workers (%.0f%% of %d logical cores).\n",
              num_cores, 100 * CORE_FRACTION, parallel::detectCores()))
  
  evaluate_combo <- function(vars) {
    f <- as.formula(paste("Y ~", paste(vars, collapse = " + ")))
    model <- try(lm(f, data = model_df), silent = TRUE)
    
    if (inherits(model, "try-error")) return(NULL)
    
    p_values <- try(summary(model)$coefficients[-1, "Pr(>|t|)"], silent = TRUE)
    
    if (inherits(p_values, "try-error") || any(is.na(p_values))) {
      return(NULL)
    }
    
    if (any(p_values > P_VALUE_THRESHOLD)) {
      return(NULL)
    }
    
    if (!is.null(expected_signs)) {
      coefs <- coef(model)[-1]
      
      for (i in seq_along(coefs)) {
        base_var_name <- sub("_lag[0-9]+$", "", names(coefs)[i])
        want <- expected_signs[[base_var_name]]
        
        if (!is.null(want) && want != 0) {
          if (sign(coefs[i]) != want) return(NULL)
        }
      }
    }
    
    list(
      vars = vars,
      aic = AIC(model),
      bic = BIC(model)
    )
  }
  
  par_lapply <- if (.Platform$OS.type == "windows") {
    function(X, FUN) lapply(X, FUN)
  } else {
    function(X, FUN) parallel::mclapply(X, FUN, mc.cores = num_cores)
  }
  
  BATCH_SIZE <- max(2000, min(25000, ceiling(total_models / (num_cores * 6))))
  idx_all <- seq_len(total_models)
  batches <- split(idx_all, ceiling(idx_all / BATCH_SIZE))
  
  cat(sprintf(
    "Processing %d batch(es) of about %d models...\n",
    length(batches),
    BATCH_SIZE
  ))
  
  flush.console()
  
  start_time <- Sys.time()
  count_done <- 0L
  results <- vector("list", total_models)
  
  for (b in seq_along(batches)) {
    ids <- batches[[b]]
    combs_batch <- combs[ids]
    
    t_b0 <- Sys.time()
    ans <- par_lapply(combs_batch, evaluate_combo)
    t_b1 <- Sys.time()
    
    results[ids] <- ans
    count_done <- count_done + length(ids)
    
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    speed <- if (elapsed > 0) count_done / elapsed else NA_real_
    remaining <- if (is.finite(speed) && speed > 0) {
      (total_models - count_done) / speed
    } else {
      NA_real_
    }
    pct <- 100 * count_done / total_models
    
    cat(sprintf(
      "[%.1f%%] %d/%d done | elapsed: %s | ETA: %s | batch: %.2fs\n",
      pct,
      count_done,
      total_models,
      fmt_dur(elapsed),
      fmt_dur(remaining),
      as.numeric(difftime(t_b1, t_b0, units = "secs"))
    ))
    
    flush.console()
  }
  
  bad <- vapply(results, function(x) inherits(x, "try-error"), logical(1))
  
  if (any(bad)) {
    results[bad] <- list(NULL)
  }
  
  valid_models <- Filter(Negate(is.null), results)
  
  if (length(valid_models) == 0) {
    stop("No model found for ", model_tag, " satisfying all restrictions.")
  }
  
  best_idx <- which.min(sapply(valid_models, `[[`, "aic"))
  best_model_info <- valid_models[[best_idx]]
  
  best_model <- lm(
    as.formula(paste("Y ~", paste(best_model_info$vars, collapse = " + "))),
    data = model_df
  )
  
  valid_models_table <- rbindlist(
    lapply(valid_models, function(x) {
      data.table(
        aic = x$aic,
        bic = x$bic,
        model = paste(x$vars, collapse = " + "),
        n_vars = length(x$vars)
      )
    }),
    use.names = TRUE,
    fill = TRUE
  )
  
  valid_models_table[, delta_aic := aic - min(aic)]
  valid_models_table[, weight_aic := exp(-0.5 * delta_aic)]
  valid_models_table[, weight_aic := weight_aic / sum(weight_aic)]
  valid_models_table <- valid_models_table[order(aic)]
  
  cat(sprintf(
    "Best model for %s | rho=%.8f | AIC=%.5f | valid models=%d\n",
    model_tag,
    rho_value,
    best_model_info$aic,
    nrow(valid_models_table)
  ))
  
  rm(results, valid_models)
  gc(verbose = FALSE)
  
  list(
    best_model = best_model,
    best_model_info = best_model_info,
    valid_models_table = valid_models_table,
    total_models = total_models,
    n_valid_models = nrow(valid_models_table),
    rho_value = rho_value,
    rho_source = rho_source,
    model_tag = model_tag
  )
}

# =============================================================================
# Center-X version for display only
# =============================================================================

make_centerX_model <- function(model_df, best_model) {
  
  vars <- setdiff(names(coef(best_model)), "(Intercept)")
  
  model_df_centerX <- copy(as.data.table(model_df))
  
  x_means <- sapply(vars, function(v) mean(model_df_centerX[[v]], na.rm = TRUE))
  
  for (v in vars) {
    model_df_centerX[[v]] <- model_df_centerX[[v]] - x_means[[v]]
  }
  
  model_df_centerX <- as.data.frame(model_df_centerX)
  
  f_center <- as.formula(paste("Y ~", paste(vars, collapse = " + ")))
  best_model_centerX <- lm(f_center, data = model_df_centerX)
  
  list(
    model = best_model_centerX,
    model_df = model_df_centerX,
    x_means = x_means
  )
}

# =============================================================================
# Regularization
# =============================================================================

run_regularization <- function(model_df, output_dir, best_model, model_tag, rho_value, p_ttc) {
  
  cat("--- Ridge and Lasso robustness for ", model_tag, " ---\n", sep = "")
  
  stopifnot(all(c("Date", "Y") %in% names(model_df)))
  
  y_vec <- model_df$Y
  X_all <- as.matrix(model_df[, setdiff(names(model_df), c("Date", "Y")), drop = FALSE])
  
  nzv2 <- drop_zero_variance(X_all)
  X_mat <- nzv2$X
  kept_cols <- colnames(X_mat)
  
  build_limits <- function(col_names, expected_signs, use_constraints = FALSE) {
    p <- length(col_names)
    ll <- rep(-Inf, p)
    ul <- rep(Inf, p)
    
    if (!use_constraints || is.null(expected_signs)) {
      return(list(ll = ll, ul = ul))
    }
    
    base_of <- function(nm) sub("_lag[0-9]+$", "", nm)
    
    for (j in seq_len(p)) {
      base <- base_of(col_names[j])
      want <- expected_signs[[base]]
      
      if (!is.null(want) && want != 0) {
        if (want == 1) ll[j] <- 0
        if (want == -1) ul[j] <- 0
      }
    }
    
    list(ll = ll, ul = ul)
  }
  
  lims <- build_limits(kept_cols, expected_signs, constraint_sign)
  
  K_FOLDS <- 5
  foldid <- make_blocked_foldid(nrow(X_mat), K = K_FOLDS)
  
  reg_dir <- file.path(output_dir, "regularization")
  dir.create(reg_dir, showWarnings = FALSE, recursive = TRUE)
  
  write_nonzero_coefs <- function(fit_cv, which_lambda, path_csv) {
    cf <- as.matrix(coef(fit_cv, s = which_lambda))
    
    out <- data.frame(
      variable = rownames(cf),
      coef = as.numeric(cf[, 1]),
      row.names = NULL
    )
    
    out$nonzero <- out$coef != 0
    out <- out[order(-abs(out$coef)), ]
    
    write.csv(out, path_csv, row.names = FALSE)
    invisible(out)
  }
  
  compute_fit_stats <- function(y, yhat, df, name) {
    n <- length(y)
    rss <- sum((y - yhat)^2)
    mse <- rss / n
    rmse <- sqrt(mse)
    aic <- n * log(mse) + 2 * df
    bic <- n * log(mse) + log(n) * df
    
    data.frame(
      model = name,
      n = n,
      df = df,
      RMSE = rmse,
      MSE = mse,
      RSS = rss,
      AIC = aic,
      BIC = bic
    )
  }
  
  take_lambda <- function(cvfit, use) {
    if (use == "lambda.min") cvfit$lambda.min else cvfit$lambda.1se
  }
  
  idx_for_lambda <- function(cvfit, lam) {
    which.min(abs(cvfit$lambda - lam))
  }
  
  metrics_list <- list()
  
  if (RUN_RIDGE) {
    cat("• Estimating RIDGE (alpha = 0)...\n")
    
    cv_ridge <- cv.glmnet(
      x = X_mat,
      y = y_vec,
      family = "gaussian",
      alpha = 0,
      standardize = TRUE,
      foldid = foldid,
      lower.limits = lims$ll,
      upper.limits = lims$ul,
      type.measure = "mse",
      nfolds = length(unique(foldid))
    )
    
    saveRDS(cv_ridge, file = file.path(reg_dir, "ridge_cvglmnet.rds"))
    
    png(file.path(reg_dir, "ridge_cv_curve.png"), width = 1000, height = 600)
    plot(cv_ridge)
    title("Ridge CV (MSE)")
    dev.off()
    
    write_nonzero_coefs(
      cv_ridge,
      USE_LAMBDA,
      file.path(reg_dir, paste0("ridge_coefs_", USE_LAMBDA, ".csv"))
    )
    
    yhat_ridge <- as.numeric(predict(cv_ridge, newx = X_mat, s = USE_LAMBDA))
    lam_r <- take_lambda(cv_ridge, USE_LAMBDA)
    idx_r <- idx_for_lambda(cv_ridge, lam_r)
    df_ridge <- cv_ridge$glmnet.fit$df[idx_r]
    
    metrics_list[["ridge"]] <- compute_fit_stats(
      y_vec,
      yhat_ridge,
      df = df_ridge,
      name = paste0("RIDGE_", USE_LAMBDA)
    )
  }
  
  if (RUN_LASSO) {
    cat("• Estimating LASSO (alpha = 1)...\n")
    
    cv_lasso <- cv.glmnet(
      x = X_mat,
      y = y_vec,
      family = "gaussian",
      alpha = 1,
      standardize = TRUE,
      foldid = foldid,
      lower.limits = lims$ll,
      upper.limits = lims$ul,
      type.measure = "mse",
      nfolds = length(unique(foldid))
    )
    
    saveRDS(cv_lasso, file = file.path(reg_dir, "lasso_cvglmnet.rds"))
    
    png(file.path(reg_dir, "lasso_cv_curve.png"), width = 1000, height = 600)
    plot(cv_lasso)
    title("Lasso CV (MSE)")
    dev.off()
    
    write_nonzero_coefs(
      cv_lasso,
      USE_LAMBDA,
      file.path(reg_dir, paste0("lasso_coefs_", USE_LAMBDA, ".csv"))
    )
    
    yhat_lasso <- as.numeric(predict(cv_lasso, newx = X_mat, s = USE_LAMBDA))
    lam_l <- take_lambda(cv_lasso, USE_LAMBDA)
    idx_l <- idx_for_lambda(cv_lasso, lam_l)
    df_lasso <- cv_lasso$glmnet.fit$df[idx_l]
    
    metrics_list[["lasso"]] <- compute_fit_stats(
      y_vec,
      yhat_lasso,
      df = df_lasso,
      name = paste0("LASSO_", USE_LAMBDA)
    )
  }
  
  if (length(metrics_list)) {
    metrics <- do.call(rbind, metrics_list)
    write.csv(metrics, file.path(reg_dir, "regularization_metrics.csv"), row.names = FALSE)
    print(metrics)
  }
  
  txt_path <- file.path(reg_dir, "regularization_summary.txt")
  con <- file(txt_path, open = "wt")
  
  wln <- function(...) {
    cat(paste0(..., collapse = ""), "\n", file = con, append = TRUE)
  }
  
  wln("Regularized regressions summary — Z-factor USA")
  wln("Generated on: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
  wln("")
  wln("Model tag: ", model_tag)
  wln("p_ttc = ", sprintf("%.8f", p_ttc))
  wln("rho   = ", sprintf("%.8f", rho_value))
  wln("")
  wln("NB_LAGS = ", NB_LAGS)
  wln("MAX_VARIABLES_IN_MODEL = ", MAX_VARIABLES_IN_MODEL)
  wln("USE_LAMBDA = ", USE_LAMBDA)
  wln("")
  
  yhat_ols <- as.numeric(predict(best_model, newdata = model_df))
  rmse_ols <- sqrt(mean((model_df$Y - yhat_ols)^2))
  
  wln("=== OLS selected model ===")
  wln(
    "AIC = ", sprintf("%.3f", AIC(best_model)),
    " | BIC = ", sprintf("%.3f", BIC(best_model)),
    " | in-sample RMSE = ", sprintf("%.5f", rmse_ols)
  )
  wln("OLS variables: ", paste(setdiff(names(coef(best_model)), "(Intercept)"), collapse = ", "))
  wln("")
  
  close(con)
  
  cat("Ridge/Lasso outputs written to ", reg_dir, ".\n\n", sep = "")
  
  invisible(TRUE)
}

# =============================================================================
# Post-OLS outputs
# =============================================================================

run_post_ols_outputs <- function(
    model_df,
    best_model,
    best_model_centerX,
    centerX_means,
    model_info,
    out_dir,
    model_tag,
    rho_value,
    p_ttc,
    mean_Z_raw,
    var_Z_raw
) {
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Main model summary.
  summary_path <- file.path(out_dir, "model_summary.txt")
  
  sink(summary_path)
  cat("Selected Z-factor satellite model — USA\n\n")
  cat("Model type: raw model used for GIRF/PD pipeline\n")
  cat("model_tag =", model_tag, "\n")
  cat("p_ttc     =", p_ttc, "\n")
  cat("rho       =", rho_value, "\n")
  cat("mean_Z_raw =", mean_Z_raw, "\n")
  cat("var_Z_raw  =", var_Z_raw, "\n\n")
  cat("NB_LAGS =", NB_LAGS, " | MAX_VARIABLES_IN_MODEL =", MAX_VARIABLES_IN_MODEL, "\n\n")
  cat("Candidate variables:", paste(vars_to_lag, collapse = ", "), "\n\n")
  cat("Sign constraints:", constraint_sign, "\n\n")
  cat("Selection criterion: lowest AIC among valid models.\n\n")
  cat("Selected model AIC:", model_info$aic, "\n\n")
  print(summary(best_model))
  sink()
  
  # Center-X display model summary.
  summary_center_path <- file.path(out_dir, "model_summary_centerX.txt")
  
  sink(summary_center_path)
  cat("Selected Z-factor satellite model — USA\n\n")
  cat("Model type: centered-X display version only\n")
  cat("Important: this model is NOT used for GIRF/PD pipeline.\n")
  cat("The dependent variable Z is NOT centered.\n")
  cat("Only the selected regressors are demeaned for display.\n\n")
  cat("model_tag =", model_tag, "\n")
  cat("p_ttc     =", p_ttc, "\n")
  cat("rho       =", rho_value, "\n")
  cat("mean_Z_raw =", mean_Z_raw, "\n")
  cat("var_Z_raw  =", var_Z_raw, "\n\n")
  cat("Variables demeaned around the following means:\n")
  print(centerX_means)
  cat("\n")
  print(summary(best_model_centerX))
  sink()
  
  plot_path <- file.path(out_dir, "model_fit_diagnostic.png")
  
  png(plot_path, width = 1000, height = 600)
  fitted_values <- predict(best_model, newdata = model_df)
  plot(
    model_df$Date,
    model_df$Y,
    type = "l",
    col = "grey50",
    lwd = 2,
    ylab = "Systematic factor Z",
    xlab = "Date",
    main = paste0("Z-factor satellite fit — ", model_tag),
    sub = paste("AIC =", round(model_info$aic, 3), "| rho =", round(rho_value, 4))
  )
  lines(model_df$Date, fitted_values, col = "darkgreen", lwd = 2)
  legend(
    "topleft",
    bty = "n",
    legend = c("Observed Z", "Fitted model"),
    col = c("grey50", "darkgreen"),
    lwd = 2,
    lty = 1
  )
  dev.off()
  
  hac <- coeftest(
    best_model,
    vcov = NeweyWest(best_model, lag = NB_LAGS, prewhite = TRUE, adjust = TRUE)
  )
  
  hc3 <- coeftest(
    best_model,
    vcov = vcovHC(best_model, type = "HC3")
  )
  
  write.table(
    capture.output(hac),
    file.path(out_dir, "hac_newey_west.txt"),
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
  
  write.table(
    capture.output(hc3),
    file.path(out_dir, "hc3_robust_se.txt"),
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
  
  # Out-of-sample evaluation is intentionally NOT performed for the EBA
  # application: the corporate default-rate sample is too short for a
  # meaningful rolling pseudo-out-of-sample exercise. The satellite is
  # estimated on the full sample, and model-selection uncertainty is handled
  # by Akaike model averaging (see output/model_averaging/eba/).
  inf <- influence.measures(best_model)
  infl_idx <- which(apply(inf$is.inf, 1, any))
  
  write.csv(
    data.frame(obs = infl_idx),
    file.path(out_dir, "influential_points.csv"),
    row.names = FALSE
  )
  
  fit_rlm <- MASS::rlm(formula(best_model), data = model_df, psi = psi.huber)
  
  write.table(
    capture.output(summary(fit_rlm)),
    file.path(out_dir, "robust_rlm_summary.txt"),
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
  
  invisible(data.table(
    model_tag = model_tag,
    rho = rho_value,
    p_ttc = p_ttc,
    mean_Z_raw = mean_Z_raw,
    var_Z_raw = var_Z_raw,
    AIC = AIC(best_model),
    BIC = BIC(best_model),
    RMSE_insample = sqrt(mean((model_df$Y - fitted_values)^2)),
    rolling_RMSE = NA_real_,  # OOS evaluation disabled for EBA (short sample)
    rolling_MAE = NA_real_,   # OOS evaluation disabled for EBA (short sample)
    n_obs = nrow(model_df),
    variables = paste(setdiff(names(coef(best_model)), "(Intercept)"), collapse = " + ")
  ))
}

# =============================================================================
# Akaike model averaging diagnostic for the Z-factor satellite
# -----------------------------------------------------------------------------
# Reduces model-selection uncertainty by averaging the best valid OLS satellite
# specifications with Akaike weights. This is a Z-stage diagnostic written to a
# dedicated, clearly-labelled folder; the GIRF/PD pipeline's averaged satellite
# is produced separately by the shared engine (see output/model_averaging/eba/).
# =============================================================================

write_zfactor_model_averaging <- function(valid_models_table, model_df,
                                          model_tag, out_dir) {
  weight_scheme <- if (exists("BMA_WEIGHT")) BMA_WEIGHT else "aic"
  rule <- if (exists("BMA_RULE")) BMA_RULE else "cum_weight"
  cum_threshold <- if (exists("BMA_CUM_THRESHOLD")) BMA_CUM_THRESHOLD else 0.95
  top_n <- if (exists("BMA_TOP_N")) BMA_TOP_N else 10L
  delta_max <- if (exists("BMA_DELTA_MAX")) BMA_DELTA_MAX else 2

  vt <- as.data.table(copy(valid_models_table))
  ic_col <- if (weight_scheme == "bic" && "bic" %in% names(vt)) "bic" else "aic"
  setorderv(vt, ic_col)

  d_all <- vt[[ic_col]] - min(vt[[ic_col]])
  w_all <- if (weight_scheme == "equal") {
    rep(1 / nrow(vt), nrow(vt))
  } else {
    ww <- exp(-0.5 * d_all)
    ww / sum(ww)
  }

  keep <- switch(
    rule,
    cum_weight = seq_len(max(1L, which(cumsum(w_all) >= cum_threshold)[1])),
    top_n      = seq_len(min(as.integer(top_n), nrow(vt))),
    delta      = which(d_all <= delta_max),
    seq_len(max(1L, which(cumsum(w_all) >= cum_threshold)[1]))
  )

  vt_top <- vt[keep]
  d <- vt_top[[ic_col]] - min(vt_top[[ic_col]])
  w <- if (weight_scheme == "equal") {
    rep(1 / nrow(vt_top), nrow(vt_top))
  } else {
    ww <- exp(-0.5 * d)
    ww / sum(ww)
  }
  vt_top[, weight := w]

  # Refit each retained model and accumulate weighted-average coefficients.
  avg <- list()
  incl <- list()
  for (i in seq_len(nrow(vt_top))) {
    vars <- strsplit(vt_top$model[i], " \\+ ")[[1]]
    fit <- lm(as.formula(paste("Y ~", paste(vars, collapse = " + "))), data = model_df)
    cf <- coef(fit)
    for (nm in names(cf)) {
      avg[[nm]] <- (if (is.null(avg[[nm]])) 0 else avg[[nm]]) + w[i] * unname(cf[nm])
      if (nm != "(Intercept)") {
        incl[[nm]] <- (if (is.null(incl[[nm]])) 0 else incl[[nm]]) + w[i]
      }
    }
  }

  averaged_dt <- data.table(
    term = names(avg),
    averaged_coefficient = as.numeric(unlist(avg)),
    inclusion_prob = vapply(names(avg), function(nm) {
      if (nm == "(Intercept)") 1 else as.numeric(incl[[nm]])
    }, numeric(1))
  )
  setorder(averaged_dt, -inclusion_prob)

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  stem <- paste0("eba_zfactor_", model_tag, "_")

  fwrite(vt_top, file.path(out_dir, paste0(stem, "top_models.csv")))
  fwrite(averaged_dt, file.path(out_dir, paste0(stem, "averaged_coefficients.csv")))

  fwrite(
    data.table(
      application = "eba",
      stage = "zfactor_satellite",
      model_tag = model_tag,
      weight_scheme = weight_scheme,
      selection_rule = rule,
      cum_threshold = cum_threshold,
      n_models_averaged = nrow(vt_top),
      n_valid_models = nrow(vt)
    ),
    file.path(out_dir, paste0(stem, "summary.csv"))
  )

  invisible(averaged_dt)
}

# =============================================================================
# Full pipeline for one rho
# =============================================================================

run_pipeline_for_rho <- function(spec, run_design_diag = FALSE) {
  
  model_tag <- spec$model_tag
  rho_value <- spec$rho_value
  rho_source <- spec$rho_source
  
  cat("\n\n################################################################\n")
  cat("PIPELINE FOR ", model_tag, " | rho = ", sprintf("%.8f", rho_value), "\n", sep = "")
  cat("################################################################\n")
  
  Z_out <- f_Z_estimation(
    data_risk_usa$Corpo_DR_WA,
    rho_fixed = rho_value
  )
  
  Z_vec <- Z_out$Z
  
  if (all(is.na(Z_vec))) {
    stop("Z-factor computation failed for ", model_tag)
  }
  
  Y_df <- data.frame(
    Date = floor_date(dates_risk, unit = "month"),
    Y = Z_vec
  )
  
  model_df <- inner_join(Y_df, data_lags, by = "Date") %>%
    na.omit()
  
  if (nrow(model_df) < 20) {
    stop("Not enough valid observations (<20) for ", model_tag)
  }
  
  rho_out_dir <- file.path(z_by_rho_dir, model_tag)
  dir.create(rho_out_dir, showWarnings = FALSE, recursive = TRUE)
  
  model_out_dir <- file.path(models_by_rho_dir, model_tag)
  dir.create(model_out_dir, showWarnings = FALSE, recursive = TRUE)
  
  z_dt <- data.table(
    Date = dates_risk,
    country = "United States",
    default_rate = data_risk_usa$Corpo_DR_WA,
    Z = Z_vec,
    rho = rho_value,
    rho_source = rho_source,
    model_tag = model_tag,
    p_ttc = p_ttc_usa,
    mean_Z_raw = Z_out$mean_Z_raw,
    var_Z_raw = Z_out$var_Z_raw
  )
  
  fwrite(z_dt, file.path(rho_out_dir, "z_factor.csv"))
  
  if (model_tag == "rho_estimated") {
    fwrite(z_dt, file.path(output_dir, "z_factor_usa.csv"))
  }
  
  if (run_design_diag) {
    run_collinearity_diagnostics(model_df, output_dir)
  }
  
  selected <- select_best_ols_model(
    model_df = model_df,
    model_tag = model_tag,
    rho_value = rho_value,
    rho_source = rho_source
  )
  
  best_model <- selected$best_model
  best_model_info <- selected$best_model_info
  
  centerX <- make_centerX_model(model_df, best_model)
  best_model_centerX <- centerX$model
  
  saveRDS(best_model, file = file.path(model_out_dir, "best_model.rds"))
  saveRDS(best_model_centerX, file = file.path(model_out_dir, "best_model_centerX.rds"))
  
  saveRDS(
    list(
      best_model = best_model,
      best_model_centerX = best_model_centerX,
      best_model_info = best_model_info,
      model_tag = model_tag,
      rho_value = rho_value,
      rho_source = rho_source,
      p_ttc = p_ttc_usa,
      mean_Z_raw = Z_out$mean_Z_raw,
      var_Z_raw = Z_out$var_Z_raw,
      centerX_means = centerX$x_means,
      model_df = model_df,
      model_df_centerX = centerX$model_df
    ),
    file = file.path(model_out_dir, "model_bundle.rds")
  )
  
  fwrite(
    selected$valid_models_table,
    file.path(selection_by_rho_dir, paste0(model_tag, "_valid_models.csv"))
  )
  
  if (model_tag == "rho_estimated") {
    saveRDS(best_model, file = file.path(models_dir, "best_model.rds"))
    saveRDS(best_model_centerX, file = file.path(models_dir, "best_model_centerX.rds"))
    cat("Baseline models also saved to 'output/models/'.\n")
  }
  
  post_stats <- run_post_ols_outputs(
    model_df = model_df,
    best_model = best_model,
    best_model_centerX = best_model_centerX,
    centerX_means = centerX$x_means,
    model_info = best_model_info,
    out_dir = rho_out_dir,
    model_tag = model_tag,
    rho_value = rho_value,
    p_ttc = p_ttc_usa,
    mean_Z_raw = Z_out$mean_Z_raw,
    var_Z_raw = Z_out$var_Z_raw
  )

  # Akaike model-averaging diagnostic over the valid satellite specifications.
  model_avg_dir <- if (exists("DIR_MODEL_AVERAGING_EBA")) {
    DIR_MODEL_AVERAGING_EBA
  } else {
    file.path("output", "model_averaging", "eba")
  }

  tryCatch(
    write_zfactor_model_averaging(
      valid_models_table = selected$valid_models_table,
      model_df = model_df,
      model_tag = model_tag,
      out_dir = model_avg_dir
    ),
    error = function(e) {
      cat("Model-averaging diagnostic skipped for ", model_tag, ": ",
          conditionMessage(e), "\n", sep = "")
    }
  )

  run_reg <- if (model_tag == "rho_estimated") {
    RUN_REGULARIZATION_FOR_BASELINE
  } else {
    RUN_REGULARIZATION_FOR_FIXED_RHO
  }
  
  if (run_reg) {
    run_regularization(
      model_df = model_df,
      output_dir = rho_out_dir,
      best_model = best_model,
      model_tag = model_tag,
      rho_value = rho_value,
      p_ttc = p_ttc_usa
    )
    
    if (model_tag == "rho_estimated") {
      reg_src <- file.path(rho_out_dir, "regularization")
      reg_dst <- file.path(output_dir, "regularization")
      
      if (dir.exists(reg_src)) {
        if (dir.exists(reg_dst)) unlink(reg_dst, recursive = TRUE)
        dir.create(dirname(reg_dst), showWarnings = FALSE, recursive = TRUE)
        file.copy(reg_src, dirname(reg_dst), recursive = TRUE)
      }
    }
  }
  
  out <- list(
    model_tag = model_tag,
    rho_value = rho_value,
    rho_source = rho_source,
    best_model = best_model,
    best_model_centerX = best_model_centerX,
    best_model_info = best_model_info,
    centerX_means = centerX$x_means,
    n_valid_models = selected$n_valid_models,
    total_models = selected$total_models,
    post_stats = post_stats
  )
  
  rm(selected, model_df, centerX)
  gc(verbose = FALSE)
  
  out
}

# =============================================================================
# Run all rho-specific pipelines
# =============================================================================

all_results <- list()

for (i in seq_len(nrow(rho_specs))) {
  spec_i <- rho_specs[i]
  
  all_results[[spec_i$model_tag]] <- run_pipeline_for_rho(
    spec = spec_i,
    run_design_diag = (spec_i$model_tag == "rho_estimated")
  )
}

# =============================================================================
# Compact multi-rho exports
# =============================================================================

model_comparison <- rbindlist(
  lapply(all_results, function(x) x$post_stats),
  use.names = TRUE,
  fill = TRUE
)

fwrite(
  model_comparison,
  file.path(output_dir, "model_comparison_by_rho.csv")
)

saveRDS(
  all_results,
  file.path(output_dir, "zfactor_models_by_rho_results.rds")
)

coef_by_rho <- rbindlist(
  lapply(all_results, function(x) {
    cf <- coef(x$best_model)
    data.table(
      model_tag = x$model_tag,
      rho = x$rho_value,
      rho_source = x$rho_source,
      model_type = "raw",
      variable = names(cf),
      coefficient = as.numeric(cf)
    )
  }),
  use.names = TRUE,
  fill = TRUE
)

coef_centerX_by_rho <- rbindlist(
  lapply(all_results, function(x) {
    cf <- coef(x$best_model_centerX)
    data.table(
      model_tag = x$model_tag,
      rho = x$rho_value,
      rho_source = x$rho_source,
      model_type = "centerX_display",
      variable = names(cf),
      coefficient = as.numeric(cf)
    )
  }),
  use.names = TRUE,
  fill = TRUE
)

fwrite(
  coef_by_rho,
  file.path(output_dir, "coefficients_by_rho.csv")
)

fwrite(
  coef_centerX_by_rho,
  file.path(output_dir, "coefficients_centerX_by_rho.csv")
)

fwrite(
  rbindlist(list(coef_by_rho, coef_centerX_by_rho), use.names = TRUE, fill = TRUE),
  file.path(output_dir, "coefficients_raw_and_centerX_by_rho.csv")
)

n_valid_by_rho <- rbindlist(
  lapply(all_results, function(x) {
    data.table(
      model_tag = x$model_tag,
      rho = x$rho_value,
      rho_source = x$rho_source,
      total_models = x$total_models,
      n_valid_models = x$n_valid_models
    )
  }),
  use.names = TRUE,
  fill = TRUE
)

fwrite(
  n_valid_by_rho,
  file.path(output_dir, "n_valid_models_by_rho.csv")
)

merton_params_enriched <- merge(
  rho_specs,
  model_comparison[, .(
    model_tag,
    mean_Z_raw,
    var_Z_raw,
    AIC,
    BIC,
    RMSE_insample,
    rolling_RMSE,
    rolling_MAE,
    variables
  )],
  by = "model_tag",
  all.x = TRUE
)

merton_params_enriched <- merge(
  merton_params_enriched,
  n_valid_by_rho[, .(
    model_tag,
    total_models,
    n_valid_models
  )],
  by = "model_tag",
  all.x = TRUE
)

fwrite(
  merton_params_enriched,
  file.path(output_dir, "merton_params_by_rho_usa.csv")
)

baseline_stats <- merton_params_enriched[model_tag == "rho_estimated"]

fwrite(
  baseline_stats,
  file.path(output_dir, "merton_params_usa.csv")
)

cat("\n================================================================\n")
cat("Analysis completed.\n")
cat("Models by rho saved in: ", models_by_rho_dir, "\n", sep = "")
cat("Z-factor outputs by rho saved in: ", z_by_rho_dir, "\n", sep = "")
cat("Model-selection tables saved in: ", selection_by_rho_dir, "\n", sep = "")
cat("Model comparison: output/zfactor/model_comparison_by_rho.csv\n")
cat("Coefficients by rho: output/zfactor/coefficients_by_rho.csv\n")
cat("Centered-X coefficients by rho: output/zfactor/coefficients_centerX_by_rho.csv\n")
cat("Merton parameters by rho: output/zfactor/merton_params_by_rho_usa.csv\n")
cat("Backward-compatible baseline model: output/models/best_model.rds\n")
cat("Baseline centered-X display model: output/models/best_model_centerX.rds\n")
cat("================================================================\n")