# =============================================================================
#  code/shared/_satellite_regularization.R — Ridge/Lasso robustness for the
#  satellite (Z) regression.
#
#  Portable version of the regularization step used in the EBA z-factor script,
#  so that both applications (EBA and DRALACBN), and every VAR information set,
#  obtain the same shrinkage-based robustness check on the satellite coefficients.
#
#  Ridge (alpha = 0) shrinks all coefficients and is stable under the strong
#  collinearity of lagged macro regressors; it probes the robustness of the
#  coefficient magnitudes. Lasso (alpha = 1) sets coefficients to zero and probes
#  the robustness of the selected set. An optional elastic net (0 < alpha < 1) is
#  available as a compromise.
#
#  Requires: glmnet (guarded by requireNamespace), data.table.
# =============================================================================

# Cross-validation folds that respect the time ordering (contiguous blocks).
.blocked_foldid <- function(n, K = 5) {
  K <- max(3, min(K, n))
  rep(rep(seq_len(K), each = ceiling(n / K)), length.out = n)
}

.drop_zero_variance <- function(X, tol = 1e-12) {
  sds <- apply(X, 2, stats::sd, na.rm = TRUE)
  keep <- is.finite(sds) & sds > tol
  list(X = X[, keep, drop = FALSE], kept = colnames(X)[keep],
       dropped = colnames(X)[!keep])
}

.reg_fit_stats <- function(y, yhat, df, name) {
  n <- length(y)
  rss <- sum((y - yhat)^2)
  mse <- rss / n
  data.frame(
    model = name, n = n, df = df,
    RMSE = sqrt(mse), MSE = mse, RSS = rss,
    AIC = n * log(mse) + 2 * df,
    BIC = n * log(mse) + log(n) * df,
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# Run ridge and/or lasso (and optional elastic net) on the satellite design and
# write coefficients, CV curves, fit metrics and a summary to out_dir.
# -----------------------------------------------------------------------------
run_satellite_regularization <- function(model_df, out_dir,
                                         label = "",
                                         run_ridge = TRUE,
                                         run_lasso = TRUE,
                                         run_elastic_net = FALSE,
                                         elastic_alpha = 0.5,
                                         use_lambda = "lambda.1se",
                                         k_folds = 5,
                                         seed = 12345) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    cat("    Ridge/Lasso robustness skipped: package glmnet unavailable.\n")
    return(invisible(FALSE))
  }

  stopifnot(all(c("Date", "Y") %in% names(model_df)))
  model_df <- as.data.frame(model_df)

  y_vec <- model_df$Y
  X_all <- as.matrix(model_df[, setdiff(names(model_df), c("Date", "Y")),
                              drop = FALSE])
  nzv <- .drop_zero_variance(X_all)
  X_mat <- nzv$X

  if (ncol(X_mat) < 2L || nrow(X_mat) < 5L) {
    cat("    Ridge/Lasso robustness skipped for ", label,
        ": design too small.\n", sep = "")
    return(invisible(FALSE))
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  foldid <- .blocked_foldid(nrow(X_mat), K = k_folds)
  n_folds <- length(unique(foldid))

  write_nonzero <- function(fit_cv, path_csv) {
    cf <- as.matrix(coef(fit_cv, s = use_lambda))
    out <- data.frame(variable = rownames(cf), coef = as.numeric(cf[, 1]),
                      row.names = NULL)
    out$nonzero <- out$coef != 0
    out <- out[order(-abs(out$coef)), ]
    write.csv(out, path_csv, row.names = FALSE)
  }

  take_lambda <- function(cvfit) {
    if (use_lambda == "lambda.min") cvfit$lambda.min else cvfit$lambda.1se
  }

  metrics <- list()

  fit_one <- function(alpha, tag, title) {
    set.seed(seed)
    cvfit <- glmnet::cv.glmnet(
      x = X_mat, y = y_vec, family = "gaussian", alpha = alpha,
      standardize = TRUE, foldid = foldid, type.measure = "mse",
      nfolds = n_folds
    )
    saveRDS(cvfit, file = file.path(out_dir, paste0(tag, "_cvglmnet.rds")))

    grDevices::png(file.path(out_dir, paste0(tag, "_cv_curve.png")),
                   width = 1000, height = 600)
    plot(cvfit)
    graphics::title(title)
    grDevices::dev.off()

    write_nonzero(cvfit, file.path(out_dir,
                                   paste0(tag, "_coefs_", use_lambda, ".csv")))

    yhat <- as.numeric(predict(cvfit, newx = X_mat, s = use_lambda))
    lam <- take_lambda(cvfit)
    idx <- which.min(abs(cvfit$lambda - lam))
    df_eff <- cvfit$glmnet.fit$df[idx]
    metrics[[tag]] <<- .reg_fit_stats(y_vec, yhat, df = df_eff,
                                      name = paste0(toupper(tag), "_", use_lambda))
  }

  if (run_ridge) fit_one(0, "ridge", "Ridge CV (MSE)")
  if (run_lasso) fit_one(1, "lasso", "Lasso CV (MSE)")
  if (run_elastic_net) {
    fit_one(elastic_alpha, "elastic_net",
            sprintf("Elastic net CV (MSE), alpha=%.2f", elastic_alpha))
  }

  if (length(metrics)) {
    met <- do.call(rbind, metrics)
    write.csv(met, file.path(out_dir, "regularization_metrics.csv"),
              row.names = FALSE)
  }

  con <- file(file.path(out_dir, "regularization_summary.txt"), open = "wt")
  on.exit(close(con), add = TRUE)
  wln <- function(...) cat(paste0(..., collapse = ""), "\n", file = con,
                           append = TRUE)
  wln("Satellite regularization robustness")
  wln("Label: ", label)
  wln("Generated on: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
  wln("n_obs = ", nrow(X_mat), " | n_predictors = ", ncol(X_mat))
  if (length(nzv$dropped)) {
    wln("Zero-variance columns dropped: ", paste(nzv$dropped, collapse = ", "))
  }
  wln("Lambda rule = ", use_lambda, " | CV folds (blocked) = ", n_folds)
  wln("Ridge (alpha=0): stable shrinkage under collinearity.")
  wln("Lasso (alpha=1): sparse selection of surviving regressors.")

  invisible(TRUE)
}

# -----------------------------------------------------------------------------
# Robust standard errors for the selected (single best) satellite: HAC
# (Newey-West, prewhitened) and HC3. Frequentist diagnostic written alongside
# the model-averaged Bayesian satellite. Pure side-effect (writes files).
# -----------------------------------------------------------------------------
write_satellite_robust_se <- function(sat_model, out_dir, nb_lags,
                                      app_label, info_set_label,
                                      app_name, info_set_name) {
  if (!(requireNamespace("lmtest", quietly = TRUE) &&
        requireNamespace("sandwich", quietly = TRUE))) {
    cat("    HAC/HC3 robustness skipped: packages lmtest/sandwich unavailable.\n")
    return(invisible(FALSE))
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  tryCatch({
    hac_ct <- lmtest::coeftest(
      sat_model,
      vcov. = sandwich::NeweyWest(sat_model, lag = nb_lags,
                                  prewhite = TRUE, adjust = TRUE)
    )
    hc3_ct <- lmtest::coeftest(
      sat_model, vcov. = sandwich::vcovHC(sat_model, type = "HC3")
    )

    writeLines(capture.output(hac_ct),
               file.path(out_dir, "satellite_hac_newey_west.txt"))
    writeLines(capture.output(hc3_ct),
               file.path(out_dir, "satellite_hc3_robust_se.txt"))
    fwrite(as.data.table(as.data.frame(unclass(hac_ct)), keep.rownames = "term"),
           file.path(out_dir, "satellite_hac_newey_west.csv"))
    fwrite(as.data.table(as.data.frame(unclass(hc3_ct)), keep.rownames = "term"),
           file.path(out_dir, "satellite_hc3_robust_se.csv"))

    write_coeftest_table_tex(
      hac_ct,
      file.path(out_dir, "satellite_hac_newey_west.tex"),
      caption = paste0(app_label,
        ": selected satellite, HAC (Newey-West) standard errors (",
        info_set_label, " information set)"),
      label = paste0("tab:", app_name, "_", info_set_name, "_satellite_hac"),
      notes = paste0(
        "Notes: OLS point estimates of the single best (AIC-selected) ",
        "satellite with heteroskedasticity- and autocorrelation-consistent ",
        "(Newey-West, ", nb_lags, " lags, prewhitened) standard errors. ",
        "Significance: $^{*}p<0.05$, $^{**}p<0.01$, $^{***}p<0.001$.")
    )
  }, error = function(e) {
    cat("    HAC/HC3 robustness skipped for ", info_set_name, ": ",
        conditionMessage(e), "\n", sep = "")
  })

  invisible(TRUE)
}
