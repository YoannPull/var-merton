# =============================================================================
#  code/extra/coralacbn/03_coralacbn_direct_channel.R
#
#  CORALACBN application — relaxing the satellite exogeneity assumption
#  (eta_t independent of VAR innovations) via a control-function term.
#
#  Model:
#      Z_t = beta0 + beta' Y_t^(s) + lambda * e_t^GPR + eps_t,
#  where e_t^GPR = u_{g,t} / sqrt(sigma_gg) is the STANDARDIZED reduced-form
#  GPR innovation from the macro-financial VAR. The baseline specification of
#  the paper is the restriction lambda = 0 (no direct geopolitical channel to
#  credit beyond the macro-financial block).
#
#  Identification note: because Psi_0 = I_n, lambda is identified only for
#  innovations of variables EXCLUDED from the contemporaneous satellite
#  selection. The GPR index never enters the satellite, so lambda_GPR is
#  identified (up to multicollinearity with correlated innovations).
#
#  Estimation is hierarchical: for each retained BVAR posterior draw m,
#  (i) the innovation series e^(m) is recomputed from (B_m, Sigma_m),
#  (ii) the augmented satellite is re-estimated by conjugate Bayesian linear
#  regression (Jeffreys prior) on the overlap sample, and one posterior draw
#  (beta0, beta, lambda, sigma_eps^2)_m is taken. Generated-regressor
#  uncertainty is therefore propagated automatically.
#
#  Closed-form GIRF corrections (exact under joint Gaussianity):
#    psi_Z(0)            += lambda * delta_sd        (delta_sd = shock in s.d. units)
#    s2(h), all h        += lambda^2 + 2 * lambda * c0,   c0 = beta' S_0 Sigma e_g / sqrt(sigma_gg)
#    s2_delta(h), h >= 1 += lambda^2 + 2 * lambda * c0    (h = 0: e_t is fixed by conditioning)
#  All Merton-Vasicek formulas (mean PD response) are then applied unchanged.
#
#  This script is STANDALONE. It does not modify any existing pipeline file
#  and is NOT part of run_all.R. Run from the repository root:
#      source("code/extra/coralacbn/03_coralacbn_direct_channel.R")
#
#  Outputs: output/applications/coralacbn/direct_channel/
# =============================================================================

source("code/00_setup.R")
source("code/shared/_helpers_var_merton.R")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# -----------------------------------------------------------------------------
# 0. Parameters and output tree
# -----------------------------------------------------------------------------

DC_INFO_SET     <- if (exists("DELINQUENCY_BASELINE_INFO_SET")) DELINQUENCY_BASELINE_INFO_SET else "real_side"
DC_QUICK_TEST   <- FALSE          # TRUE: use only DC_TEST_DRAWS posterior draws
DC_TEST_DRAWS   <- 500L

p_lags       <- if (exists("P_LAGS")) P_LAGS else 2
nrep         <- if (exists("NREP")) NREP else 30000
M_target     <- if (exists("M_TARGET")) M_TARGET else 10000L
H            <- if (exists("HORIZON")) HORIZON else 12
seed         <- if (exists("SEED_ROBUSTNESS")) SEED_ROBUSTNESS else 123
seed_sat     <- if (exists("SEED_ZFACTOR")) SEED_ZFACTOR else 12345
impulse_name <- if (exists("IMPULSE_NAME")) IMPULSE_NAME else "log_GPRD"
NB_LAGS_SAT  <- if (exists("ALT_SPEC_NB_LAGS_SAT")) ALT_SPEC_NB_LAGS_SAT else 4

dc_out_dir <- file.path("output", "applications", "coralacbn", "direct_channel")
for (sub in c("", "tables", "figures", "irf")) {
  dir.create(file.path(dc_out_dir, sub), recursive = TRUE, showWarnings = FALSE)
}

cat("\n========== CORALACBN direct-channel (control-function) exercise ==========\n")
cat("Information set: ", DC_INFO_SET, " | impulse: ", impulse_name, "\n", sep = "")

# -----------------------------------------------------------------------------
# 1. Reconstruct the CORALACBN systematic factor (identical to wrapper 02)
# -----------------------------------------------------------------------------

delinquency_path <- file.path(
  if (exists("DIR_RAW")) DIR_RAW else "data/raw", "default", "CORALACBN.csv"
)
if (!file.exists(delinquency_path)) {
  stop("Missing charge-off series: ", delinquency_path, call. = FALSE)
}

delinq_raw <- fread(delinquency_path)
date_col <- intersect(c("DATE", "Date", "date", "observation_date"), names(delinq_raw))[1]
val_col  <- intersect(c("CORALACBN", "value", "VALUE", "Value"), names(delinq_raw))[1]
if (is.na(date_col)) date_col <- names(delinq_raw)[1]
if (is.na(val_col))  val_col  <- setdiff(names(delinq_raw), date_col)[1]

delinq <- data.table(
  Date = as.Date(delinq_raw[[date_col]]),
  value = as.numeric(delinq_raw[[val_col]])
)
delinq <- delinq[is.finite(value)]
setorder(delinq, Date)
if (exists("DATA_END_DATE")) delinq <- delinq[Date <= DATA_END_DATE]

z_out   <- f_Z_estimation(delinq$value / 100)
APP_RHO   <- z_out$rho
APP_P_TTC <- z_out$p_ttc

Z_DT <- data.table(
  Date = as.Date(format(delinq$Date + 62, "%Y-%m-01")),
  Z = z_out$Z
)

cat("Merton parameters: p_ttc = ", sprintf("%.6f", APP_P_TTC),
    ", rho = ", sprintf("%.6f", APP_RHO), "\n", sep = "")

# -----------------------------------------------------------------------------
# 2. Load data and the existing BVAR kernel (reused, never re-estimated here)
# -----------------------------------------------------------------------------

DT_raw <- safe_fread(if (exists("PATH_DATA_VAR")) PATH_DATA_VAR else "data/processed/data_var_for_model.csv")
DT_raw[, Date := as.Date(Date)]

kern <- get_or_estimate_kernel(
  spec_name = DC_INFO_SET,
  vars_extra = VAR_SPECS[[DC_INFO_SET]],
  DT_raw = DT_raw,
  p_lags = p_lags, nrep = nrep, H = H, seed = seed, M_target = M_target,
  kernel_dir = if (exists("DIR_VAR_KERNELS")) DIR_VAR_KERNELS else file.path("output", "shared", "var_kernels"),
  impulse_name = impulse_name
)

variables  <- kern$variables
impulse_ix <- match(impulse_name, variables)
if (is.na(impulse_ix)) stop("Impulse variable not found in VAR kernel.", call. = FALSE)

M_full <- dim(kern$Psi_draws)[4]
draw_ids <- seq_len(M_full)
if (isTRUE(DC_QUICK_TEST) && M_full > DC_TEST_DRAWS) {
  set.seed(seed)
  draw_ids <- sort(sample.int(M_full, DC_TEST_DRAWS))
}
M <- length(draw_ids)

Psi_draws  <- kern$Psi_draws[, , , draw_ids, drop = FALSE]
B_kept     <- kern$B_kept[, , draw_ids, drop = FALSE]
Sigma_kept <- kern$Sigma_kept[, , draw_ids, drop = FALSE]
kern_sub <- kern
kern_sub$Psi_draws  <- Psi_draws
kern_sub$B_kept     <- B_kept
kern_sub$Sigma_kept <- Sigma_kept

cat("Posterior draws used: ", M, " (of ", M_full, ")\n", sep = "")

# -----------------------------------------------------------------------------
# 3. Effective dates of the VAR sample and per-draw GPR innovations
# -----------------------------------------------------------------------------
# Replicates the row selection of estimate_bvar_kernel(): complete cases of the
# VAR information set, in order. Effective residual rows are (p+1):n.

cc_mask <- stats::complete.cases(as.data.frame(DT_raw[, ..variables]))
dates_kept <- as.Date(DT_raw$Date)[cc_mask]
if (length(dates_kept) != nrow(kern$DT)) {
  stop("Date alignment failed: complete-cases rows (", length(dates_kept),
       ") do not match kernel rows (", nrow(kern$DT), ").", call. = FALSE)
}

Y_mat <- as.matrix(kern$DT[, ..variables])
n_var <- nrow(Y_mat)
Y_t <- Y_mat[(p_lags + 1):n_var, , drop = FALSE]
X_t <- do.call(
  cbind,
  lapply(1:p_lags, function(l) Y_mat[(p_lags + 1 - l):(n_var - l), , drop = FALSE])
)
X_t <- cbind(1, X_t)

dates_eff <- dates_kept[(p_lags + 1):n_var]
dates_eff_shifted <- as.Date(format(dates_eff + 62, "%Y-%m-01"))

# Standardized GPR innovation per posterior draw: e^(m) = u_g^(m) / sqrt(sigma_gg^(m)).
E_innov <- matrix(NA_real_, nrow = nrow(Y_t), ncol = M)
for (m in seq_len(M)) {
  U_m <- Y_t - X_t %*% B_kept[, , m]
  sgg <- Sigma_kept[impulse_ix, impulse_ix, m]
  if (!is.finite(sgg) || sgg <= 0) next
  E_innov[, m] <- U_m[, impulse_ix] / sqrt(sgg)
}

# -----------------------------------------------------------------------------
# 4. Satellite design: best-BIC specification of the main CORALACBN run
# -----------------------------------------------------------------------------

best_sat_path <- file.path("output", "applications", "coralacbn",
                           "main_specification", "satellite", "best_satellite.rds")
if (!file.exists(best_sat_path)) {
  stop("Missing ", best_sat_path,
       "\nRun the main CORALACBN application first (run_all.R).", call. = FALSE)
}
best_sat <- readRDS(best_sat_path)
sat_terms <- setdiff(names(coef(best_sat)), "(Intercept)")

parsed   <- parse_satellite_terms(best_sat, variables)
terms_df <- parsed$terms        # term / base / lag / beta(OLS) / col

vars_sat <- setdiff(variables, impulse_name)
sat_df <- prepare_lagged_satellite_data(
  z_input = Z_DT, DT_raw = DT_raw, vars_sat = vars_sat, nb_lags = NB_LAGS_SAT
)
sat_df <- as.data.table(sat_df)

# Align satellite rows with VAR innovation dates (contemporaneous match).
row_map <- match(sat_df$Date, dates_eff_shifted)
keep <- !is.na(row_map)
if (any(!keep)) {
  cat("Dropping ", sum(!keep),
      " satellite rows without a matching VAR innovation date.\n", sep = "")
}
sat_df  <- sat_df[keep]
row_map <- row_map[keep]
n_sat_obs <- nrow(sat_df)

cat("Overlap sample for the augmented satellite: ", n_sat_obs, " quarters (",
    format(min(sat_df$Date)), " to ", format(max(sat_df$Date)), ")\n", sep = "")
cat("Satellite design (best-BIC spec): ",
    paste(sat_terms, collapse = " + "), " + gpr_innov\n", sep = "")

y_vec  <- as.numeric(sat_df$Y)
X_base <- cbind(`(Intercept)` = 1, as.matrix(sat_df[, ..sat_terms]))
k_base <- ncol(X_base)

# -----------------------------------------------------------------------------
# 5. Hierarchical conjugate estimation, one satellite draw per BVAR draw
# -----------------------------------------------------------------------------

draw_conjugate_once <- function(X, y) {
  n <- nrow(X); p <- ncol(X); df <- n - p
  XtX_inv <- tryCatch(chol2inv(chol(crossprod(X))), error = function(e) NULL)
  if (is.null(XtX_inv)) return(NULL)
  beta_hat <- as.numeric(XtX_inv %*% crossprod(X, y))
  rss <- sum((y - X %*% beta_hat)^2)
  s2_hat <- rss / df
  sigma2 <- df * s2_hat / stats::rchisq(1L, df = df)
  R <- chol(XtX_inv)
  beta <- beta_hat + sqrt(sigma2) * as.numeric(crossprod(R, stats::rnorm(p)))
  list(beta = beta, beta_hat = beta_hat, sigma2 = sigma2, s2_hat = s2_hat, rss = rss)
}

set.seed(seed_sat)

# --- Augmented satellite (with the standardized GPR innovation) --------------
beta0_aug  <- numeric(M)
betaT_aug  <- matrix(NA_real_, nrow = M, ncol = length(sat_terms),
                     dimnames = list(NULL, sat_terms))
lambda_aug <- numeric(M)
sigma2_aug <- numeric(M)
lambda_hat <- numeric(M)        # per-draw OLS point estimate (diagnostic)
r2_aug     <- numeric(M)

# --- Restricted satellite (lambda = 0), same sample, same machinery ----------
beta0_res  <- numeric(M)
betaT_res  <- matrix(NA_real_, nrow = M, ncol = length(sat_terms),
                     dimnames = list(NULL, sat_terms))
sigma2_res <- numeric(M)

tss <- sum((y_vec - mean(y_vec))^2)

for (m in seq_len(M)) {
  e_m <- E_innov[row_map, m]
  if (anyNA(e_m)) next

  X_aug <- cbind(X_base, gpr_innov = e_m)
  d_aug <- draw_conjugate_once(X_aug, y_vec)
  d_res <- draw_conjugate_once(X_base, y_vec)
  if (is.null(d_aug) || is.null(d_res)) next

  beta0_aug[m]   <- d_aug$beta[1L]
  betaT_aug[m, ] <- d_aug$beta[2:k_base]
  lambda_aug[m]  <- d_aug$beta[k_base + 1L]
  sigma2_aug[m]  <- d_aug$sigma2
  lambda_hat[m]  <- d_aug$beta_hat[k_base + 1L]
  r2_aug[m]      <- 1 - d_aug$rss / tss

  beta0_res[m]   <- d_res$beta[1L]
  betaT_res[m, ] <- d_res$beta[2:k_base]
  sigma2_res[m]  <- d_res$sigma2
}

ok <- is.finite(lambda_aug) & is.finite(sigma2_aug) & sigma2_aug > 0
if (!all(ok)) {
  cat("Dropping ", sum(!ok), " draws with degenerate designs.\n", sep = "")
  keep_m <- which(ok)
} else {
  keep_m <- seq_len(M)
}

# Subset every per-draw object consistently.
Psi_draws  <- Psi_draws[, , , keep_m, drop = FALSE]
B_kept     <- B_kept[, , keep_m, drop = FALSE]
Sigma_kept <- Sigma_kept[, , keep_m, drop = FALSE]
kern_sub$Psi_draws  <- Psi_draws
kern_sub$B_kept     <- B_kept
kern_sub$Sigma_kept <- Sigma_kept
beta0_aug  <- beta0_aug[keep_m];  betaT_aug  <- betaT_aug[keep_m, , drop = FALSE]
lambda_aug <- lambda_aug[keep_m]; sigma2_aug <- sigma2_aug[keep_m]
lambda_hat <- lambda_hat[keep_m]; r2_aug     <- r2_aug[keep_m]
beta0_res  <- beta0_res[keep_m];  betaT_res  <- betaT_res[keep_m, , drop = FALSE]
sigma2_res <- sigma2_res[keep_m]
M <- length(keep_m)

# bsat-like objects consumed by the existing GIRF machinery.
# n_draws = M ensures IDENTITY pairing in pair_satellite_draw(): satellite
# draw m is the one estimated conditional on BVAR draw m (hierarchical scheme).
bsat_aug <- list(
  terms = terms_df,
  beta0_draws = beta0_aug,
  beta_term_draws = betaT_aug,
  sigma2_draws = sigma2_aug,
  n_draws = M
)
bsat_res <- list(
  terms = terms_df,
  beta0_draws = beta0_res,
  beta_term_draws = betaT_res,
  sigma2_draws = sigma2_res,
  n_draws = M
)

# -----------------------------------------------------------------------------
# 6. Posterior of lambda and exclusion test (lambda = 0)
# -----------------------------------------------------------------------------

qs <- stats::quantile(lambda_aug, c(0.05, 0.16, 0.50, 0.84, 0.95))
lambda_summary <- data.table(
  parameter = "lambda (direct GPR channel, per 1-s.d. innovation)",
  mean = mean(lambda_aug), sd = stats::sd(lambda_aug),
  q05 = qs[[1]], q16 = qs[[2]], median = qs[[3]], q84 = qs[[4]], q95 = qs[[5]],
  prob_negative = mean(lambda_aug < 0),
  prob_positive = mean(lambda_aug > 0),
  n_obs = n_sat_obs,
  n_draws = M,
  mean_r2_augmented = mean(r2_aug)
)
fwrite(lambda_summary, file.path(dc_out_dir, "tables", "lambda_posterior_summary.csv"))

# Frequentist diagnostic at the posterior-median innovation series (HAC).
e_med <- apply(E_innov[row_map, keep_m, drop = FALSE], 1, stats::median)
ols_df <- data.frame(Y = y_vec, X_base[, -1, drop = FALSE], gpr_innov = e_med,
                     check.names = FALSE)
ols_aug <- stats::lm(
  stats::as.formula(paste("Y ~", paste(c(sat_terms, "gpr_innov"), collapse = " + "))),
  data = ols_df
)
hac_row <- tryCatch({
  if (requireNamespace("sandwich", quietly = TRUE) &&
      requireNamespace("lmtest", quietly = TRUE)) {
    ct <- lmtest::coeftest(ols_aug, vcov. = sandwich::NeweyWest(ols_aug, lag = NB_LAGS_SAT))
    data.table(
      estimate = ct["gpr_innov", 1], hac_se = ct["gpr_innov", 2],
      t_value = ct["gpr_innov", 3], p_value = ct["gpr_innov", 4],
      note = "OLS with gpr_innov fixed at its posterior-median series; Newey-West HAC."
    )
  } else {
    sm <- summary(ols_aug)$coefficients
    data.table(
      estimate = sm["gpr_innov", 1], hac_se = sm["gpr_innov", 2],
      t_value = sm["gpr_innov", 3], p_value = sm["gpr_innov", 4],
      note = "sandwich/lmtest unavailable; classical OLS standard errors."
    )
  }
}, error = function(e) data.table(estimate = NA_real_, hac_se = NA_real_,
                                  t_value = NA_real_, p_value = NA_real_,
                                  note = conditionMessage(e)))
fwrite(hac_row, file.path(dc_out_dir, "tables", "lambda_ols_hac_diagnostic.csv"))

# Coefficient comparison: restricted vs augmented (posterior means and sds).
comp <- data.table(
  term = c("(Intercept)", sat_terms, "gpr_innov"),
  mean_restricted = c(mean(beta0_res), colMeans(betaT_res), NA_real_),
  sd_restricted   = c(stats::sd(beta0_res), apply(betaT_res, 2, stats::sd), NA_real_),
  mean_augmented  = c(mean(beta0_aug), colMeans(betaT_aug), mean(lambda_aug)),
  sd_augmented    = c(stats::sd(beta0_aug), apply(betaT_aug, 2, stats::sd),
                      stats::sd(lambda_aug))
)
fwrite(comp, file.path(dc_out_dir, "tables", "satellite_restricted_vs_augmented.csv"))

# -----------------------------------------------------------------------------
# 7. GIRFs: factor moments, psi_Z and PD, with exact closed-form corrections
# -----------------------------------------------------------------------------

# Lag-0 cross-term c0_m = beta_m' S_0 Sigma_m e_g / sqrt(sigma_gg,m).
lag0_rows <- which(terms_df$lag == 0L)
c0_draws <- numeric(M)
if (length(lag0_rows) > 0L) {
  for (m in seq_len(M)) {
    Sigma_m <- 0.5 * (Sigma_kept[, , m] + t(Sigma_kept[, , m]))
    sgg <- Sigma_m[impulse_ix, impulse_ix]
    if (!is.finite(sgg) || sgg <= 0) next
    c0_draws[m] <- sum(betaT_aug[m, lag0_rows] *
                         Sigma_m[terms_df$col[lag0_rows], impulse_ix]) / sqrt(sgg)
  }
}

bands_from_draws <- function(draw_mat, extra = list()) {
  qs <- c(0.05, 0.16, 0.50, 0.84, 0.95)
  qmat <- t(apply(draw_mat, 1, stats::quantile, probs = qs, na.rm = TRUE))
  out <- data.table(
    horizon = 0:(nrow(draw_mat) - 1L),
    lower90 = qmat[, 1], lower68 = qmat[, 2], median = qmat[, 3],
    upper68 = qmat[, 4], upper90 = qmat[, 5]
  )
  for (nm in names(extra)) out[, (nm) := extra[[nm]]]
  out[]
}

# Conditional moments do not depend on the shock size: compute them once.
cat("Computing conditional factor moments (restricted) ...\n")
mom_res <- compute_factor_moment_draws_bayes(kern_sub, bsat_res,
                                             impulse_ix, p_lags, H)
cat("Computing conditional factor moments (augmented) ...\n")
mom_aug <- compute_factor_moment_draws_bayes(kern_sub, bsat_aug,
                                             impulse_ix, p_lags, H)

# Exact closed-form lambda corrections to the conditional variances.
corr <- lambda_aug^2 + 2 * lambda_aug * c0_draws        # length M
s2_corr  <- pmax(sweep(mom_aug$s2_draws, 2, corr, `+`), 0)
s2d_corr <- mom_aug$s2_delta_draws
if (H >= 1) {
  s2d_corr[2:(H + 1L), ] <- pmax(sweep(s2d_corr[2:(H + 1L), , drop = FALSE],
                                       2, corr, `+`), 0)
}

run_one_shock <- function(shock_scale, shock_lab) {
  # Restricted model (lambda = 0): existing machinery, unchanged.
  psiZ_res <- inject_girf_into_Z_bayes(Psi_draws, Sigma_kept, bsat_res,
                                       impulse_ix, shock_scale = shock_scale)
  pd_res <- compute_pd_girf(psiZ_res$draws, mom_res$mu_draws,
                            mom_res$s2_draws, mom_res$s2_delta_draws,
                            p = APP_P_TTC, rho = APP_RHO)

  # Augmented model: beta-channel through the same machinery, plus the
  # direct-channel impact correction at h = 0.
  psiZ_aug <- inject_girf_into_Z_bayes(Psi_draws, Sigma_kept, bsat_aug,
                                       impulse_ix, shock_scale = shock_scale)
  psiZ_aug_draws <- psiZ_aug$draws
  psiZ_aug_draws[1L, ] <- psiZ_aug_draws[1L, ] + lambda_aug * shock_scale

  pd_aug <- compute_pd_girf(psiZ_aug_draws, mom_aug$mu_draws,
                            s2_corr, s2d_corr,
                            p = APP_P_TTC, rho = APP_RHO)

  list(
    shock = shock_lab,
    psiZ_res_bands = bands_from_draws(psiZ_res$draws, list(model = "restricted", shock = shock_lab)),
    psiZ_aug_bands = bands_from_draws(psiZ_aug_draws, list(model = "augmented", shock = shock_lab)),
    pd_res_bands = bands_from_draws(pd_res$draws, list(model = "restricted", shock = shock_lab)),
    pd_aug_bands = bands_from_draws(pd_aug$draws, list(model = "augmented", shock = shock_lab))
  )
}

shock_2001 <- compute_structural_e1_median(
  kernel = kern_sub, dates_raw = dates_kept,
  target_quarter = "2001Q3", impulse_idx = impulse_ix
)

cat("Computing GIRFs (restricted and augmented), 1-s.d. shock ...\n")
res_1sd <- run_one_shock(1, "one_sd")
cat("Computing GIRFs (restricted and augmented), 2001Q3 shock ...\n")
res_2001 <- run_one_shock(shock_2001$shock_scale, "2001Q3")

all_psiZ <- rbindlist(list(res_1sd$psiZ_res_bands, res_1sd$psiZ_aug_bands,
                           res_2001$psiZ_res_bands, res_2001$psiZ_aug_bands))
all_pd <- rbindlist(list(res_1sd$pd_res_bands, res_1sd$pd_aug_bands,
                         res_2001$pd_res_bands, res_2001$pd_aug_bands))

fwrite(all_psiZ, file.path(dc_out_dir, "irf", "psiZ_bands_restricted_vs_augmented.csv"))
fwrite(all_pd,   file.path(dc_out_dir, "irf", "pd_girf_bands_restricted_vs_augmented.csv"))

saveRDS(
  list(
    info_set = DC_INFO_SET, impulse = impulse_name,
    n_obs = n_sat_obs, n_draws = M,
    merton = list(p_ttc = APP_P_TTC, rho = APP_RHO),
    sat_terms = sat_terms,
    lambda_draws = lambda_aug, lambda_hat_ols = lambda_hat,
    c0_draws = c0_draws,
    lambda_summary = lambda_summary, hac_diagnostic = hac_row,
    coefficients_comparison = comp,
    shock_2001 = shock_2001$target,
    psiZ_bands = all_psiZ, pd_bands = all_pd
  ),
  file.path(dc_out_dir, "direct_channel_results.rds")
)

# -----------------------------------------------------------------------------
# 8. Figures
# -----------------------------------------------------------------------------

plot_lambda_density <- function(lam) {
  ggplot(data.frame(lambda = lam), aes(x = lambda)) +
    geom_density(fill = "#AB4A7D", alpha = 0.35, color = "#6F1732", linewidth = 0.9) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey35", linewidth = 0.4) +
    geom_vline(xintercept = stats::median(lam), color = "#6F1732", linewidth = 0.6) +
    labs(x = expression(lambda ~ "(response of Z to a 1-s.d. GPR innovation, direct channel)"),
         y = "Posterior density") +
    theme_robustness()
}

plot_model_overlay <- function(bands, ylab_txt) {
  dt <- as.data.table(copy(bands))
  dt[, model_label := fifelse(model == "augmented",
                              "Augmented (direct channel)",
                              "Restricted (baseline)")]
  ggplot(dt, aes(x = horizon, y = median, color = model_label,
                 fill = model_label, linetype = model_label)) +
    geom_ribbon(aes(ymin = lower68, ymax = upper68), alpha = 0.18, color = NA) +
    geom_hline(yintercept = 0, color = "grey45", linewidth = 0.35, linetype = "dashed") +
    geom_line(linewidth = 1.0) +
    scale_color_manual(values = c("Augmented (direct channel)" = "#6F1732",
                                  "Restricted (baseline)" = "#6B6B6B")) +
    scale_fill_manual(values = c("Augmented (direct channel)" = "#AB4A7D",
                                 "Restricted (baseline)" = "#6B6B6B")) +
    scale_linetype_manual(values = c("Augmented (direct channel)" = "solid",
                                     "Restricted (baseline)" = "longdash")) +
    scale_x_continuous(
      breaks = function(x) seq(max(0, ceiling(x[1])), floor(x[2]), by = 1),
      minor_breaks = NULL
    ) +
    labs(x = "Horizon (quarters)", y = ylab_txt,
         caption = "Shaded bands: 68% posterior intervals.") +
    theme_robustness() +
    theme(plot.caption = element_blank())
}

ggsave_both(plot_lambda_density(lambda_aug),
            file.path(dc_out_dir, "figures"), "lambda_posterior_density",
            width = 6.5, height = 4.2)

ggsave_both(plot_model_overlay(all_pd[shock == "one_sd"], expression(Delta ~ PD)),
            file.path(dc_out_dir, "figures"), "pd_girf_overlay_1sd",
            width = 7, height = 4.4)
ggsave_both(plot_model_overlay(all_pd[shock == "2001Q3"], expression(Delta ~ PD)),
            file.path(dc_out_dir, "figures"), "pd_girf_overlay_2001Q3",
            width = 7, height = 4.4)
ggsave_both(plot_model_overlay(all_psiZ[shock == "one_sd"], expression(psi[Z](h))),
            file.path(dc_out_dir, "figures"), "z_girf_overlay_1sd",
            width = 7, height = 4.4)
ggsave_both(plot_model_overlay(all_psiZ[shock == "2001Q3"], expression(psi[Z](h))),
            file.path(dc_out_dir, "figures"), "z_girf_overlay_2001Q3",
            width = 7, height = 4.4)

# -----------------------------------------------------------------------------
# 9. LaTeX table for the paper
# -----------------------------------------------------------------------------

tex_path <- file.path(dc_out_dir, "tables", "lambda_posterior_summary.tex")
con <- file(tex_path, open = "wt")
cat("\\begin{table}[!htbp]\n\\centering\n", file = con)
cat("\\caption{Direct geopolitical channel in the satellite equation (CORALACBN)}\n", file = con)
cat("\\label{tab:coralacbn_direct_channel}\n", file = con)
cat("\\begin{tabular}{lr}\n\\toprule\n", file = con)
cat(sprintf("Posterior mean of $\\lambda$ & %.4f \\\\\n", lambda_summary$mean), file = con)
cat(sprintf("Posterior s.d. & %.4f \\\\\n", lambda_summary$sd), file = con)
cat(sprintf("68\\%% credible interval & [%.4f, %.4f] \\\\\n",
            lambda_summary$q16, lambda_summary$q84), file = con)
cat(sprintf("90\\%% credible interval & [%.4f, %.4f] \\\\\n",
            lambda_summary$q05, lambda_summary$q95), file = con)
cat(sprintf("$\\Pr(\\lambda < 0 \\mid \\text{data})$ & %.3f \\\\\n",
            lambda_summary$prob_negative), file = con)
cat(sprintf("Observations (overlap sample) & %d \\\\\n", n_sat_obs), file = con)
cat(sprintf("Posterior draws & %d \\\\\n", M), file = con)
cat("\\bottomrule\n\\end{tabular}\n", file = con)
cat("\\begin{tablenotes}\n\\small\n", file = con)
cat("\\item Notes: $\\lambda$ is the coefficient on the standardized reduced-form GPR innovation $e^{GPR}_t = u_{g,t}/\\sqrt{\\sigma_{gg}}$ added to the best-BIC satellite specification (control-function relaxation of the exclusion restriction $\\eta_t \\perp u_t$). Estimation is hierarchical: the innovation series and the conjugate Normal--Inverse-Gamma satellite posterior are recomputed for each retained BVAR posterior draw, so generated-regressor uncertainty is fully propagated. The baseline specification of the paper corresponds to $\\lambda = 0$.\n", file = con)
cat("\\end{tablenotes}\n\\end{table}\n", file = con)
close(con)

# -----------------------------------------------------------------------------
# 10. README and console summary
# -----------------------------------------------------------------------------

readme <- c(
  "# CORALACBN direct-channel (control-function) exercise",
  "",
  "Relaxation of the satellite exogeneity assumption eta_t independent of VAR innovations.",
  "The satellite is augmented with the standardized reduced-form GPR innovation:",
  "",
  "    Z_t = beta0 + beta' Y_t^(s) + lambda * e_t^GPR + eps_t",
  "",
  "The baseline of the paper is the restriction lambda = 0. Estimation is hierarchical",
  "(per BVAR posterior draw), so generated-regressor uncertainty is propagated.",
  "GIRF corrections are exact and closed-form: psi_Z(0) += lambda * delta_sd;",
  "s2(h) += lambda^2 + 2*lambda*c0 (baseline: all h; shocked: h >= 1, since e_t is",
  "fixed by the conditioning at h = 0), with c0 = beta' S_0 Sigma e_g / sqrt(sigma_gg).",
  "",
  "Generated by code/extra/coralacbn/03_coralacbn_direct_channel.R (standalone; not in run_all.R).",
  "",
  "Contents:",
  "  tables/lambda_posterior_summary.csv|.tex   posterior of lambda + exclusion test",
  "  tables/lambda_ols_hac_diagnostic.csv       frequentist HAC check (innovation at posterior median)",
  "  tables/satellite_restricted_vs_augmented.csv  coefficient comparison",
  "  irf/psiZ_bands_restricted_vs_augmented.csv    psi_Z bands, both models, both shocks",
  "  irf/pd_girf_bands_restricted_vs_augmented.csv Delta PD bands, both models, both shocks",
  "  figures/                                    lambda density + overlays (Z and PD; 1sd and 2001Q3)",
  "  direct_channel_results.rds                  all objects (draws, bands, summaries)"
)
writeLines(readme, file.path(dc_out_dir, "README.md"))

pk_res <- res_1sd$pd_res_bands[which.max(abs(median))]
pk_aug <- res_1sd$pd_aug_bands[which.max(abs(median))]

cat("\n--- Direct-channel summary -------------------------------------------\n")
cat(sprintf("lambda: posterior mean = %.4f, sd = %.4f, 90%% CI = [%.4f, %.4f]\n",
            lambda_summary$mean, lambda_summary$sd,
            lambda_summary$q05, lambda_summary$q95))
cat(sprintf("Pr(lambda < 0 | data) = %.3f   (negative = direct adverse GPR effect on credit)\n",
            lambda_summary$prob_negative))
cat(sprintf("HAC diagnostic: estimate = %.4f, t = %.2f, p = %.4f\n",
            hac_row$estimate, hac_row$t_value, hac_row$p_value))
cat(sprintf("Peak Delta PD (1 s.d.): restricted = %.3e (h=%d) | augmented = %.3e (h=%d)\n",
            pk_res$median, pk_res$horizon, pk_aug$median, pk_aug$horizon))
cat("Outputs written to: ", dc_out_dir, "\n", sep = "")
cat("------------------------------------------------------------------------\n")

invisible(TRUE)
