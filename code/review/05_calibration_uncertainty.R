# =============================================================================
#  code/review/05_calibration_uncertainty.R
#
#  Coarse-review point: the reported bands condition on the calibrated Merton
#  parameters (p_hat, rho_hat) and treat the reconstructed factor Z as observed,
#  so calibration error and generated-regressor error are excluded.
#
#  This script propagates both, holding the macro VAR at its point estimate (to
#  isolate the credit block). A moving-block bootstrap of the delinquency series
#  delivers the sampling distribution of (p_hat, rho_hat) (order-invariant
#  statistics p_hat = mean(d), rho_hat = V/(1+V), V = Var(Phi^{-1}(d))). Within
#  each draw Z is reconstructed with the bootstrapped calibration on the original,
#  date-aligned d, the satellite is re-estimated on the regenerated Z, and its
#  coefficient posterior is drawn. The result is compared to a calibration-fixed
#  band (satellite posterior only) for BOTH the projected PD levels and the
#  responses (mean, PD-aR, expected shortfall).
#
#  Finding: every reported conditional object -- the projected PD level path AND
#  the mean / PD-aR / expected-shortfall responses -- is EXACTLY invariant to the
#  calibration (p_hat, rho_hat). The reconstructed factor is affine in
#  Phi^{-1}(d), Z = a + b w with a = Phi^{-1}(p)/sqrt(rho), b = -sqrt(1-rho)/sqrt(rho);
#  re-estimating the satellite rescales the fitted mean and the conditional sd by
#  b, after which
#       Phi^{-1}(p) - sqrt(rho) mu = sqrt(1-rho) w_hat,
#       (1-rho) + rho s^2          = (1-rho)(1 + s_hat^2),
#  so the level Phi(w_hat / sqrt(1 + s_hat^2)) and (by differencing) the responses
#  are free of (p, rho). The calibration is a pure reparametrisation of the latent
#  factor, absorbed by the satellite re-estimation; the band sets coincide and
#  there is nothing to propagate.
#
#  Scope of the claim: this is invariance to the CALIBRATION CHOICE on the same
#  data (re-estimating the satellite). It is distinct from the fixed-loadings
#  comparative static "a portfolio with higher asset correlation amplifies the
#  tail more", which holds the factor loadings fixed and varies rho. The
#  unconditional through-the-cycle level (= p_hat by construction) is not the
#  object here; every CONDITIONAL reported object is invariant.
#
#  STANDALONE: not in run_all.R. Run from the repository root:
#      source("code/review/05_calibration_uncertainty.R")
#  Outputs:
#      output/review/figures/calibration_uncertainty_girf.pdf   (levels + responses)
#      output/review/figures/calibration_bootstrap.pdf          ((p_hat, rho_hat) draws)
#      output/review/tables/calibration_uncertainty_peaks.csv
# =============================================================================

source("code/00_setup.R")
source("code/shared/_helpers_var_merton.R")
source("code/shared/_pd_var.R")
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
set.seed(20260629)

REVIEW_FIG <- file.path(DIR_OUTPUT, "review", "figures")
REVIEW_TAB <- file.path(DIR_OUTPUT, "review", "tables")
dir.create(REVIEW_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(REVIEW_TAB, recursive = TRUE, showWarnings = FALSE)

H      <- if (exists("HORIZON")) HORIZON else 12L
P      <- P_LAGS
ALPHAS <- c(0.99, 0.999)
ND     <- 2000L                # posterior / bootstrap draws per band set
BLOCK  <- 12L                  # moving-block length for the delinquency bootstrap
VARS6  <- c("log_GPRD", "log_inv_pc", "log_gdp_pc",
            "log_private_pc", "log_oil_real", "infl_yoy_pct")
G_IDX  <- 1L

# -----------------------------------------------------------------------------
# 1. Fixed macro VAR(P) (OLS point estimate) and MA coefficients
# -----------------------------------------------------------------------------
vd <- fread(file.path(DIR_PROCESSED, "data_var_for_model.csv"))
Y  <- as.matrix(vd[, ..VARS6]); n <- nrow(Y); k <- ncol(Y)
Xd <- cbind(1, do.call(cbind, lapply(seq_len(P), function(l) Y[(P - l + 1):(n - l), ])))
Yt <- Y[(P + 1):n, ]
B  <- qr.solve(Xd, Yt)
Sigma <- crossprod(Yt - Xd %*% B) / (nrow(Yt) - ncol(Xd)); Sigma <- 0.5 * (Sigma + t(Sigma))
Psi    <- compute_ma_coefficients(B, P, H)
Y_fore <- forecast_baseline_path(B, Y[(n - P + 1):n, , drop = FALSE], P, H)

sgg      <- Sigma[G_IDX, G_IDX]
delta    <- Sigma[, G_IDX] / sqrt(sgg)
Sig_cond <- Sigma - tcrossprod(Sigma[, G_IDX], Sigma[, G_IDX]) / sgg
Sig_cond <- 0.5 * (Sig_cond + t(Sig_cond))

# -----------------------------------------------------------------------------
# 2. Delinquency proxy and the best-BIC satellite design
# -----------------------------------------------------------------------------
dr <- fread(file.path(DIR_RAW, "default", "DRALACBN.csv"))
setnames(dr, c("observation_date", "DRALACBN"), c("date", "value")); dr[, date := as.Date(date)]
dr <- dr[date >= as.Date("1986-01-01") & date <= as.Date("2024-12-31")]
d  <- pmin(pmax(dr$value / 100, 1e-8), 1 - 1e-8); Td <- length(d)

terms_df <- data.table(
  base = c("log_inv_pc", "log_inv_pc", "log_gdp_pc", "log_oil_real", "log_oil_real", "infl_yoy_pct"),
  lag  = c(0L, 2L, 2L, 0L, 4L, 4L))
terms_df[, col := match(base, VARS6)]
mterm <- nrow(terms_df); Lmax <- max(terms_df$lag)
sel   <- build_selection_from_terms(terms_df, VARS6); S_list <- sel$S_list

Xs <- cbind(1, do.call(cbind, lapply(seq_len(mterm), function(r)
  Y[(Lmax + 1 - terms_df$lag[r]):(Td - terms_df$lag[r]), terms_df$col[r]])))
nsat <- nrow(Xs); XtXinv <- solve(crossprod(Xs)); dfree <- nsat - (mterm + 1L)
Gload <- lapply(0:H, function(hh) lapply(0:hh, function(q) build_G_hq(hh, q, S_list, Lmax, Psi)))

reconstruct_Z <- function(p, rho) (qnorm(p) - sqrt(1 - rho) * qnorm(d)) / sqrt(rho)

# -----------------------------------------------------------------------------
# 3. Level path and responses from a calibration / satellite draw
# -----------------------------------------------------------------------------
objects_from <- function(p, rho, beta0, bvec, sigma_eta2) {
  hist_tail <- Y[(Td - Lmax):Td, , drop = FALSE]
  muZ <- rep(beta0, H + 1)
  for (r in seq_len(mterm)) {
    cr <- terms_df$col[r]; lr <- terms_df$lag[r]
    for (hh in 0:H) {
      val <- if (hh >= lr) Y_fore[hh - lr + 1, cr] else hist_tail[(Lmax + 1) - (lr - hh), cr]
      muZ[hh + 1] <- muZ[hh + 1] + bvec[r] * val
    }
  }
  s2 <- numeric(H + 1); s2d <- numeric(H + 1); psiZ <- numeric(H + 1)
  for (hh in 0:H) {
    for (q in 0:hh) {
      Bhq <- as.numeric(crossprod(bvec, Gload[[hh + 1]][[q + 1]]))
      s2[hh + 1]  <- s2[hh + 1]  + as.numeric(Bhq %*% Sigma %*% Bhq)
      s2d[hh + 1] <- s2d[hh + 1] + as.numeric(Bhq %*% (if (q == 0) Sig_cond else Sigma) %*% Bhq)
    }
    s2[hh + 1] <- s2[hh + 1] + sigma_eta2; s2d[hh + 1] <- s2d[hh + 1] + sigma_eta2
    psiZ[hh + 1] <- as.numeric(crossprod(bvec, Gload[[hh + 1]][[1]]) %*% delta)
  }
  col1 <- function(x) matrix(x, ncol = 1)
  # projected (baseline) PD level path, E[pi(Z_{t+h})]
  level <- pnorm((qnorm(p) - sqrt(rho) * muZ) / sqrt((1 - rho) + rho * s2))
  out <- list(level = level,
              mean = as.numeric(compute_pd_girf(col1(psiZ), col1(muZ), col1(s2), col1(s2d), p, rho)$draws))
  for (a in ALPHAS) {
    out[[sprintf("VaR%g", a)]] <- as.numeric(compute_pd_var_girf(col1(psiZ), col1(muZ), col1(s2), col1(s2d), p, rho, a)$draws)
    out[[sprintf("ES%g",  a)]] <- as.numeric(compute_pd_es_girf (col1(psiZ), col1(muZ), col1(s2), col1(s2d), p, rho, a)$draws)
  }
  out
}
draw_satellite <- function(Z) {
  y <- Z[(Lmax + 1):Td]; beta <- XtXinv %*% crossprod(Xs, y)
  s2 <- sum((y - Xs %*% beta)^2) / dfree
  sig <- 1 / rgamma(1, dfree / 2, rate = dfree * s2 / 2)
  bdraw <- as.numeric(beta) + as.numeric(crossprod(chol(sig * XtXinv), rnorm(mterm + 1L)))
  list(beta0 = bdraw[1], bvec = bdraw[-1], sigma_eta2 = sig)
}

# -----------------------------------------------------------------------------
# 4. Calibration point estimate and the two band sets
# -----------------------------------------------------------------------------
Vhat <- var(qnorm(d)); rho_hat <- Vhat / (1 + Vhat); p_hat <- mean(d)
cat(sprintf("p_hat = %.4f, rho_hat = %.4f\n", p_hat, rho_hat))

OBJS <- c("level", "mean", sprintf("VaR%g", ALPHAS), sprintf("ES%g", ALPHAS))
emptymat <- function() setNames(lapply(OBJS, function(o) matrix(NA_real_, ND, H + 1)), OBJS)

# (i) calibration fixed: only the satellite posterior is drawn
Z0 <- reconstruct_Z(p_hat, rho_hat); base <- emptymat()
for (i in seq_len(ND)) {
  s <- draw_satellite(Z0); r <- objects_from(p_hat, rho_hat, s$beta0, s$bvec, s$sigma_eta2)
  for (o in OBJS) base[[o]][i, ] <- r[[o]]
}
# (ii) calibration + generated regressor: moving-block bootstrap of d
prop <- emptymat(); pr <- matrix(NA_real_, ND, 2)
nblk <- ceiling(Td / BLOCK)
for (i in seq_len(ND)) {
  starts <- sample.int(Td - BLOCK + 1L, nblk, replace = TRUE)
  db <- unlist(lapply(starts, function(s) d[s:(s + BLOCK - 1L)]))[1:Td]
  pb <- mean(db); Vb <- var(qnorm(db)); rb <- Vb / (1 + Vb); pr[i, ] <- c(pb, rb)
  Zb <- reconstruct_Z(pb, rb); s <- draw_satellite(Zb)
  r <- objects_from(pb, rb, s$beta0, s$bvec, s$sigma_eta2)
  for (o in OBJS) prop[[o]][i, ] <- r[[o]]
}

# -----------------------------------------------------------------------------
# 5. Summary table: peak / end-of-horizon and 90% band widths
# -----------------------------------------------------------------------------
pk  <- which.max(abs(apply(base$mean, 2, median)))            # response peak horizon
w90 <- function(M, h) diff(quantile(M[, h] * 100, c(0.05, 0.95)))
summ <- rbindlist(lapply(OBJS, function(o) {
  h <- if (o == "level") H + 1L else pk          # level reported at the 3-year horizon
  data.table(object = o, horizon = h - 1L,
             median_pp_or_pct = 100 * median(base[[o]][, h]),
             band90_fixed = w90(base[[o]], h), band90_propagated = w90(prop[[o]], h),
             ratio = w90(prop[[o]], h) / w90(base[[o]], h))
}))
fwrite(summ, file.path(REVIEW_TAB, "calibration_uncertainty_peaks.csv")); print(summ)

# -----------------------------------------------------------------------------
# 6. Figure A: level + responses, calibration fixed (shaded) vs propagated (dashed)
# -----------------------------------------------------------------------------
labs <- c(level = "PD level (%)", mean = "Mean response (pp)",
          `VaR0.99` = "PD-aR 99% (pp)", `VaR0.999` = "PD-aR 99.9% (pp)",
          `ES0.99` = "ES 99% (pp)", `ES0.999` = "ES 99.9% (pp)")
qb <- function(M, p) apply(M * 100, 2, quantile, probs = p)
dt <- rbindlist(lapply(OBJS, function(o) data.table(
  horizon = 0:H, panel = labs[[o]],
  med_fix = qb(base[[o]], .5), lo_fix = qb(base[[o]], .05), hi_fix = qb(base[[o]], .95),
  lo_pro = qb(prop[[o]], .05), hi_pro = qb(prop[[o]], .95))))
dt[, panel := factor(panel, levels = labs)]
ggA <- ggplot(dt, aes(horizon)) +
  geom_ribbon(aes(ymin = lo_fix, ymax = hi_fix), fill = "#6F1732", alpha = 0.20) +
  geom_line(aes(y = med_fix), color = "#6F1732", linewidth = 0.8) +
  geom_line(aes(y = lo_pro), color = "#1F6FB2", linetype = "dashed", linewidth = 0.5) +
  geom_line(aes(y = hi_pro), color = "#1F6FB2", linetype = "dashed", linewidth = 0.5) +
  facet_wrap(~panel, scales = "free_y") +
  labs(x = "Horizon (quarters)", y = NULL,
       subtitle = "Shaded: calibration fixed (90%). Dashed: calibration + generated-regressor (90%).") +
  theme_minimal(base_size = 11)
ggsave(file.path(REVIEW_FIG, "calibration_uncertainty_girf.pdf"), ggA, width = 11, height = 6)

# -----------------------------------------------------------------------------
# 7. Figure B: bootstrap of (p_hat, rho_hat) -- the calibration does vary
# -----------------------------------------------------------------------------
prdt <- data.table(p = pr[, 1] * 100, rho = pr[, 2])
ggB <- ggplot(prdt, aes(p, rho)) +
  geom_point(color = "#555555", alpha = 0.35, size = 1) +
  geom_vline(xintercept = 100 * p_hat, color = "#6F1732") +
  geom_hline(yintercept = rho_hat, color = "#6F1732") +
  labs(x = expression(hat(p) ~ "(%)"), y = expression(hat(rho)),
       title = "Moving-block bootstrap of the Merton calibration") +
  theme_minimal(base_size = 11)
ggsave(file.path(REVIEW_FIG, "calibration_bootstrap.pdf"), ggB, width = 5, height = 4)

cat("\nCalibration-uncertainty propagation written to ", REVIEW_FIG, "\n", sep = "")
cat("Both band sets coincide (ratios near 1) for the level and every response:\n",
    "the projected levels and the responses are invariant to (p_hat, rho_hat),\n",
    "which the satellite re-estimation absorbs. There is nothing to propagate.\n", sep = "")
