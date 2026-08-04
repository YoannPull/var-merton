# =============================================================================
#  code/review/04_simulation_benchmark.R
#
#  Coarse-review point: "No numerical benchmark of the closed form against
#  default-event simulation."
#
#  This script confronts the closed-form responses of Section 2 (mean PD,
#  Corollary PD-at-Risk, Corollary expected shortfall) with a brute-force
#  forward Monte Carlo simulation of the system, at a fixed parameter vector.
#  The point is not to confirm correctness (the proofs already establish it) but
#  to demonstrate, and quantify, the payoff: the closed form returns the entire
#  conditional PD distribution, including the 99.9% tail, at no cost, whereas the
#  simulation needs a very large number of joint innovation draws and still
#  carries Monte Carlo error in the deep tail.
#
#  Design (rigorous, internally consistent).
#  * Fixed parameters theta: the VAR is the OLS / posterior-mean point estimate
#    of the real-side VAR(P) (GPR ordered first); the satellite uses the
#    posterior-mean coefficients of the Bayesian satellite. The benchmark is a
#    closed-form-vs-simulation comparison AT THIS theta, so it is agnostic to
#    posterior uncertainty (which is reported elsewhere).
#  * Both the closed form and the simulation are built from the SAME objects:
#    the factor follows, exactly,
#        Z_{t+h} = mu_{t+h} + sum_{q=0}^{h} b_{h,q}' u_{t+q} + eta_{t+h},
#        b_{h,q}' = beta' G_{h,q},   (G_{h,q} from the pipeline helper),
#    with u_{t+q} ~ N(0, Sigma_u), eta ~ N(0, sigma_eta^2). The closed-form
#    moments (mu, s^2, mu^shock, s^2,shock) and responses reuse the production
#    functions compute_pd_girf / compute_pd_var_girf / compute_pd_es_girf.
#  * The shock conditions the impact innovation on u_{g,t} = delta_g (a
#    one-standard-deviation GPR innovation). Common random numbers are used: the
#    baseline and shocked paths share eta and all future innovations, and share
#    the conditional component of u_t, differing only in the g-innovation at
#    impact. This mirrors the "identical posterior draws" treatment of the paper
#    and minimises the Monte Carlo error of the response (a difference).
#
#  The simulation targets the asymptotic single-risk-factor (Gordy) default
#  rate pi(Z); a finite-pool default-event simulation would add granularity
#  noise on top, only widening the gap the closed form removes.
#
#  PART OF THE PAPER PIPELINE: run_all.R executes this script (phase
#  "numerical_benchmark"), after the DRALACBN application, because
#  make_paper_outputs.R copies the figure below into the paper deliverables as
#  figures/main/simulation_benchmark.pdf (Figure of Section 2.5).
#  It can also be run on its own, from the repository root, once the DRALACBN
#  main specification exists:
#      source("code/review/04_simulation_benchmark.R")
#  Outputs: output/review/figures/simulation_benchmark.pdf
#           output/review/tables/simulation_benchmark.csv
# =============================================================================

source("code/00_setup.R")
source("code/shared/_helpers_var_merton.R")
source("code/shared/_pd_var.R")
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

set.seed(20260629)

REVIEW_FIG <- file.path(DIR_OUTPUT, "review", "figures")
REVIEW_TAB <- file.path(DIR_OUTPUT, "review", "tables")
dir.create(REVIEW_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(REVIEW_TAB, recursive = TRUE, showWarnings = FALSE)

H        <- if (exists("HORIZON")) HORIZON else 12L
P        <- P_LAGS
P_TTC    <- if (exists("APP_P_TTC")) APP_P_TTC else 0.032
RHO      <- if (exists("APP_RHO")) APP_RHO else 0.051
ALPHAS   <- c(0.99, 0.999)
N_SIM    <- 1e6L          # draws for the main comparison (raise for a tighter tail)
N_GRID   <- c(1e3, 1e4, 1e5, 1e6)   # convergence ladder for the 99.9% peak
VARS6    <- c("log_GPRD", "log_inv_pc", "log_gdp_pc",
             "log_private_pc", "log_oil_real", "infl_yoy_pct")
G_IDX    <- 1L            # GPR ordered first

# -----------------------------------------------------------------------------
# 1. Fixed parameters: OLS VAR(P) and posterior-mean satellite
# -----------------------------------------------------------------------------
vd <- fread(file.path(DIR_PROCESSED, "data_var_for_model.csv"))
Y  <- as.matrix(vd[, ..VARS6]); n <- nrow(Y); k <- ncol(Y)

# Design [1, Y_{t-1}, ..., Y_{t-P}]; B is (1 + kP) x k with intercept in row 1,
# matching the layout expected by build_companion / compute_ma_coefficients.
Xd <- cbind(1, do.call(cbind, lapply(seq_len(P), function(l) Y[(P - l + 1):(n - l), ])))
Yt <- Y[(P + 1):n, ]
B  <- qr.solve(Xd, Yt)                         # (1+kP) x k
U_resid <- Yt - Xd %*% B
Sigma   <- crossprod(U_resid) / (nrow(Yt) - ncol(Xd))
Sigma   <- 0.5 * (Sigma + t(Sigma))

Psi   <- compute_ma_coefficients(B, P, H)      # [H+1, k, k]
Y_hist <- Y[(n - P + 1):n, , drop = FALSE]
Y_fore <- forecast_baseline_path(B, Y_hist, P, H)

bsat  <- readRDS(file.path(DIR_DRALACBN_APP, "main_specification",
                           "satellite", "bayesian_satellite.rds"))
beta0 <- mean(bsat$beta0_draws)
beta  <- colMeans(bsat$beta_term_draws)        # posterior-mean mixture coefficients
sigma_eta2 <- mean(bsat$sigma2_draws)
terms_df <- as.data.table(bsat$terms)
terms_df[, beta := beta]

# -----------------------------------------------------------------------------
# 2. Closed-form conditional moments and satellite-weighted MA loadings
# -----------------------------------------------------------------------------
sel    <- build_selection_from_terms(terms_df, VARS6)
S_list <- sel$S_list; Lmax <- sel$Lmax
beta_vec <- terms_df$beta

# Baseline mean of Z: beta0 + sum_r beta_r * (forecast or history value).
hist_tail <- Y[(n - Lmax):n, , drop = FALSE]
muZ <- rep(beta0, H + 1)
for (r in seq_len(nrow(terms_df))) {
  cr <- terms_df$col[r]; lr <- terms_df$lag[r]; br <- terms_df$beta[r]
  for (hh in 0:H) {
    val <- if (hh >= lr) Y_fore[hh - lr + 1, cr] else hist_tail[(Lmax + 1) - (lr - hh), cr]
    muZ[hh + 1] <- muZ[hh + 1] + br * val
  }
}

# Conditional innovation covariance given u_{g}, and the impact mean shift delta.
sgg     <- Sigma[G_IDX, G_IDX]
delta   <- Sigma[, G_IDX] / sqrt(sgg)                 # E[u_t | u_gt = sqrt(sigma_gg)]
Sig_cond <- Sigma - tcrossprod(Sigma[, G_IDX], Sigma[, G_IDX]) / sgg
Sig_cond <- 0.5 * (Sig_cond + t(Sig_cond))

# Loadings b_{h,q} = beta' G_{h,q}, variances, and the factor GIRF psi_Z(h).
bload <- vector("list", H + 1)                        # bload[[h+1]][[q+1]] = numeric(k)
s2 <- numeric(H + 1); s2d <- numeric(H + 1); psiZ <- numeric(H + 1)
for (hh in 0:H) {
  bload[[hh + 1]] <- vector("list", hh + 1)
  for (q in 0:hh) {
    Bhq <- as.numeric(crossprod(beta_vec, build_G_hq(hh, q, S_list, Lmax, Psi)))  # length k
    bload[[hh + 1]][[q + 1]] <- Bhq
    s2[hh + 1]  <- s2[hh + 1]  + as.numeric(Bhq %*% Sigma %*% Bhq)
    s2d[hh + 1] <- s2d[hh + 1] + as.numeric(Bhq %*% (if (q == 0) Sig_cond else Sigma) %*% Bhq)
  }
  s2[hh + 1]  <- s2[hh + 1]  + sigma_eta2
  s2d[hh + 1] <- s2d[hh + 1] + sigma_eta2
  psiZ[hh + 1] <- as.numeric(bload[[hh + 1]][[1]] %*% delta)
}
s <- sqrt(s2); sd_ <- sqrt(s2d)

# -----------------------------------------------------------------------------
# 3. Closed-form responses (reuse production functions; single parameter column)
# -----------------------------------------------------------------------------
col1 <- function(x) matrix(x, ncol = 1)
psiZ_m <- col1(psiZ); mu_m <- col1(muZ); s2_m <- col1(s2); s2d_m <- col1(s2d)
cf_mean <- as.numeric(compute_pd_girf(psiZ_m, mu_m, s2_m, s2d_m, P_TTC, RHO)$draws)
cf <- list(mean = cf_mean)
for (a in ALPHAS) {
  cf[[sprintf("VaR%g", a)]] <-
    as.numeric(compute_pd_var_girf(psiZ_m, mu_m, s2_m, s2d_m, P_TTC, RHO, a)$draws)
  cf[[sprintf("ES%g", a)]] <-
    as.numeric(compute_pd_es_girf(psiZ_m, mu_m, s2_m, s2d_m, P_TTC, RHO, a)$draws)
}

# -----------------------------------------------------------------------------
# 4. Forward Monte Carlo simulation of the system (common random numbers)
# -----------------------------------------------------------------------------
Lf <- t(chol(Sigma))
Lc <- t(chol(Sig_cond + 1e-12 * diag(k)))
a_g <- Sigma[, G_IDX] / sgg                            # u_t = a_g * u_gt + eps

pi_map <- function(z) stats::pnorm((stats::qnorm(P_TTC) - sqrt(RHO) * z) / sqrt(1 - RHO))

# Returns the (H+1)-vector responses for mean / VaR_a / ES_a at a given N, plus
# a moment check at the peak horizon. Chunked to bound memory.
simulate_responses <- function(N, chunk = 2e5L) {
  qpi <- stats::qnorm(P_TTC)
  # accumulators: store PD columns to take quantiles at the end
  PDb <- matrix(0, N, H + 1); PDs <- matrix(0, N, H + 1)
  done <- 0L
  while (done < N) {
    nb <- min(chunk, N - done); idx <- (done + 1L):(done + nb)
    eps0 <- matrix(stats::rnorm(nb * k), nb, k) %*% t(Lc)      # shared conditional impact noise
    zb   <- stats::rnorm(nb, 0, sqrt(sgg))                     # baseline g-innovation
    eta  <- matrix(stats::rnorm(nb * (H + 1), 0, sqrt(sigma_eta2)), nb, H + 1)
    Uq   <- lapply(1:H, function(q) matrix(stats::rnorm(nb * k), nb, k) %*% t(Lf))  # shared futures
    u0b  <- outer(zb, a_g) + eps0                             # baseline impact ~ N(0, Sigma)
    Zb <- matrix(rep(muZ, each = nb), nb, H + 1)
    Zs <- matrix(rep(muZ + psiZ, each = nb), nb, H + 1)        # shock mean shift baked in
    for (hh in 0:H) {
      B0 <- bload[[hh + 1]][[1]]
      Zb[, hh + 1] <- Zb[, hh + 1] + as.numeric(u0b %*% B0)
      Zs[, hh + 1] <- Zs[, hh + 1] + as.numeric(eps0 %*% B0)   # shock impact: conditional noise only
      if (hh >= 1) for (q in 1:hh) {
        Bhq <- bload[[hh + 1]][[q + 1]]
        contrib <- as.numeric(Uq[[q]] %*% Bhq)
        Zb[, hh + 1] <- Zb[, hh + 1] + contrib
        Zs[, hh + 1] <- Zs[, hh + 1] + contrib
      }
      Zb[, hh + 1] <- Zb[, hh + 1] + eta[, hh + 1]
      Zs[, hh + 1] <- Zs[, hh + 1] + eta[, hh + 1]
    }
    PDb[idx, ] <- pi_map(Zb); PDs[idx, ] <- pi_map(Zs)
    done <- done + nb
  }
  out <- list()
  out$mean <- colMeans(PDs) - colMeans(PDb)
  out$pd_base_mean <- colMeans(PDb)    # simulated baseline mean PD level (for the check)
  for (a in ALPHAS) {
    qa_s <- apply(PDs, 2, stats::quantile, probs = a)
    qa_b <- apply(PDb, 2, stats::quantile, probs = a)
    out[[sprintf("VaR%g", a)]] <- qa_s - qa_b
    es_s <- sapply(1:(H + 1), function(h) mean(PDs[PDs[, h] >= qa_s[h], h]))
    es_b <- sapply(1:(H + 1), function(h) mean(PDb[PDb[, h] >= qa_b[h], h]))
    out[[sprintf("ES%g", a)]] <- es_s - es_b
  }
  out
}

# Closed-form baseline mean PD level, E[pi(Z_base)], for the consistency check.
pd_mean_level <- function(mu, sv) {
  stats::pnorm((stats::qnorm(P_TTC) - sqrt(RHO) * mu) / sqrt((1 - RHO) + RHO * sv^2))
}

# Consistency check (single run): the simulated baseline mean PD level must
# match the closed-form E[pi(Z_base)], validating the forward construction.
chk <- simulate_responses(2e5L)
stopifnot(max(abs(pd_mean_level(muZ, s) - chk$pd_base_mean)) < 5e-4)
cat(sprintf("Moment check OK: max |sim - closed-form| baseline mean PD level = %.2e\n",
            max(abs(pd_mean_level(muZ, s) - chk$pd_base_mean))))

# -----------------------------------------------------------------------------
# 5. Replication study across the draw ladder (Monte Carlo mean and SE)
#    Each N is simulated R_REPS times independently; we report, per object, the
#    mean response across replications (unbiasedness) and the standard deviation
#    across replications (Monte Carlo error of a single simulation of size N).
# -----------------------------------------------------------------------------
R_REPS <- 100L
pk <- which.max(abs(cf$mean))           # peak horizon index (1-based)
objs <- names(cf)

t0 <- Sys.time()
grid_res <- list()
for (N in N_GRID) {
  mats <- setNames(lapply(objs, function(nm) matrix(NA_real_, R_REPS, H + 1L)), objs)
  for (rr in seq_len(R_REPS)) {
    s <- simulate_responses(as.integer(N))
    for (nm in objs) mats[[nm]][rr, ] <- s[[nm]]
  }
  grid_res[[as.character(N)]] <- list(
    mean = setNames(lapply(objs, function(nm) colMeans(mats[[nm]])), objs),
    se   = setNames(lapply(objs, function(nm) apply(mats[[nm]], 2, sd)), objs)
  )
}
runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
big <- names(grid_res)[length(grid_res)]   # key of the largest grid N (robust to formatting)
stopifnot(!is.null(grid_res[[big]]))

# -----------------------------------------------------------------------------
# 6. Comparison at the peak: closed form vs simulation mean, with Monte Carlo SE
# -----------------------------------------------------------------------------
cmp <- rbindlist(lapply(objs, function(nm) data.table(
  object         = nm,
  horizon        = pk - 1L,
  closed_form_pp = 100 * cf[[nm]][pk],
  sim_mean_pp    = 100 * grid_res[[big]]$mean[[nm]][pk],
  mc_se_pp       = 100 * grid_res[[big]]$se[[nm]][pk],
  bias_pp        = 100 * (grid_res[[big]]$mean[[nm]][pk] - cf[[nm]][pk])
)))
cmp[, `:=`(n_sim = N_SIM, reps = R_REPS, runtime_sec = round(runtime, 1))]
fwrite(cmp, file.path(REVIEW_TAB, "simulation_benchmark.csv"))
print(cmp)

# -----------------------------------------------------------------------------
# 7. Convergence of the 99.9% peak response: Monte Carlo SE and RMSE vs N
# -----------------------------------------------------------------------------
conv <- rbindlist(lapply(N_GRID, function(N) {
  Nc <- as.character(N)
  bias_v <- grid_res[[Nc]]$mean[["VaR0.999"]][pk] - cf[["VaR0.999"]][pk]
  bias_e <- grid_res[[Nc]]$mean[["ES0.999"]][pk]  - cf[["ES0.999"]][pk]
  se_v <- grid_res[[Nc]]$se[["VaR0.999"]][pk]
  se_e <- grid_res[[Nc]]$se[["ES0.999"]][pk]
  data.table(N = N,
             mc_se_VaR999_pp = 100 * se_v,
             rmse_VaR999_pp  = 100 * sqrt(bias_v^2 + se_v^2),
             mc_se_ES999_pp  = 100 * se_e,
             rmse_ES999_pp   = 100 * sqrt(bias_e^2 + se_e^2))
}))
fwrite(conv, file.path(REVIEW_TAB, "simulation_benchmark_convergence.csv"))
print(conv)

# -----------------------------------------------------------------------------
# 7b. Timing: closed-form bands over B posterior draws versus the nested
#     simulation analogue. Producing the bands evaluates the response at every
#     posterior draw; the closed form does so analytically (B vectorised
#     evaluations), whereas a simulation-based band nests an inner Monte Carlo of
#     N draws inside each of the B draws, an O(B*N) computation.
# -----------------------------------------------------------------------------
B_DRAWS <- 1e4L                       # representative number of posterior draws
mu_B  <- matrix(muZ,  H + 1, B_DRAWS); s2_B  <- matrix(s2,  H + 1, B_DRAWS)
s2d_B <- matrix(s2d,  H + 1, B_DRAWS); psi_B <- matrix(psiZ, H + 1, B_DRAWS)
t_cf <- system.time({
  invisible(compute_pd_girf(psi_B, mu_B, s2_B, s2d_B, P_TTC, RHO))
  for (a in ALPHAS) {
    invisible(compute_pd_var_girf(psi_B, mu_B, s2_B, s2d_B, P_TTC, RHO, a))
    invisible(compute_pd_es_girf (psi_B, mu_B, s2_B, s2d_B, P_TTC, RHO, a))
  }
})[["elapsed"]]
t_inner   <- system.time(invisible(simulate_responses(N_SIM)))[["elapsed"]]
nested_s  <- B_DRAWS * t_inner        # B independent inner simulations
fwrite(data.table(B = B_DRAWS, N = N_SIM,
                  closed_form_bands_sec = t_cf,
                  one_inner_sim_sec = t_inner,
                  nested_sim_sec = nested_s,
                  speedup = nested_s / max(t_cf, 1e-6)),
       file.path(REVIEW_TAB, "simulation_benchmark_timing.csv"))
cat(sprintf(
  "\nTiming (B=%s posterior draws, N=%s):\n  closed-form bands:        %.2f s\n  one inner simulation:     %.2f s\n  nested simulation (B x N): %.0f s (~%.1f h)\n  speedup factor:           %.0fx\n",
  format(B_DRAWS, big.mark = ","), format(N_SIM, big.mark = ","),
  t_cf, t_inner, nested_s, nested_s / 3600, nested_s / max(t_cf, 1e-6)))

# -----------------------------------------------------------------------------
# 8. Figure: closed-form response paths (lines) vs simulation mean (points)
#    with +/- 2 Monte Carlo SE error bars, at the display N.
# -----------------------------------------------------------------------------
mk <- function(nm, lab) data.table(
  horizon = 0:H, closed = 100 * cf[[nm]],
  sim = 100 * grid_res[[big]]$mean[[nm]],
  se  = 100 * grid_res[[big]]$se[[nm]], series = lab)
plot_dt <- rbindlist(list(
  mk("mean",    "Mean"),
  mk("ES0.99",  "ES 99%"),
  mk("ES0.999", "ES 99.9%")
))
plot_dt[, series := factor(series, levels = c("Mean", "ES 99%", "ES 99.9%"))]
gg <- ggplot(plot_dt, aes(horizon, color = series)) +
  geom_line(aes(y = closed), linewidth = 0.9) +
  geom_errorbar(aes(ymin = sim - 2 * se, ymax = sim + 2 * se), width = 0.25) +
  geom_point(aes(y = sim), shape = 1, size = 1.8, stroke = 0.7) +
  scale_color_manual(values = c("Mean" = "#6F1732", "ES 99%" = "#1F6FB2",
                                "ES 99.9%" = "#E08214")) +
  labs(x = "Horizon (quarters)", y = expression(Delta ~ PD ~ "(pp)"), color = NULL) +
  theme_minimal(base_size = 11) + theme(legend.position = "top")
ggsave(file.path(REVIEW_FIG, "simulation_benchmark.pdf"), gg, width = 7, height = 4.2)

cat(sprintf("\nForward simulation benchmark written to %s\n", file.path(DIR_OUTPUT, "review")))
cat(sprintf("Grid %s, %d reps each; runtime %.1f s. The closed form returns every object exactly, at zero variance.\n",
            paste(format(N_GRID, scientific = TRUE), collapse = ", "), R_REPS, runtime))
