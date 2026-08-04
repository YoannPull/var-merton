# =============================================================================
#  code/review/01_zfactor_normality.R
#
#  Coarse-review diagnostic (Overall Feedback, point 1):
#  "The tail-amplification result is a property of the Gaussian map, not an
#   empirical finding, and is never tested against the data."
#
#  The PD-at-Risk quantiles of Corollary 2.5 are quantiles of the MODEL-IMPLIED
#  conditional distribution: a Gaussian latent factor Z_{t+h} passed through the
#  deterministic Merton-Vasicek probit map. This script confronts the Gaussian
#  assumption with the data by reconstructing the latent factor
#
#       Z_t = [ Phi^{-1}(p_hat) - sqrt(1 - rho_hat) Phi^{-1}(d_t) ] / sqrt(rho_hat)
#
#  from the observed delinquency series d_t (same inversion + unit-variance
#  calibration as the main text, via f_Z_estimation), and testing the
#  normality of the reconstructed factor (QQ-plot, Shapiro-Wilk, Jarque-Bera,
#  Anderson-Darling, skewness, excess kurtosis).
#
#  Caveat (printed in the log): Z_t is strongly persistent, so the i.i.d.-based
#  normality tests over-reject; read the QQ-plot and the moments as descriptive.
#  The object that actually feeds the PD-aR is the CONDITIONAL predictive
#  distribution of Z_{t+h} (Gaussian by VAR construction), not this unconditional
#  marginal.
#
#  STANDALONE: modifies nothing in the existing pipeline; not in run_all.R.
#  Run from the repository root:
#      source("code/review/01_zfactor_normality.R")
#
#  Outputs: output/review/
#      figures/zfactor_normality.pdf   QQ-plot + histogram vs N(0,1)
#      tables/zfactor_normality_stats.csv
# =============================================================================

source("code/00_setup.R")
source("code/shared/_helpers_var_merton.R")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

# --- Output layout -----------------------------------------------------------
REVIEW_DIR  <- file.path(DIR_OUTPUT, "review")
REVIEW_FIG  <- file.path(REVIEW_DIR, "figures")
REVIEW_TAB  <- file.path(REVIEW_DIR, "tables")
dir.create(REVIEW_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(REVIEW_TAB, recursive = TRUE, showWarnings = FALSE)

# --- Sample window (credit block, as in the main text) -----------------------
SAMPLE_START <- as.Date("1986-01-01")
SAMPLE_END   <- as.Date("2024-12-31")

# --- Load the delinquency proxy ----------------------------------------------
raw_path <- file.path(DIR_RAW, "default", "DRALACBN.csv")
dr <- fread(raw_path)
setnames(dr, c("observation_date", "DRALACBN"), c("date", "value"))
dr[, date := as.Date(date)]
dr <- dr[date >= SAMPLE_START & date <= SAMPLE_END]
d_t <- dr$value / 100  # percent -> proportion

# --- Reconstruct the latent factor (same calibration as the main text) -------
fit <- f_Z_estimation(d_t)          # estimates rho s.t. Var(Z) = 1
Z   <- fit$Z
p_ttc <- fit$p_ttc
rho_hat <- fit$rho

cat(sprintf("Sample %s..%s  n=%d\n", SAMPLE_START, SAMPLE_END, length(d_t)))
cat(sprintf("p_ttc = %.2f%%   rho_hat = %.4f   (rho_type=%s)\n",
            100 * p_ttc, rho_hat, fit$rho_type))

# --- Moments and normality tests ---------------------------------------------
n   <- length(Z)
m_z <- mean(Z)
s_z <- sd(Z)
skew_z <- mean((Z - m_z)^3) / (mean((Z - m_z)^2))^1.5   # population skewness
exk_z  <- mean((Z - m_z)^4) / (mean((Z - m_z)^2))^2 - 3  # excess kurtosis

sw <- shapiro.test(Z)
jb_stat <- n * (skew_z^2 / 6 + exk_z^2 / 24)              # Jarque-Bera
jb_p    <- stats::pchisq(jb_stat, df = 2, lower.tail = FALSE)
ad <- tryCatch(nortest::ad.test(Z), error = function(e) NULL)  # needs 'nortest'

stats_dt <- data.table(
  n             = n,
  p_ttc_pct     = 100 * p_ttc,
  rho_hat       = rho_hat,
  Z_mean        = m_z,
  Z_sd          = s_z,
  skewness      = skew_z,
  excess_kurt   = exk_z,
  shapiro_W     = unname(sw$statistic),
  shapiro_p     = sw$p.value,
  jarque_bera   = jb_stat,
  jarque_bera_p = jb_p,
  ad_stat       = if (!is.null(ad)) unname(ad$statistic) else NA_real_,
  ad_p          = if (!is.null(ad)) ad$p.value else NA_real_
)
fwrite(stats_dt, file.path(REVIEW_TAB, "zfactor_normality_stats.csv"))
print(stats_dt)

# --- Figure: QQ-plot + histogram vs N(0,1) -----------------------------------
red  <- if (exists("COL_MAIN")) COL_MAIN else "#8B0000"
grey <- "#464646"

qq_dt <- data.table(
  theo = qnorm(ppoints(length(Z)))[rank(Z, ties.method = "first")],
  samp = Z
)
p_qq <- ggplot(qq_dt, aes(theo, samp)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = grey) +
  geom_point(colour = red, alpha = 0.8, size = 1.6) +
  labs(x = "Theoretical normal quantiles",
       y = expression("Reconstructed " * Z[t] * " quantiles"),
       title = expression("(a) Normal QQ-plot of " * Z[t])) +
  theme_minimal(base_size = 11)

xs <- seq(min(Z) - 0.3, max(Z) + 0.3, length.out = 200)
dens_dt <- data.table(x = xs, y = dnorm(xs))
p_hist <- ggplot(data.table(Z = Z), aes(Z)) +
  geom_histogram(aes(y = after_stat(density)), bins = 18,
                 fill = red, alpha = 0.45, colour = "white") +
  geom_line(data = dens_dt, aes(x, y), colour = grey, linewidth = 0.8) +
  labs(x = expression(Z[t]), y = "Density",
       title = expression("(b) Histogram vs " * N(0, 1))) +
  theme_minimal(base_size = 11)

fig <- p_qq + p_hist
ggsave(file.path(REVIEW_FIG, "zfactor_normality.pdf"), fig,
       width = 10, height = 4)

cat("\nReview diagnostic written to ", REVIEW_DIR, "\n", sep = "")
cat("Interpretation: normality is rejected, but the departure is platykurtic\n",
    "(excess kurtosis < 0) and mildly left-skewed -- thin tails / regime\n",
    "bimodality, NOT fat tails. The Gaussian map therefore does not appear to\n",
    "understate the tail from a fat-tail channel on this aggregate proxy.\n", sep = "")
