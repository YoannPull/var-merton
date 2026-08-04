# =============================================================================
#  code/review/02_conditional_normality.R
#
#  Coarse-review diagnostic, conditional version (companion to 01).
#
#  The closed-form PD-aR needs the CONDITIONAL Gaussianity of Z_{t+h}, which is
#  delivered by (i) the satellite error eta_t and (ii) the VAR reduced-form
#  innovations u_t -- NOT by the unconditional marginal of the reconstructed
#  factor (tested in 01, and biased there by persistence). This script tests the
#  right objects:
#
#     (a) eta_t : residual of the best-BIC satellite bridge  Z_t = X_t beta + eta_t
#     (b) u_t   : reduced-form innovations of the VAR(P) on the real-side set,
#                 full sample and excluding the 2020 COVID quarters.
#
#  Findings (full sample):
#     - eta_t is approximately Gaussian (Jarque-Bera does not reject; mild
#       platykurtosis only) -> the credit-side building block is fine.
#     - The real-activity VAR innovations (GDP p.c., private employment,
#       investment) are strongly non-normal ON THE FULL SAMPLE, but this is an
#       artefact of the 2020 pandemic quarters: excluding 2020Q1-Q4 they become
#       essentially Gaussian (JB p > 0.10). This is the well-known COVID-VAR
#       outlier problem (Lenza & Primiceri, 2022), addressable by COVID-volatility
#       scaling / dummies.
#     - GPR, oil and inflation innovations remain non-Gaussian even ex-COVID
#       (right-skew / fat tails intrinsic to those series). The shock of interest
#       (GPR) is treated at fixed sizes (1 s.d. and historical episodes), so its
#       innovation non-normality does not bias the response to a GIVEN shock.
#
#  STANDALONE: not in run_all.R. Run from the repository root:
#      source("code/review/02_conditional_normality.R")
#  Outputs: output/review/figures/conditional_normality.pdf
#           output/review/tables/conditional_normality_stats.csv
# =============================================================================

source("code/00_setup.R")
source("code/shared/_helpers_var_merton.R")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

REVIEW_FIG <- file.path(DIR_OUTPUT, "review", "figures")
REVIEW_TAB <- file.path(DIR_OUTPUT, "review", "tables")
dir.create(REVIEW_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(REVIEW_TAB, recursive = TRUE, showWarnings = FALSE)

P            <- P_LAGS
SAMPLE_START <- as.Date("1986-01-01")
SAMPLE_END   <- as.Date("2024-12-31")
COVID_Q      <- c(2020.0, 2020.25, 2020.5, 2020.75)

VARS6 <- c("log_GPRD", "log_inv_pc", "log_gdp_pc",
           "log_private_pc", "log_oil_real", "infl_yoy_pct")

# --- Normality helpers -------------------------------------------------------
moments_row <- function(name, x) {
  x  <- x[is.finite(x)]
  m  <- mean(x); s2 <- mean((x - m)^2)
  sk <- mean((x - m)^3) / s2^1.5
  ek <- mean((x - m)^4) / s2^2 - 3
  jb <- length(x) * (sk^2 / 6 + ek^2 / 24)
  data.table(series = name, n = length(x), skew = sk, excess_kurt = ek,
             shapiro_p = shapiro.test(x)$p.value,
             jb_p = stats::pchisq(jb, df = 2, lower.tail = FALSE))
}

# --- Reconstruct Z and the satellite residual (best-BIC bridge) --------------
dr <- fread(file.path(DIR_RAW, "default", "DRALACBN.csv"))
setnames(dr, c("observation_date", "DRALACBN"), c("date", "value"))
dr[, date := as.Date(date)]
dr <- dr[date >= SAMPLE_START & date <= SAMPLE_END]
Z  <- f_Z_estimation(dr$value / 100)$Z

vd <- fread(file.path(DIR_PROCESSED, "data_var_for_model.csv"))
V  <- as.matrix(vd[, ..VARS6])
n  <- nrow(V)
dq <- vd$Date_quarter

lagcol <- function(name, L) {
  x <- V[, name]
  if (L == 0) return(x)
  c(rep(NA_real_, L), x[seq_len(n - L)])
}
# best-BIC satellite specification (direct-channel restricted model)
X <- cbind(1,
           lagcol("log_inv_pc", 0), lagcol("log_inv_pc", 2),
           lagcol("log_gdp_pc", 2),
           lagcol("log_oil_real", 0), lagcol("log_oil_real", 4),
           lagcol("infl_yoy_pct", 4))
ok  <- stats::complete.cases(X)
eta <- as.numeric(Z[ok] - X[ok, ] %*% qr.solve(X[ok, ], Z[ok]))

# --- VAR(P) reduced-form innovations (GPR ordered first) ---------------------
Xv <- cbind(1, do.call(cbind, lapply(seq_len(P), function(l) V[(P - l + 1):(n - l), ])))
U  <- V[(P + 1):n, ] - Xv %*% qr.solve(Xv, V[(P + 1):n, ])
dqU    <- dq[(P + 1):n]
covid  <- round(dqU, 2) %in% COVID_Q

# --- Stats table -------------------------------------------------------------
stat_rows <- list(moments_row("satellite_eta", eta))
for (nm in VARS6) stat_rows[[paste0("full_", nm)]]   <- moments_row(paste0("u_", nm, "_full"),   U[, nm])
for (nm in VARS6) stat_rows[[paste0("ex20_", nm)]]   <- moments_row(paste0("u_", nm, "_ex2020"), U[!covid, nm])
stats_dt <- rbindlist(stat_rows)
fwrite(stats_dt, file.path(REVIEW_TAB, "conditional_normality_stats.csv"))
print(stats_dt)

# --- Figure: eta + the two COVID-driven innovations --------------------------
red <- "#8B0000"; grey <- "#464646"; blue <- "#1f4e79"
qq_panel <- function(x, title, highlight = NULL) {
  x <- as.numeric(x)
  o <- order(x)
  theo <- qnorm(ppoints(length(x)))
  d <- data.table(theo = theo, samp = sort(x),
                  hl = if (is.null(highlight)) FALSE else highlight[o])
  ggplot(d, aes(theo, samp)) +
    geom_abline(slope = sd(x), intercept = mean(x), linetype = "dashed", colour = grey) +
    geom_point(aes(colour = hl, size = hl)) +
    scale_colour_manual(values = c(`FALSE` = red, `TRUE` = blue), guide = "none") +
    scale_size_manual(values = c(`FALSE` = 1.4, `TRUE` = 2.6), guide = "none") +
    labs(x = "Normal quantiles", y = "Sample quantiles", title = title) +
    theme_minimal(base_size = 11)
}
fig <- qq_panel(eta, expression("(a) Satellite residual " * eta[t])) +
  qq_panel(U[, "log_gdp_pc"], "(b) VAR innovation: GDP p.c.", covid) +
  qq_panel(U[, "log_private_pc"], "(c) VAR innovation: Priv. employment", covid)
ggsave(file.path(REVIEW_FIG, "conditional_normality.pdf"), fig, width = 13, height = 4)

cat("\nConditional-normality diagnostic written to ", file.path(DIR_OUTPUT, "review"), "\n", sep = "")
