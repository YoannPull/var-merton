# =============================================================================
#  tests/test_pd_es_closedform.R — Expected-shortfall closed form unit test.
#
#  Run from the repository root:
#     Rscript tests/test_pd_es_closedform.R
#
#  Cross-checks the Simpson-quadrature implementation .pd_es_level() (used by
#  compute_pd_es_girf in code/shared/_pd_var.R) against:
#    (i)  the bivariate-normal closed form of the Expected-shortfall corollary,
#         ES_alpha = Phi2(m, qnorm(1-alpha); rho_star) / (1-alpha),
#         with m = (qnorm(p) - sqrt(rho) mu)/sqrt((1-rho)+rho s^2),
#              rho_star = sqrt(rho) s / sqrt((1-rho)+rho s^2),
#         evaluated with mvtnorm::pmvnorm;
#    (ii) a brute-force Monte Carlo of E[pi(Z) | Z <= z_alpha].
#
#  Exits with status 1 if any case disagrees beyond tolerance.
# =============================================================================

if (!file.exists("config.R")) stop("Run from the repository root.")
source("config.R")
source("code/shared/_pd_var.R")
suppressPackageStartupMessages(library(mvtnorm))

es_phi2 <- function(mu, s, p, rho, alpha) {
  den <- sqrt((1 - rho) + rho * s^2)
  m   <- (qnorm(p) - sqrt(rho) * mu) / den
  rs  <- sqrt(rho) * s / den
  cr  <- matrix(c(1, rs, rs, 1), 2)
  as.numeric(mvtnorm::pmvnorm(upper = c(m, qnorm(1 - alpha)), corr = cr)) / (1 - alpha)
}

es_mc <- function(mu, s, p, rho, alpha, N = 4e6) {
  Z  <- rnorm(N, mu, s)
  PD <- pnorm((qnorm(p) - sqrt(rho) * Z) / sqrt(1 - rho))
  zt <- mu + s * qnorm(1 - alpha)         # P(Z <= zt) = 1 - alpha
  mean(PD[Z <= zt])
}

set.seed(1)
cases <- list(
  c(p = 0.032, rho = 0.051, mu =  0.0, s = 1.0, alpha = 0.99),
  c(p = 0.032, rho = 0.051, mu = -0.3, s = 0.8, alpha = 0.99),
  c(p = 0.032, rho = 0.051, mu =  0.2, s = 1.2, alpha = 0.999),
  c(p = 0.050, rho = 0.150, mu =  0.0, s = 1.0, alpha = 0.975),
  c(p = 0.032, rho = 0.051, mu = -0.5, s = 0.6, alpha = 0.999)
)

n_fail <- 0L
for (cc in cases) {
  quad <- .pd_es_level(cc["mu"], cc["s"], cc["p"], cc["rho"], cc["alpha"])
  phi2 <- es_phi2(cc["mu"], cc["s"], cc["p"], cc["rho"], cc["alpha"])
  mc   <- es_mc(cc["mu"], cc["s"], cc["p"], cc["rho"], cc["alpha"])
  d_phi2 <- abs(quad - phi2)
  d_mc   <- abs(quad - mc)
  ok <- (d_phi2 < 1e-6) && (d_mc < 5e-3)
  if (!ok) n_fail <- n_fail + 1L
  cat(sprintf("%s  quad=%.6f  phi2=%.6f (d=%.1e)  MC=%.6f (d=%.1e)\n",
              if (ok) " ok  " else "FAIL ", quad, phi2, d_phi2, mc, d_mc))
}

# Limit checks: ES at alpha->0 equals the mean PD; ES >= VaR.
mu <- -0.2; s <- 0.9; p <- 0.032; rho <- 0.051
mean_pd <- pnorm((qnorm(p) - sqrt(rho) * mu) / sqrt((1 - rho) + rho * s^2))
es0 <- .pd_es_level(mu, s, p, rho, 1e-9)
var99 <- pnorm((qnorm(p) - sqrt(rho) * mu + sqrt(rho) * qnorm(0.99) * s) / sqrt(1 - rho))
es99  <- .pd_es_level(mu, s, p, rho, 0.99)
cat(sprintf("%s  ES(alpha->0)=%.6f vs meanPD=%.6f\n",
            if (abs(es0 - mean_pd) < 1e-6) " ok  " else "FAIL ", es0, mean_pd))
cat(sprintf("%s  ES99=%.6f >= VaR99=%.6f\n",
            if (es99 >= var99 - 1e-9) " ok  " else "FAIL ", es99, var99))
if (abs(es0 - mean_pd) >= 1e-6) n_fail <- n_fail + 1L
if (es99 < var99 - 1e-9) n_fail <- n_fail + 1L

if (n_fail > 0L) { cat("FAILED:", n_fail, "\n"); quit(status = 1L) }
cat("All expected-shortfall closed-form checks passed.\n")
