# =============================================================================
#  code/extra/dralacbn/11_dralacbn_perfect_foresight_gulf_war.R
#
#  The cost of perfect foresight, evaluated at the most adverse geopolitical
#  state: the Gulf War (1990:Q3).
#
#  This script combines two exercises of the paper:
#    * the perfect-foresight bias (08_dralacbn_perfect_foresight.R), which
#      compares the exact closed forms of Proposition 3 to the supervisory
#      convention of treating the macro path as known; and
#    * the state dependence of the PD response (10_dralacbn_state_dependence.R),
#      which conditions the response on the credit state preceding a historical
#      episode.
#
#  The macro block is linear, so psi_Z and the conditional variances are
#  history-invariant; only the baseline mean path mu_{t+h} changes with the
#  conditioning history (cf. 10_*). We therefore condition mu on the history
#  ending the quarter BEFORE 1990:Q3 (the stressed credit state, delinquency
#  ~5.3%) and propagate the episode's OWN historical GPR shock through the
#  three variance variants:
#    exact     full closed form (Prop. 3): s2 = macro forecast variance + eta
#    pf_macro  macro path known, satellite noise retained: s2 = sigma_eta^2
#              (closest to the EBA perfect-foresight convention)
#    pf_full   full plug-in: s2 = 0 (PD evaluated at the mean path)
#
#  Object of interest: the SHOCKED PD LEVEL path along the Gulf-War scenario,
#  i.e. the absolute PD a supervisor would book under the adverse scenario
#  (not the response). pi is convex for small p, so perfect foresight
#  UNDERSTATES the expected PD level; the distortion grows with the horizon.
#  This is the most policy-relevant comparison because the Gulf-War state
#  carries the largest responses in the sample.
#
#  STANDALONE: modifies nothing in the existing pipeline; not in run_all.R
#  unless added. Requires the main DRALACBN application (02) to have run.
#  Run from the repository root:
#      source("code/extra/dralacbn/11_dralacbn_perfect_foresight_gulf_war.R")
#
#  Outputs: output/applications/dralacbn/perfect_foresight_gulf_war/
# =============================================================================

source("code/00_setup.R")
source("code/shared/_helpers_var_merton.R")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

# -----------------------------------------------------------------------------
# 0. Parameters, episode, output tree
# -----------------------------------------------------------------------------

PF_INFO_SET  <- if (exists("DELINQUENCY_BASELINE_INFO_SET")) DELINQUENCY_BASELINE_INFO_SET else "real_side"
p_lags       <- if (exists("P_LAGS")) P_LAGS else 2
nrep         <- if (exists("NREP")) NREP else 30000
M_target     <- if (exists("M_TARGET")) M_TARGET else 10000L
H            <- if (exists("HORIZON")) HORIZON else 12
seed         <- if (exists("SEED_ROBUSTNESS")) SEED_ROBUSTNESS else 123
impulse_name <- if (exists("IMPULSE_NAME")) IMPULSE_NAME else "log_GPRD"

GULF_QUARTER <- "1990Q3"
GULF_LABEL   <- "Gulf War (1990:Q3)"

pf_out_dir <- file.path("output", "applications", "dralacbn",
                        "perfect_foresight_gulf_war")
for (sub in c("", "tables", "figures", "irf")) {
  dir.create(file.path(pf_out_dir, sub), recursive = TRUE, showWarnings = FALSE)
}

VARIANT_LEVELS <- c("Exact (closed form)",
                    "PF: macro path known",
                    "PF: full plug-in")
VARIANT_COLORS <- c("Exact (closed form)" = "#6F1732",
                    "PF: macro path known" = "#AB4A7D",
                    "PF: full plug-in" = "#6B6B6B")
VARIANT_LTY <- c("Exact (closed form)" = "solid",
                 "PF: macro path known" = "dashed",
                 "PF: full plug-in" = "dotted")

theme_paper <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      axis.title = element_text(color = "#1A1A1A"),
      axis.text = element_text(color = "grey20"),
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linetype = "dotted", linewidth = 0.25),
      plot.caption = element_blank()
    )
}

save_fig <- function(plot, stem, width, height, dpi = 300) {
  ggsave(file.path(pf_out_dir, "figures", paste0(stem, ".png")),
         plot, width = width, height = height, dpi = dpi)
  tryCatch(
    ggsave(file.path(pf_out_dir, "figures", paste0(stem, ".pdf")),
           plot, width = width, height = height, device = grDevices::pdf),
    error = function(e) tryCatch(
      ggsave(file.path(pf_out_dir, "figures", paste0(stem, ".pdf")),
             plot, width = width, height = height),
      error = function(e2) message("PDF export failed for ", stem)
    )
  )
}

fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)

cat("\n===== DRALACBN: perfect foresight at the Gulf-War state (1990:Q3) =====\n")

# -----------------------------------------------------------------------------
# 1. Inputs: Merton parameters, data, kernel, Bayesian satellite
# -----------------------------------------------------------------------------

delinquency_path <- file.path(
  if (exists("DIR_RAW")) DIR_RAW else "data/raw", "default", "DRALACBN.csv"
)
delinq_raw <- fread(delinquency_path)
date_col <- intersect(c("DATE", "Date", "date", "observation_date"), names(delinq_raw))[1]
val_col  <- intersect(c("DRALACBN", "value", "VALUE", "Value"), names(delinq_raw))[1]
if (is.na(date_col)) date_col <- names(delinq_raw)[1]
if (is.na(val_col))  val_col  <- setdiff(names(delinq_raw), date_col)[1]
delinq <- data.table(Date = as.Date(delinq_raw[[date_col]]),
                     value = as.numeric(delinq_raw[[val_col]]))
delinq <- delinq[is.finite(value)]
setorder(delinq, Date)
if (exists("DATA_END_DATE")) delinq <- delinq[Date <= DATA_END_DATE]
delinq[, quarter := make_quarter(Date)]

z_out <- f_Z_estimation(delinq$value / 100)
P_TTC <- z_out$p_ttc
RHO   <- z_out$rho
cat("Merton parameters: p_ttc = ", sprintf("%.6f", P_TTC),
    ", rho = ", sprintf("%.6f", RHO), "\n", sep = "")

DT_raw <- safe_fread(if (exists("PATH_DATA_VAR")) PATH_DATA_VAR else "data/processed/data_var_for_model.csv")
DT_raw[, Date := as.Date(Date)]

kern <- get_or_estimate_kernel(
  spec_name = PF_INFO_SET,
  vars_extra = VAR_SPECS[[PF_INFO_SET]],
  DT_raw = DT_raw,
  p_lags = p_lags, nrep = nrep, H = H, seed = seed, M_target = M_target,
  kernel_dir = if (exists("DIR_VAR_KERNELS")) DIR_VAR_KERNELS else file.path("output", "shared", "var_kernels"),
  impulse_name = impulse_name
)
variables  <- kern$variables
impulse_ix <- match(impulse_name, variables)
if (is.na(impulse_ix)) stop("Impulse variable not found in VAR kernel.", call. = FALSE)
M <- dim(kern$Psi_draws)[4]

bsat_path <- file.path("output", "applications", "dralacbn",
                       "main_specification", "satellite", "bayesian_satellite.rds")
if (!file.exists(bsat_path)) {
  stop("Missing ", bsat_path,
       "\nRun the main DRALACBN application (02) first.", call. = FALSE)
}
bsat <- readRDS(bsat_path)
terms_df <- bsat$terms
n_sat <- nrow(bsat$beta_term_draws)
Lmax  <- max(terms_df$lag)
cat("Bayesian satellite loaded (", n_sat, " draws, main specification).\n", sep = "")

# Row dates of the kernel (complete cases, in order) -> quarter alignment.
cc_mask <- stats::complete.cases(as.data.frame(DT_raw[, ..variables]))
dates_kept <- as.Date(DT_raw$Date)[cc_mask]
if (length(dates_kept) != nrow(kern$DT)) {
  stop("Date alignment failed between DT_raw and the kernel.", call. = FALSE)
}
quarters_kept <- make_quarter(dates_kept)
Y_all <- as.matrix(kern$DT[, ..variables])
Tn <- nrow(Y_all)

# -----------------------------------------------------------------------------
# 2. History-invariant objects (computed once on the full posterior)
#    s2 / s2_delta (conditional variances) and psi_Z do NOT depend on the
#    conditioning history; only the baseline mean path mu does (cf. 10_*).
# -----------------------------------------------------------------------------

cat("Computing conditional variances and psi_Z (history-invariant) ...\n")
mom <- compute_factor_moment_draws_bayes(kern, bsat, impulse_ix, p_lags, H)
psiZ_1sd <- inject_girf_into_Z_bayes(kern$Psi_draws, kern$Sigma_kept, bsat,
                                     impulse_ix, shock_scale = 1)

# Satellite-error variance per BVAR draw (engine pairing convention).
sigma2_paired <- bsat$sigma2_draws[vapply(seq_len(M), pair_satellite_draw,
                                          integer(1), n_sat = n_sat)]
zero_mat <- matrix(0.0, nrow = H + 1, ncol = M)
eta_mat  <- matrix(rep(sigma2_paired, each = H + 1), nrow = H + 1, ncol = M)

variants <- list(
  exact    = list(s2 = mom$s2_draws,  s2d = mom$s2_delta_draws, label = VARIANT_LEVELS[1]),
  pf_macro = list(s2 = eta_mat,       s2d = eta_mat,            label = VARIANT_LEVELS[2]),
  pf_full  = list(s2 = zero_mat,      s2d = zero_mat,           label = VARIANT_LEVELS[3])
)

# -----------------------------------------------------------------------------
# 3. Gulf-War state: conditioning history and own historical shock size
# -----------------------------------------------------------------------------

history_end_row <- function(shock_quarter) {
  idx <- match(shock_quarter, quarters_kept)
  if (is.na(idx)) stop("Quarter not in sample: ", shock_quarter, call. = FALSE)
  if (idx - 1L < max(p_lags, Lmax + 1L)) {
    stop("Not enough history before ", shock_quarter, call. = FALSE)
  }
  idx - 1L
}

compute_mu_draws_at <- function(end_row) {
  Y_hist <- Y_all[(end_row - p_lags + 1L):end_row, , drop = FALSE]
  hist_tail <- Y_all[(end_row - Lmax):end_row, , drop = FALSE]
  mu_draws <- matrix(NA_real_, nrow = H + 1L, ncol = M)
  for (m in seq_len(M)) {
    sm <- pair_satellite_draw(m, n_sat)
    beta_m <- as.numeric(bsat$beta_term_draws[sm, terms_df$term, drop = TRUE])
    beta0_m <- bsat$beta0_draws[sm]
    Y_fore <- forecast_baseline_path(kern$B_kept[, , m], Y_hist, p_lags, H)
    mu_h <- rep(beta0_m, H + 1L)
    for (r in seq_len(nrow(terms_df))) {
      base_j <- terms_df$col[r]
      lag_r <- terms_df$lag[r]
      for (hh in 0:H) {
        val <- if (hh >= lag_r) {
          Y_fore[hh - lag_r + 1L, base_j]
        } else {
          hist_tail[(Lmax + 1L) - (lag_r - hh), base_j]
        }
        mu_h[hh + 1L] <- mu_h[hh + 1L] + beta_m[r] * val
      }
    }
    mu_draws[, m] <- mu_h
  }
  mu_draws
}

gw_end_row <- history_end_row(GULF_QUARTER)
mu_gw <- compute_mu_draws_at(gw_end_row)

gw_shock <- compute_structural_e1_median(
  kernel = kern, dates_raw = dates_kept,
  target_quarter = GULF_QUARTER, impulse_idx = impulse_ix
)
gw_size <- gw_shock$shock_scale
psiZ_gw <- psiZ_1sd$draws * gw_size

d_obs_gw <- delinq[quarter == GULF_QUARTER, value]
cat("Gulf-War state: history ends ", quarters_kept[gw_end_row],
    "; observed delinquency = ",
    if (length(d_obs_gw)) sprintf("%.2f%%", d_obs_gw[1]) else "NA",
    "; own shock size = ", sprintf("%.3f s.d.", gw_size), "\n", sep = "")

# -----------------------------------------------------------------------------
# 4. Shocked PD level and PD GIRF under each variant (Gulf-War scenario, pp)
# -----------------------------------------------------------------------------

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

# Shocked PD LEVEL along the Gulf-War scenario (same pd_shock convention as
# compute_pd_girf): pd_shock = Phi((qpi - sqrt(rho)*(mu + psiZ)) / s_shock).
pd_shock_level_draws <- function(mu_mat, psiZ_mat, s2d_mat) {
  qpi <- qnorm(P_TTC)
  mu_shock <- mu_mat + psiZ_mat
  denom <- sqrt((1 - RHO) + RHO * pmax(s2d_mat, 0))
  100 * pnorm((qpi - sqrt(RHO) * mu_shock) / denom)   # percentage points
}

cat("Computing shocked PD levels and PD GIRFs under the three variants ...\n")

shocked_levels_bands <- rbindlist(lapply(names(variants), function(v) {
  vv <- variants[[v]]
  bands_from_draws(pd_shock_level_draws(mu_gw, psiZ_gw, vv$s2d),
                   list(variant = vv$label))
}))
shocked_levels_bands[, variant := factor(variant, levels = VARIANT_LEVELS)]

girf_bands <- rbindlist(lapply(names(variants), function(v) {
  vv <- variants[[v]]
  pd <- compute_pd_girf(psiZ_gw, mu_gw, vv$s2, vv$s2d, p = P_TTC, rho = RHO)
  bands_from_draws(100 * pd$draws, list(variant = vv$label))
}))
girf_bands[, variant := factor(variant, levels = VARIANT_LEVELS)]

fwrite(shocked_levels_bands,
       file.path(pf_out_dir, "irf", "pd_shocked_level_bands_by_variant_gulf_war_pp.csv"))
fwrite(girf_bands,
       file.path(pf_out_dir, "irf", "pd_girf_bands_by_variant_gulf_war_pp.csv"))

# -----------------------------------------------------------------------------
# 5. Bias summaries (relative to the exact variant)
# -----------------------------------------------------------------------------

lev <- dcast(shocked_levels_bands[, .(variant, horizon, median)],
             horizon ~ variant, value.var = "median")
setnames(lev, VARIANT_LEVELS, c("exact", "pf_macro", "pf_full"))
lev[, bias_pf_macro_pct := 100 * (pf_macro - exact) / exact]
lev[, bias_pf_full_pct  := 100 * (pf_full  - exact) / exact]

med <- dcast(girf_bands[, .(variant, horizon, median)],
             horizon ~ variant, value.var = "median")
setnames(med, VARIANT_LEVELS, c("exact", "pf_macro", "pf_full"))
med[, bias_pf_macro_pct := fifelse(abs(exact) > 1e-8,
                                   100 * (pf_macro - exact) / abs(exact), NA_real_)]
med[, bias_pf_full_pct  := fifelse(abs(exact) > 1e-8,
                                   100 * (pf_full  - exact) / abs(exact), NA_real_)]

fwrite(lev, file.path(pf_out_dir, "tables", "shocked_level_bias_by_horizon.csv"))
fwrite(med, file.path(pf_out_dir, "tables", "girf_bias_by_horizon.csv"))

peak_row <- med[which.max(abs(exact))]

tex_path <- file.path(pf_out_dir, "tables", "perfect_foresight_gulf_war_bias.tex")
con <- file(tex_path, open = "wt")
cat("\\begin{table}[!htbp]\n\\centering\n", file = con)
cat("\\caption{The cost of perfect foresight at the Gulf-War state (1990:Q3)}\n", file = con)
cat("\\label{tab:perfect_foresight_gulf_war_bias}\n", file = con)
cat("\\begin{tabular}{llrrrrr}\n\\toprule\n", file = con)
cat(" & & Exact & PF macro & PF plug-in & Bias macro (\\%) & Bias plug-in (\\%) \\\\\n\\midrule\n", file = con)
for (hh in c(4L, 8L, 12L)) {
  r <- lev[horizon == hh]
  cat(sprintf("Shocked PD level (pp) & $h=%d$ & %s & %s & %s & %s & %s \\\\\n",
              hh, fmt(r$exact), fmt(r$pf_macro), fmt(r$pf_full),
              fmt(r$bias_pf_macro_pct, 1), fmt(r$bias_pf_full_pct, 1)),
      file = con)
}
cat("\\midrule\n", file = con)
{
  r <- peak_row
  cat(sprintf("Peak $\\Delta$PD (pp) & $h=%d$ & %s & %s & %s & %s & %s \\\\\n",
              as.integer(r$horizon),
              fmt(r$exact), fmt(r$pf_macro), fmt(r$pf_full),
              fmt(r$bias_pf_macro_pct, 1), fmt(r$bias_pf_full_pct, 1)),
      file = con)
}
cat("\\bottomrule\n\\end{tabular}\n", file = con)
cat("\\begin{tablenotes}\n\\small\n", file = con)
cat(sprintf("\\item Notes: Gulf-War scenario (own historical GPR shock, %.2f s.d.), conditional on the credit state preceding 1990:Q3. Posterior medians on identical posterior draws; the variants differ only in the predictive variance entering the closed forms. ``Exact'' is Proposition~3; ``PF macro'' treats the macro-financial path as known but retains the satellite-error variance (the closest counterpart to the perfect-foresight convention of supervisory stress tests); ``PF plug-in'' evaluates the Merton--Vasicek map at the conditional-mean path. The shocked PD level is the absolute PD along the scenario; biases are relative to the exact variant, with negative level biases indicating understatement of the scenario PD.\n", gw_size), file = con)
cat("\\end{tablenotes}\n\\end{table}\n", file = con)
close(con)

# -----------------------------------------------------------------------------
# 6. Figure: shocked PD level path, exact vs the two perfect-foresight variants
# -----------------------------------------------------------------------------

cap_txt <- sprintf(
  "PD level along the Gulf-War scenario (own historical GPR shock, %.2f s.d.; credit state preceding 1990:Q3): exact closed form versus perfect foresight (macro path known; full plug-in). Shaded area: 68%% band of the exact shocked path.",
  gw_size)

p_shock_lev <- ggplot(shocked_levels_bands,
                      aes(x = horizon, y = median,
                          color = variant, linetype = variant)) +
  geom_ribbon(data = shocked_levels_bands[variant == VARIANT_LEVELS[1]],
              aes(x = horizon, ymin = lower68, ymax = upper68),
              inherit.aes = FALSE, fill = "#AB4A7D", alpha = 0.14) +
  geom_line(linewidth = 0.95) +
  scale_color_manual(values = VARIANT_COLORS) +
  scale_linetype_manual(values = VARIANT_LTY) +
  scale_x_continuous(breaks = seq(0, 12, by = 1), minor_breaks = NULL) +
  labs(
    x = "Horizon (quarters)",
    y = "Shocked PD path (percentage points)",
    caption = cap_txt
  ) +
  theme_paper()

save_fig(p_shock_lev, "fig_perfect_foresight_gulf_war_shocked_levels",
         width = 7.2, height = 4.4)

# -----------------------------------------------------------------------------
# 7. Save and recap
# -----------------------------------------------------------------------------

saveRDS(
  list(info_set = PF_INFO_SET, p_ttc = P_TTC, rho = RHO, n_draws = M,
       gulf_quarter = GULF_QUARTER, gulf_shock_size = gw_size,
       gulf_delinquency_obs = if (length(d_obs_gw)) d_obs_gw[1] else NA_real_,
       shocked_level_bands_pp = shocked_levels_bands, girf_bands_pp = girf_bands,
       shocked_level_bias = lev, girf_bias = med),
  file.path(pf_out_dir, "perfect_foresight_gulf_war_results.rds")
)

readme <- c(
  "# DRALACBN -- perfect foresight at the Gulf-War state (1990:Q3)",
  "",
  "Combines the perfect-foresight bias (08_*) with the state-dependence",
  "conditioning (10_*): the exact closed forms vs the two perfect-foresight",
  "variants, evaluated along the Gulf-War scenario (own historical GPR shock,",
  "conditional on the stressed credit state preceding 1990:Q3 -- the largest",
  "responses in the sample). Same posterior draws in all variants; only the",
  "predictive variance differs:",
  "  exact     Proposition 3 (full conditional variance)",
  "  pf_macro  macro path known, satellite-error variance retained",
  "  pf_full   plug-in: PD evaluated at the conditional-mean path",
  "",
  "Generated by code/extra/dralacbn/11_dralacbn_perfect_foresight_gulf_war.R (standalone).",
  "",
  "Contents:",
  "  figures/fig_perfect_foresight_gulf_war_shocked_levels  shocked PD path: exact vs PF",
  "  tables/perfect_foresight_gulf_war_bias.tex|csv         shocked-level & peak biases",
  "  irf/*.csv                                              bands by variant (pp)",
  "  perfect_foresight_gulf_war_results.rds                 all objects"
)
writeLines(readme, file.path(pf_out_dir, "README.md"))

cat("\n--- Perfect foresight at the Gulf-War state (posterior medians) ---------\n")
cat("Shocked PD level (pp) h=12: exact = ", fmt(lev[horizon == 12, exact]),
    " | PF plug-in = ", fmt(lev[horizon == 12, pf_full]),
    " (bias ", fmt(lev[horizon == 12, bias_pf_full_pct], 1), "%)\n", sep = "")
cat("Peak dPD (pp) h=", as.integer(peak_row$horizon),
    ": exact = ", fmt(peak_row$exact),
    " | PF plug-in = ", fmt(peak_row$pf_full),
    " (bias ", fmt(peak_row$bias_pf_full_pct, 1), "%)\n", sep = "")
cat("Outputs written to: ", pf_out_dir, "\n", sep = "")
cat("--------------------------------------------------------------------------\n")

invisible(TRUE)
