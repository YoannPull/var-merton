# =============================================================================
#  code/dralacbn/08_dralacbn_perfect_foresight.R
#
#  The cost of perfect foresight.
#
#  Supervisory stress tests treat the macroeconomic scenario as a KNOWN path
#  (EBA 2025 EU-wide stress test, Box 1 and par. 130) and project risk
#  parameters along it. In the framework of the paper, this amounts to
#  evaluating the Merton-Vasicek map at the conditional-mean path of the
#  systematic factor, i.e. setting the conditional variances of Proposition 1
#  to zero in the closed forms. This script quantifies the resulting bias on
#  the DRALACBN baseline specification, on the SAME posterior draws as the
#  main results (estimation uncertainty is common to all variants; only
#  predictive uncertainty differs).
#
#  Three variants:
#    exact     full closed form (Prop. 3): s2 = macro forecast variance + eta
#    pf_macro  macro path known, satellite noise retained: s2 = sigma_eta^2
#              (closest to the EBA perfect-foresight convention)
#    pf_full   full plug-in: s2 = 0 (PD evaluated at the mean path)
#
#  Predictions (sign): pi is convex in the relevant region (p small), so
#  perfect foresight UNDERSTATES the expected PD level (Jensen); the variance
#  smooths the map, so perfect foresight OVERSTATES the marginal response.
#  Both biases grow with the horizon and with the size of the shock.
#
#  STANDALONE: modifies nothing in the existing pipeline; not in run_all.R
#  unless added. Requires the main DRALACBN application (02) to have run.
#  Run from the repository root:
#      source("code/dralacbn/08_dralacbn_perfect_foresight.R")
#
#  Outputs: output/applications/dralacbn/perfect_foresight/
# =============================================================================

source("code/00_setup.R")
source("code/shared/_helpers_var_merton.R")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

# -----------------------------------------------------------------------------
# 0. Parameters and output tree
# -----------------------------------------------------------------------------

PF_INFO_SET  <- if (exists("DELINQUENCY_BASELINE_INFO_SET")) DELINQUENCY_BASELINE_INFO_SET else "real_side"
p_lags       <- if (exists("P_LAGS")) P_LAGS else 2
nrep         <- if (exists("NREP")) NREP else 30000
M_target     <- if (exists("M_TARGET")) M_TARGET else 10000L
H            <- if (exists("HORIZON")) HORIZON else 12
seed         <- if (exists("SEED_ROBUSTNESS")) SEED_ROBUSTNESS else 123
impulse_name <- if (exists("IMPULSE_NAME")) IMPULSE_NAME else "log_GPRD"

pf_out_dir <- file.path("output", "applications", "dralacbn", "perfect_foresight")
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
      strip.text = element_text(face = "bold", color = "#1A1A1A", size = base_size),
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

cat("\n========== DRALACBN: the cost of perfect foresight ==========\n")

# -----------------------------------------------------------------------------
# 1. Inputs: Merton parameters, kernel, main-specification Bayesian satellite
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
n_sat <- nrow(bsat$beta_term_draws)
cat("Bayesian satellite loaded (", n_sat, " draws, main specification).\n", sep = "")

# Satellite-error variance per BVAR draw (engine pairing convention).
sigma2_paired <- bsat$sigma2_draws[vapply(seq_len(M), pair_satellite_draw,
                                          integer(1), n_sat = n_sat)]

# -----------------------------------------------------------------------------
# 2. Conditional moments (computed once) and the three variance variants
# -----------------------------------------------------------------------------

cat("Computing conditional factor moments ...\n")
mom <- compute_factor_moment_draws_bayes(kern, bsat, impulse_ix, p_lags, H)

zero_mat <- matrix(0.0, nrow = H + 1, ncol = M)
eta_mat  <- matrix(rep(sigma2_paired, each = H + 1), nrow = H + 1, ncol = M)

variants <- list(
  exact    = list(s2 = mom$s2_draws, s2d = mom$s2_delta_draws,
                  label = VARIANT_LEVELS[1]),
  pf_macro = list(s2 = eta_mat, s2d = eta_mat,
                  label = VARIANT_LEVELS[2]),
  pf_full  = list(s2 = zero_mat, s2d = zero_mat,
                  label = VARIANT_LEVELS[3])
)

shock_2001 <- compute_structural_e1_median(
  kernel = kern, dates_raw = as.Date(DT_raw$Date),
  target_quarter = "2001Q3", impulse_idx = impulse_ix
)

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

# -----------------------------------------------------------------------------
# 3. Baseline PD levels and PD GIRFs under each variant (in pp)
# -----------------------------------------------------------------------------

pd_level_draws <- function(mu_mat, s2_mat) {
  qpi <- qnorm(P_TTC)
  denom <- sqrt((1 - RHO) + RHO * pmax(s2_mat, 0))
  100 * pnorm((qpi - sqrt(RHO) * mu_mat) / denom)   # percentage points
}

levels_bands <- rbindlist(lapply(names(variants), function(v) {
  bands_from_draws(pd_level_draws(mom$mu_draws, variants[[v]]$s2),
                   list(variant = variants[[v]]$label))
}))

cat("Computing PD GIRFs under the three variants ...\n")
psiZ_1sd  <- inject_girf_into_Z_bayes(kern$Psi_draws, kern$Sigma_kept, bsat,
                                      impulse_ix, shock_scale = 1)
psiZ_2001 <- inject_girf_into_Z_bayes(kern$Psi_draws, kern$Sigma_kept, bsat,
                                      impulse_ix, shock_scale = shock_2001$shock_scale)

girf_bands <- rbindlist(lapply(names(variants), function(v) {
  vv <- variants[[v]]
  rbindlist(lapply(list(list(psiZ_1sd, "one_sd"), list(psiZ_2001, "2001Q3")),
                   function(sh) {
    pd <- compute_pd_girf(sh[[1]]$draws, mom$mu_draws, vv$s2, vv$s2d,
                          p = P_TTC, rho = RHO)
    b <- bands_from_draws(100 * pd$draws,
                          list(variant = vv$label, shock = sh[[2]]))
    b
  }))
}))

girf_bands[, variant := factor(variant, levels = VARIANT_LEVELS)]
levels_bands[, variant := factor(variant, levels = VARIANT_LEVELS)]

fwrite(girf_bands, file.path(pf_out_dir, "irf", "pd_girf_bands_by_variant_pp.csv"))
fwrite(levels_bands, file.path(pf_out_dir, "irf", "pd_level_bands_by_variant_pp.csv"))

# -----------------------------------------------------------------------------
# 3b. Shocked PD level trajectories (2001Q3 scenario), by variant
#
#     Absolute PD path ALONG the shocked scenario, i.e. the level a supervisor
#     reads off as "the PD under the adverse scenario", NOT the response:
#       pd_shock = Phi( (qpi - sqrt(rho) * (mu_base + psiZ)) / s_shock ),
#     with s_shock built from the shocked-scenario predictive variance of each
#     variant (s2_delta) -- the same convention as pd_shock in compute_pd_girf.
#     Comparing exact vs the two perfect-foresight variants shows the LEVEL
#     distortion of treating the macro path as known. 2001Q3 shock only.
# -----------------------------------------------------------------------------

pd_shock_level_draws <- function(mu_mat, psiZ_mat, s2d_mat) {
  qpi <- qnorm(P_TTC)
  mu_shock <- mu_mat + psiZ_mat
  denom <- sqrt((1 - RHO) + RHO * pmax(s2d_mat, 0))
  100 * pnorm((qpi - sqrt(RHO) * mu_shock) / denom)   # percentage points
}

shocked_levels_bands <- rbindlist(lapply(names(variants), function(v) {
  vv <- variants[[v]]
  bands_from_draws(
    pd_shock_level_draws(mom$mu_draws, psiZ_2001$draws, vv$s2d),
    list(variant = vv$label)
  )
}))
shocked_levels_bands[, variant := factor(variant, levels = VARIANT_LEVELS)]

fwrite(shocked_levels_bands,
       file.path(pf_out_dir, "irf",
                 "pd_shocked_level_bands_by_variant_2001Q3_pp.csv"))

# Shocked-level bias table (2001Q3 scenario), consistent with the shocked-level
# figure: absolute PD along the scenario, not the baseline level. Reported
# alongside the peak DeltaPD bias of the same 2001Q3 shock.
lev_shock <- dcast(shocked_levels_bands[, .(variant, horizon, median)],
                   horizon ~ variant, value.var = "median")
setnames(lev_shock, VARIANT_LEVELS, c("exact", "pf_macro", "pf_full"))
lev_shock[, bias_pf_macro_pct := 100 * (pf_macro - exact) / exact]
lev_shock[, bias_pf_full_pct  := 100 * (pf_full  - exact) / exact]
fwrite(lev_shock,
       file.path(pf_out_dir, "tables", "shocked_level_bias_by_horizon.csv"))

g2001 <- dcast(girf_bands[shock == "2001Q3", .(variant, horizon, median)],
               horizon ~ variant, value.var = "median")
setnames(g2001, VARIANT_LEVELS, c("exact", "pf_macro", "pf_full"))
g2001[, bias_pf_macro_pct := fifelse(abs(exact) > 1e-8,
                                     100 * (pf_macro - exact) / abs(exact), NA_real_)]
g2001[, bias_pf_full_pct  := fifelse(abs(exact) > 1e-8,
                                     100 * (pf_full  - exact) / abs(exact), NA_real_)]
peak2001 <- g2001[which.max(abs(exact))]

shock_tex <- file.path(pf_out_dir, "tables", "perfect_foresight_shocked_bias.tex")
con <- file(shock_tex, open = "wt")
cat("\\begin{table}[!htbp]\n\\centering\n", file = con)
cat("\\caption{The cost of perfect foresight along the September 11 scenario (2001:Q3)}\n", file = con)
cat("\\label{tab:perfect_foresight_shocked_bias}\n", file = con)
cat("\\begin{tabular}{llrrrrr}\n\\toprule\n", file = con)
cat(" & & Exact & PF macro & PF plug-in & Bias macro (\\%) & Bias plug-in (\\%) \\\\\n\\midrule\n", file = con)
for (hh in c(4L, 8L, 12L)) {
  r <- lev_shock[horizon == hh]
  cat(sprintf("Shocked PD level (pp) & $h=%d$ & %s & %s & %s & %s & %s \\\\\n",
              hh, fmt(r$exact), fmt(r$pf_macro), fmt(r$pf_full),
              fmt(r$bias_pf_macro_pct, 1), fmt(r$bias_pf_full_pct, 1)),
      file = con)
}
cat("\\midrule\n", file = con)
cat(sprintf("Peak $\\Delta$PD (pp) & $h=%d$ & %s & %s & %s & %s & %s \\\\\n",
            as.integer(peak2001$horizon),
            fmt(peak2001$exact), fmt(peak2001$pf_macro), fmt(peak2001$pf_full),
            fmt(peak2001$bias_pf_macro_pct, 1), fmt(peak2001$bias_pf_full_pct, 1)),
    file = con)
cat("\\bottomrule\n\\end{tabular}\n", file = con)
cat("\\begin{tablenotes}\n\\small\n", file = con)
cat("\\item Notes: September 11 scenario (2001Q3-sized GPR shock). Posterior medians on identical posterior draws; the variants differ only in the predictive variance entering the closed forms. ``Exact'' is Proposition~3; ``PF macro'' treats the macro-financial path as known but retains the satellite-error variance (the closest counterpart to the perfect-foresight convention of supervisory stress tests); ``PF plug-in'' evaluates the Merton--Vasicek map at the conditional-mean path. The shocked PD level is the absolute PD along the scenario; biases are relative to the exact variant, with negative level biases indicating understatement of the scenario PD.\n", file = con)
cat("\\end{tablenotes}\n\\end{table}\n", file = con)
close(con)

# -----------------------------------------------------------------------------
# 4. Bias summaries
# -----------------------------------------------------------------------------

med <- dcast(girf_bands[, .(variant, shock, horizon, median)],
             shock + horizon ~ variant, value.var = "median")
setnames(med, VARIANT_LEVELS, c("exact", "pf_macro", "pf_full"))
med[, bias_pf_macro_pct := fifelse(abs(exact) > 1e-8,
                                   100 * (pf_macro - exact) / abs(exact), NA_real_)]
med[, bias_pf_full_pct := fifelse(abs(exact) > 1e-8,
                                  100 * (pf_full - exact) / abs(exact), NA_real_)]

lev <- dcast(levels_bands[, .(variant, horizon, median)],
             horizon ~ variant, value.var = "median")
setnames(lev, VARIANT_LEVELS, c("exact", "pf_macro", "pf_full"))
lev[, bias_pf_macro_pct := 100 * (pf_macro - exact) / exact]
lev[, bias_pf_full_pct := 100 * (pf_full - exact) / exact]

fwrite(med, file.path(pf_out_dir, "tables", "girf_bias_by_horizon.csv"))
fwrite(lev, file.path(pf_out_dir, "tables", "level_bias_by_horizon.csv"))

peak_rows <- med[, .SD[which.max(abs(exact))], by = shock]

tex_path <- file.path(pf_out_dir, "tables", "perfect_foresight_bias.tex")
con <- file(tex_path, open = "wt")
cat("\\begin{table}[!htbp]\n\\centering\n", file = con)
cat("\\caption{The cost of perfect foresight: PD levels and responses}\n", file = con)
cat("\\label{tab:perfect_foresight_bias}\n", file = con)
cat("\\begin{tabular}{llrrrrr}\n\\toprule\n", file = con)
cat(" & & Exact & PF macro & PF plug-in & Bias macro (\\%) & Bias plug-in (\\%) \\\\\n\\midrule\n", file = con)
for (hh in c(4L, 8L, 12L)) {
  r <- lev[horizon == hh]
  cat(sprintf("PD level (pp) & $h=%d$ & %s & %s & %s & %s & %s \\\\\n",
              hh, fmt(r$exact), fmt(r$pf_macro), fmt(r$pf_full),
              fmt(r$bias_pf_macro_pct, 1), fmt(r$bias_pf_full_pct, 1)),
      file = con)
}
cat("\\midrule\n", file = con)
for (i in seq_len(nrow(peak_rows))) {
  r <- peak_rows[i]
  sh_lab <- if (r$shock == "one_sd") "One s.d." else "2001Q3"
  cat(sprintf("Peak $\\Delta$PD (pp), %s & $h=%d$ & %s & %s & %s & %s & %s \\\\\n",
              sh_lab, as.integer(r$horizon),
              fmt(r$exact), fmt(r$pf_macro), fmt(r$pf_full),
              fmt(r$bias_pf_macro_pct, 1), fmt(r$bias_pf_full_pct, 1)),
      file = con)
}
cat("\\bottomrule\n\\end{tabular}\n", file = con)
cat("\\begin{tablenotes}\n\\small\n", file = con)
cat("\\item Notes: Posterior medians on identical posterior draws; the variants differ only in the predictive variance entering the closed forms. ``Exact'' is Proposition~3; ``PF macro'' treats the macro-financial path as known but retains the satellite-error variance, the closest counterpart to the perfect-foresight convention of supervisory stress tests; ``PF plug-in'' evaluates the Merton--Vasicek map at the conditional-mean path. Biases are relative to the exact variant; negative level biases indicate understatement of the expected PD path.\n", file = con)
cat("\\end{tablenotes}\n\\end{table}\n", file = con)
close(con)

# -----------------------------------------------------------------------------
# 5. Figures
# -----------------------------------------------------------------------------

exact_band <- girf_bands[variant == VARIANT_LEVELS[1]]
girf_plot <- copy(girf_bands)
girf_plot[, shock_label := fifelse(shock == "one_sd",
                                   "One-s.d. GPR shock", "2001Q3-sized GPR shock")]
girf_plot[, shock_label := factor(shock_label,
                                  levels = c("One-s.d. GPR shock",
                                             "2001Q3-sized GPR shock"))]
exact_band2 <- girf_plot[variant == VARIANT_LEVELS[1]]

p_girf <- ggplot(girf_plot, aes(x = horizon, y = median,
                                color = variant, linetype = variant)) +
  geom_ribbon(data = exact_band2,
              aes(x = horizon, ymin = lower68, ymax = upper68),
              inherit.aes = FALSE, fill = "#AB4A7D", alpha = 0.14) +
  geom_hline(yintercept = 0, color = "grey45", linewidth = 0.35, linetype = "dashed") +
  geom_line(linewidth = 0.95) +
  facet_wrap(~ shock_label, nrow = 1, scales = "free_y") +
  scale_color_manual(values = VARIANT_COLORS) +
  scale_linetype_manual(values = VARIANT_LTY) +
  scale_x_continuous(breaks = seq(0, 12, by = 2), minor_breaks = NULL) +
  labs(
    x = "Horizon (quarters)",
    y = expression(Delta ~ "PD (percentage points)"),
    caption = "Posterior medians on identical draws; shaded area: 68% band of the exact response. The variants differ only in the predictive variance entering the closed forms."
  ) +
  theme_paper()

save_fig(p_girf, "fig_perfect_foresight_girf", width = 9.2, height = 4.0)

p_lev <- ggplot(levels_bands, aes(x = horizon, y = median,
                                  color = variant, linetype = variant)) +
  geom_ribbon(data = levels_bands[variant == VARIANT_LEVELS[1]],
              aes(x = horizon, ymin = lower68, ymax = upper68),
              inherit.aes = FALSE, fill = "#AB4A7D", alpha = 0.14) +
  geom_line(linewidth = 0.95) +
  scale_color_manual(values = VARIANT_COLORS) +
  scale_linetype_manual(values = VARIANT_LTY) +
  scale_x_continuous(breaks = seq(0, 12, by = 1), minor_breaks = NULL) +
  labs(
    x = "Horizon (quarters)",
    y = "Baseline PD path (percentage points)",
    caption = "Expected PD path (exact) versus PD evaluated along the conditional-mean path (perfect foresight). Shaded area: 68% band of the exact path."
  ) +
  theme_paper()

save_fig(p_lev, "fig_perfect_foresight_levels", width = 7.2, height = 4.4)

# Shocked PD level path along the 2001Q3 scenario: exact closed form versus the
# two perfect-foresight variants. Single panel (2001Q3 shock only).
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
    caption = "PD level along the 2001Q3-sized GPR scenario: exact closed form versus perfect foresight (macro path known; full plug-in). Shaded area: 68% band of the exact shocked path."
  ) +
  theme_paper()

save_fig(p_shock_lev, "fig_perfect_foresight_shocked_levels",
         width = 7.2, height = 4.4)

bias_plot <- melt(med[horizon >= 1],
                  id.vars = c("shock", "horizon"),
                  measure.vars = c("bias_pf_macro_pct", "bias_pf_full_pct"),
                  variable.name = "variant", value.name = "bias_pct")
bias_plot[, variant := fifelse(variant == "bias_pf_macro_pct",
                               VARIANT_LEVELS[2], VARIANT_LEVELS[3])]
bias_plot[, shock_label := fifelse(shock == "one_sd",
                                   "One-s.d. GPR shock", "2001Q3-sized GPR shock")]
p_bias <- ggplot(bias_plot, aes(x = horizon, y = bias_pct,
                                color = variant, linetype = variant)) +
  geom_hline(yintercept = 0, color = "grey45", linewidth = 0.35, linetype = "dashed") +
  geom_line(linewidth = 0.95) +
  geom_point(size = 1.4) +
  facet_wrap(~ shock_label, nrow = 1) +
  scale_color_manual(values = VARIANT_COLORS[2:3]) +
  scale_linetype_manual(values = VARIANT_LTY[2:3]) +
  scale_x_continuous(breaks = seq(0, 12, by = 2), minor_breaks = NULL) +
  labs(
    x = "Horizon (quarters)",
    y = expression("Bias of " * Delta * "PD relative to the exact response (%)"),
    caption = "Median PF response minus median exact response, relative to the absolute exact response."
  ) +
  theme_paper()

save_fig(p_bias, "fig_perfect_foresight_bias", width = 9.2, height = 4.0)

# -----------------------------------------------------------------------------
# 6. Save and recap
# -----------------------------------------------------------------------------

saveRDS(
  list(info_set = PF_INFO_SET, p_ttc = P_TTC, rho = RHO, n_draws = M,
       girf_bands_pp = girf_bands, level_bands_pp = levels_bands,
       shocked_level_bands_pp = shocked_levels_bands,
       girf_bias = med, level_bias = lev, shocked_level_bias = lev_shock,
       shock_2001 = shock_2001$target),
  file.path(pf_out_dir, "perfect_foresight_results.rds")
)

readme <- c(
  "# DRALACBN -- the cost of perfect foresight",
  "",
  "Quantifies the bias of the perfect-foresight convention of supervisory",
  "stress tests (macro path treated as known; EBA 2025, Box 1 / par. 130)",
  "relative to the exact closed-form responses of the paper. Same posterior",
  "draws in all variants; only the predictive variance differs:",
  "  exact     Proposition 3 (full conditional variance)",
  "  pf_macro  macro path known, satellite-error variance retained",
  "  pf_full   plug-in: PD evaluated at the conditional-mean path",
  "",
  "Generated by code/dralacbn/08_dralacbn_perfect_foresight.R (standalone).",
  "",
  "Contents:",
  "  figures/fig_perfect_foresight_girf            Delta PD: exact vs PF, both shocks",
  "  figures/fig_perfect_foresight_levels          baseline PD path: exact vs PF",
  "  figures/fig_perfect_foresight_shocked_levels  shocked PD path (2001Q3): exact vs PF",
  "  figures/fig_perfect_foresight_bias            % bias by horizon",
  "  tables/perfect_foresight_bias.tex|csv         baseline levels and peaks with % biases",
  "  tables/perfect_foresight_shocked_bias.tex      shocked-level (2001Q3) biases by horizon",
  "  irf/*.csv                               bands by variant (pp)",
  "  perfect_foresight_results.rds           all objects"
)
writeLines(readme, file.path(pf_out_dir, "README.md"))

cat("\n--- Perfect-foresight bias (posterior medians) --------------------------\n")
cat("PD level at h=12 (pp): exact = ", fmt(lev[horizon == 12, exact]),
    " | PF plug-in = ", fmt(lev[horizon == 12, pf_full]),
    " (bias ", fmt(lev[horizon == 12, bias_pf_full_pct], 1), "%)\n", sep = "")
for (i in seq_len(nrow(peak_rows))) {
  r <- peak_rows[i]
  cat("Peak dPD ", r$shock, " (pp): exact = ", fmt(r$exact),
      " | PF plug-in = ", fmt(r$pf_full),
      " (bias ", fmt(r$bias_pf_full_pct, 1), "%)\n", sep = "")
}
cat("Outputs written to: ", pf_out_dir, "\n", sep = "")
cat("--------------------------------------------------------------------------\n")

invisible(TRUE)
