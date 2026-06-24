# =============================================================================
#  code/extra/dralacbn/05_dralacbn_short_sample.R
#
#  DRALACBN application — asymmetric-sample (short credit history) exercise.
#
#  The paper's modular architecture claims that the macro-financial block
#  (long history) and the credit-risk block (short history) can be estimated
#  on mismatched samples. This script demonstrates the claim: the BVAR is kept
#  on its full 1986-2024 sample, while the DEFAULT sample used to reconstruct
#  Z and estimate the satellite is progressively truncated:
#
#      1985-2024 (full), 1995-2024, 2005-2024, 2015-2024.
#
#  Practitioner setting: within each window, the Merton parameters (p_ttc,
#  rho) and the satellite are re-estimated from scratch on the short default
#  sample only. The satellite DESIGN is fixed at the full-sample best-BIC
#  specification, so the exercise isolates the effect of the sample length
#  (not of model selection).
#
#  Expected message: posterior-median PD responses remain stable while
#  credible bands widen as the default sample shrinks -- the macro-credit
#  bridge stays estimable with as little as ~40 quarters of default data.
#
#  STANDALONE: modifies nothing in the existing pipeline; not in run_all.R.
#  Run from the repository root AFTER the main DRALACBN application:
#      source("code/extra/dralacbn/05_dralacbn_short_sample.R")
#
#  Outputs: output/applications/dralacbn/short_sample/
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

SS_INFO_SET    <- if (exists("DELINQUENCY_BASELINE_INFO_SET")) DELINQUENCY_BASELINE_INFO_SET else "real_side"
SS_START_YEARS <- c(1985, 1995, 2005, 2015)   # first year of each default window

p_lags       <- if (exists("P_LAGS")) P_LAGS else 2
nrep         <- if (exists("NREP")) NREP else 30000
M_target     <- if (exists("M_TARGET")) M_TARGET else 10000L
H            <- if (exists("HORIZON")) HORIZON else 12
seed         <- if (exists("SEED_ROBUSTNESS")) SEED_ROBUSTNESS else 123
seed_sat     <- if (exists("SEED_ZFACTOR")) SEED_ZFACTOR else 12345
impulse_name <- if (exists("IMPULSE_NAME")) IMPULSE_NAME else "log_GPRD"
NB_LAGS_SAT  <- if (exists("ALT_SPEC_NB_LAGS_SAT")) ALT_SPEC_NB_LAGS_SAT else 4

ss_out_dir <- file.path("output", "applications", "dralacbn", "short_sample")
for (sub in c("", "tables", "figures", "irf")) {
  dir.create(file.path(ss_out_dir, sub), recursive = TRUE, showWarnings = FALSE)
}

COL_WINDOWS <- c("#6F1732", "#AB4A7D", "#D78FB4", "#6B6B6B")

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
  ggsave(file.path(ss_out_dir, "figures", paste0(stem, ".png")),
         plot, width = width, height = height, dpi = dpi)
  tryCatch(
    ggsave(file.path(ss_out_dir, "figures", paste0(stem, ".pdf")),
           plot, width = width, height = height, device = grDevices::pdf),
    error = function(e) tryCatch(
      ggsave(file.path(ss_out_dir, "figures", paste0(stem, ".pdf")),
             plot, width = width, height = height),
      error = function(e2) message("PDF export failed for ", stem)
    )
  )
}

cat("\n========== DRALACBN short-sample (asymmetric-sample) exercise ==========\n")

# -----------------------------------------------------------------------------
# 1. Load default proxy, data, BVAR kernel and full-sample satellite design
# -----------------------------------------------------------------------------

delinquency_path <- file.path(
  if (exists("DIR_RAW")) DIR_RAW else "data/raw", "default", "DRALACBN.csv"
)
if (!file.exists(delinquency_path)) {
  stop("Missing delinquency series: ", delinquency_path, call. = FALSE)
}

delinq_raw <- fread(delinquency_path)
date_col <- intersect(c("DATE", "Date", "date", "observation_date"), names(delinq_raw))[1]
val_col  <- intersect(c("DRALACBN", "value", "VALUE", "Value"), names(delinq_raw))[1]
if (is.na(date_col)) date_col <- names(delinq_raw)[1]
if (is.na(val_col))  val_col  <- setdiff(names(delinq_raw), date_col)[1]

delinq <- data.table(
  Date = as.Date(delinq_raw[[date_col]]),
  value = as.numeric(delinq_raw[[val_col]])
)
delinq <- delinq[is.finite(value)]
setorder(delinq, Date)
if (exists("DATA_END_DATE")) delinq <- delinq[Date <= DATA_END_DATE]

DT_raw <- safe_fread(if (exists("PATH_DATA_VAR")) PATH_DATA_VAR else "data/processed/data_var_for_model.csv")
DT_raw[, Date := as.Date(Date)]

kern <- get_or_estimate_kernel(
  spec_name = SS_INFO_SET,
  vars_extra = VAR_SPECS[[SS_INFO_SET]],
  DT_raw = DT_raw,
  p_lags = p_lags, nrep = nrep, H = H, seed = seed, M_target = M_target,
  kernel_dir = if (exists("DIR_VAR_KERNELS")) DIR_VAR_KERNELS else file.path("output", "shared", "var_kernels"),
  impulse_name = impulse_name
)

variables  <- kern$variables
impulse_ix <- match(impulse_name, variables)
if (is.na(impulse_ix)) stop("Impulse variable not found in VAR kernel.", call. = FALSE)
M <- dim(kern$Psi_draws)[4]

best_sat_path <- file.path("output", "applications", "dralacbn",
                           "main_specification", "satellite", "best_satellite.rds")
if (!file.exists(best_sat_path)) {
  stop("Missing ", best_sat_path,
       "\nRun the main DRALACBN application first.", call. = FALSE)
}
best_sat <- readRDS(best_sat_path)
sat_formula <- formula(best_sat)
sat_terms <- setdiff(names(coef(best_sat)), "(Intercept)")

cat("VAR sample (unchanged): full macro-financial history, ", M,
    " posterior draws.\n", sep = "")
cat("Satellite design (fixed, full-sample best-BIC): ",
    paste(sat_terms, collapse = " + "), "\n", sep = "")

vars_sat <- setdiff(variables, impulse_name)

shock_2001 <- compute_structural_e1_median(
  kernel = kern,
  dates_raw = as.Date(DT_raw$Date),
  target_quarter = "2001Q3",
  impulse_idx = impulse_ix
)

# -----------------------------------------------------------------------------
# 2. One run per default window
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

run_window <- function(start_year) {
  d_w <- delinq[Date >= as.Date(paste0(start_year, "-01-01"))]
  Td <- nrow(d_w)

  z_out <- f_Z_estimation(d_w$value / 100)
  if (!is.finite(z_out$rho)) {
    stop("Unit-variance calibration of rho failed for window ", start_year,
         call. = FALSE)
  }

  Z_DT_w <- data.table(
    Date = as.Date(format(d_w$Date + 62, "%Y-%m-01")),
    Z = z_out$Z
  )

  sat_df_w <- prepare_lagged_satellite_data(
    z_input = Z_DT_w, DT_raw = DT_raw, vars_sat = vars_sat, nb_lags = NB_LAGS_SAT
  )

  fit_w <- lm(sat_formula, data = sat_df_w)
  r2_w <- summary(fit_w)$r.squared

  bsat_w <- make_bayesian_satellite(
    sat_model = fit_w,
    model_df = sat_df_w,
    variables = variables,
    n_draws = M,
    seed = seed_sat
  )

  cat(sprintf("  Window %d-2024: Td = %3d, n_sat = %3d, p_ttc = %.4f, rho = %.4f, R2 = %.3f\n",
              start_year, Td, stats::nobs(fit_w), z_out$p_ttc, z_out$rho, r2_w))

  mom <- compute_factor_moment_draws_bayes(kern, bsat_w, impulse_ix, p_lags, H)

  out_shock <- function(scale, lab) {
    psiZ <- inject_girf_into_Z_bayes(kern$Psi_draws, kern$Sigma_kept, bsat_w,
                                     impulse_ix, shock_scale = scale)
    pd <- compute_pd_girf(psiZ$draws, mom$mu_draws, mom$s2_draws,
                          mom$s2_delta_draws, p = z_out$p_ttc, rho = z_out$rho)
    list(
      psiZ = bands_from_draws(psiZ$draws, list(window = start_year, shock = lab)),
      pd = bands_from_draws(pd$draws, list(window = start_year, shock = lab))
    )
  }

  s1 <- out_shock(1, "one_sd")
  s2 <- out_shock(shock_2001$shock_scale, "2001Q3")

  list(
    start_year = start_year, Td = Td, n_sat = stats::nobs(fit_w),
    p_ttc = z_out$p_ttc, rho = z_out$rho, r2 = r2_w,
    psiZ = rbind(s1$psiZ, s2$psiZ),
    pd = rbind(s1$pd, s2$pd)
  )
}

cat("\nRunning windows (VAR block identical in all of them):\n")
runs <- lapply(SS_START_YEARS, run_window)
names(runs) <- as.character(SS_START_YEARS)

pd_all <- rbindlist(lapply(runs, function(r) r$pd))
zz_all <- rbindlist(lapply(runs, function(r) r$psiZ))

# Window labels with sample size, used everywhere for clarity.
meta <- rbindlist(lapply(runs, function(r) {
  data.table(window = r$start_year, Td = r$Td, n_sat = r$n_sat,
             p_ttc = r$p_ttc, rho = r$rho, r2 = r$r2)
}))
meta[, label := sprintf("%d-2024  (T[d] == %d)", window, Td)]
meta[, label_plain := sprintf("%d-2024 (Td = %d)", window, Td)]
win_levels <- meta$label_plain

add_labels <- function(dt) {
  dt <- merge(dt, meta[, .(window, label_plain)], by = "window", sort = FALSE)
  dt[, window_label := factor(label_plain, levels = win_levels)]
  dt
}
pd_all <- add_labels(pd_all)
zz_all <- add_labels(zz_all)

# Convert PD bands to percentage points for all outputs.
num_cols <- c("lower90", "lower68", "median", "upper68", "upper90")
pd_all[, (num_cols) := lapply(.SD, function(x) 100 * x), .SDcols = num_cols]

fwrite(pd_all, file.path(ss_out_dir, "irf", "pd_girf_bands_by_window_pp.csv"))
fwrite(zz_all, file.path(ss_out_dir, "irf", "psiZ_bands_by_window.csv"))

# -----------------------------------------------------------------------------
# 3. Summary table (CSV + LaTeX)
# -----------------------------------------------------------------------------

peak_stats <- pd_all[shock == "one_sd",
                     .SD[which.max(abs(median))],
                     by = window][, .(window, peak_pp = median, peak_h = horizon,
                                      band68_at_peak = upper68 - lower68)]
summary_dt <- merge(meta[, .(window, label_plain, Td, n_sat, p_ttc, rho, r2)],
                    peak_stats, by = "window")
full_band <- summary_dt[window == min(window), band68_at_peak]
summary_dt[, band_ratio_vs_full := band68_at_peak / full_band]
fwrite(summary_dt, file.path(ss_out_dir, "tables", "short_sample_summary.csv"))

fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)
tex_path <- file.path(ss_out_dir, "tables", "short_sample_summary.tex")
con <- file(tex_path, open = "wt")
cat("\\begin{table}[!htbp]\n\\centering\n", file = con)
cat("\\caption{Asymmetric-sample exercise: shortening the default sample, keeping the macro-financial VAR}\n", file = con)
cat("\\label{tab:dralacbn_short_sample}\n", file = con)
cat("\\begin{tabular}{lrrrrrrr}\n\\toprule\n", file = con)
cat("Default sample & $T_d$ & $\\hat p$ (\\%) & $\\hat\\rho$ & Sat. $R^2$ & Peak $\\Delta$PD (pp) & $h$ & 68\\% width ratio \\\\\n\\midrule\n", file = con)
for (i in seq_len(nrow(summary_dt))) {
  cat(sprintf("%s & %d & %s & %s & %s & %s & %d & %s \\\\\n",
              gsub("Td = ", "$T_d$ = ", summary_dt$label_plain[i]),
              summary_dt$Td[i],
              fmt(100 * summary_dt$p_ttc[i], 2),
              fmt(summary_dt$rho[i], 3),
              fmt(summary_dt$r2[i], 3),
              fmt(summary_dt$peak_pp[i], 3),
              as.integer(summary_dt$peak_h[i]),
              fmt(summary_dt$band_ratio_vs_full[i], 2)),
      file = con)
}
cat("\\bottomrule\n\\end{tabular}\n", file = con)
cat("\\begin{tablenotes}\n\\small\n", file = con)
cat("\\item Notes: The macro-financial VAR is estimated once on the full sample; only the default sample used to reconstruct the systematic factor and estimate the satellite is truncated. Within each window, $\\hat p$, $\\hat\\rho$ and the satellite coefficients are re-estimated from scratch (the satellite design is fixed at the full-sample best-BIC specification). Peak $\\Delta$PD refers to the posterior-median response to a one-s.d. GPR shock, in percentage points. The last column reports the 68\\% credible-band width at the peak relative to the full-sample window.\n", file = con)
cat("\\end{tablenotes}\n\\end{table}\n", file = con)
close(con)

# -----------------------------------------------------------------------------
# 4. Figures
# -----------------------------------------------------------------------------

pd_1sd <- pd_all[shock == "one_sd"]

# (a) MAIN FIGURE: four facets, common y-scale, bands per window.
p_facets <- ggplot(pd_1sd, aes(x = horizon, y = median)) +
  geom_ribbon(aes(ymin = lower90, ymax = upper90), fill = "#F3D6E3", alpha = 0.65) +
  geom_ribbon(aes(ymin = lower68, ymax = upper68), fill = "#AB4A7D", alpha = 0.42) +
  geom_hline(yintercept = 0, color = "grey45", linewidth = 0.35, linetype = "dashed") +
  geom_line(color = "#6F1732", linewidth = 0.95) +
  facet_wrap(~ window_label, nrow = 1) +
  scale_x_continuous(breaks = seq(0, 12, by = 3), minor_breaks = NULL) +
  labs(
    x = "Horizon (quarters)",
    y = expression(Delta ~ "PD (percentage points)"),
    caption = "One-s.d. GPR shock. Posterior medians with 68% (dark) and 90% (light) credible bands. Common y-scale across panels. The macro-financial VAR is identical in all panels; only the default sample is truncated."
  ) +
  theme_paper()

save_fig(p_facets, "fig_short_sample_pd_facets", width = 10.5, height = 3.6)

# (b) Overlay of posterior medians, full-sample band as reference.
ref_band <- pd_1sd[window == min(SS_START_YEARS)]
col_map <- stats::setNames(COL_WINDOWS[seq_along(win_levels)], win_levels)

p_overlay <- ggplot(pd_1sd, aes(x = horizon, y = median, color = window_label)) +
  geom_ribbon(
    data = ref_band,
    aes(x = horizon, ymin = lower68, ymax = upper68),
    inherit.aes = FALSE, fill = "#AB4A7D", alpha = 0.14
  ) +
  geom_hline(yintercept = 0, color = "grey45", linewidth = 0.35, linetype = "dashed") +
  geom_line(linewidth = 1.0) +
  geom_point(size = 1.5) +
  scale_color_manual(values = col_map) +
  scale_x_continuous(breaks = seq(0, 12, by = 1), minor_breaks = NULL) +
  labs(
    x = "Horizon (quarters)",
    y = expression(Delta ~ "PD (percentage points)"),
    caption = "One-s.d. GPR shock. Posterior medians by default-sample window; shaded area: 68% band of the full-sample (1985-2024) window."
  ) +
  theme_paper()

save_fig(p_overlay, "fig_short_sample_pd_overlay", width = 7.2, height = 4.4)

# (c) Band width by horizon: the cost of a short default sample, made explicit.
width_dt <- pd_1sd[, .(window_label, horizon, width68 = upper68 - lower68)]
p_width <- ggplot(width_dt, aes(x = horizon, y = width68, color = window_label)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 1.5) +
  scale_color_manual(values = col_map) +
  scale_x_continuous(breaks = seq(0, 12, by = 1), minor_breaks = NULL) +
  labs(
    x = "Horizon (quarters)",
    y = expression("Width of the 68% credible band for" ~ Delta ~ "PD (pp)"),
    caption = "One-s.d. GPR shock. Estimation uncertainty rises as the default sample shrinks, while median responses remain comparable."
  ) +
  theme_paper()

save_fig(p_width, "fig_short_sample_band_width", width = 7.2, height = 4.4)

# (d) Context: the delinquency series with window starts.
cut_dt <- data.table(x = as.Date(paste0(SS_START_YEARS[-1], "-01-01")))
p_series <- ggplot(delinq, aes(x = Date, y = value)) +
  geom_line(color = "#6F1732", linewidth = 0.8) +
  geom_vline(data = cut_dt, aes(xintercept = x),
             linetype = "dashed", color = "#6B6B6B", linewidth = 0.5) +
  annotate("text", x = cut_dt$x, y = max(delinq$value) * 0.97,
           label = paste0(SS_START_YEARS[-1], "-"),
           hjust = -0.1, size = 3.2, color = "grey30") +
  labs(
    x = "Date", y = "Delinquency rate (%)",
    caption = "DRALACBN, quarterly. Dashed lines: start of each truncated default window."
  ) +
  theme_paper()

save_fig(p_series, "fig_short_sample_windows", width = 7.2, height = 3.6)

# -----------------------------------------------------------------------------
# 5. Save everything + README
# -----------------------------------------------------------------------------

saveRDS(
  list(info_set = SS_INFO_SET, windows = SS_START_YEARS, meta = meta,
       summary = summary_dt, pd_bands_pp = pd_all, psiZ_bands = zz_all,
       shock_2001 = shock_2001$target),
  file.path(ss_out_dir, "short_sample_results.rds")
)

readme <- c(
  "# DRALACBN short-sample (asymmetric-sample) exercise",
  "",
  "Demonstrates the modular two-block architecture: the macro-financial BVAR is",
  "estimated once on the full history, while the default sample used for the",
  "credit block (Z reconstruction, Merton parameters, satellite estimation) is",
  "progressively truncated (1985-, 1995-, 2005-, 2015-2024).",
  "",
  "Within each window, p_ttc, rho and the satellite coefficients are re-estimated",
  "from scratch; the satellite DESIGN is fixed at the full-sample best-BIC",
  "specification so that the exercise isolates the sample-length effect.",
  "",
  "Generated by code/extra/dralacbn/05_dralacbn_short_sample.R (standalone).",
  "",
  "Contents:",
  "  figures/fig_short_sample_pd_facets.(png|pdf)   main figure: 4 panels, common y-scale",
  "  figures/fig_short_sample_pd_overlay.(png|pdf)  medians overlay + full-sample band",
  "  figures/fig_short_sample_band_width.(png|pdf)  68% band width by horizon",
  "  figures/fig_short_sample_windows.(png|pdf)     delinquency series + window starts",
  "  tables/short_sample_summary.(csv|tex)          Td, p, rho, R2, peak, band ratio",
  "  irf/pd_girf_bands_by_window_pp.csv             Delta PD bands (pp), both shocks",
  "  irf/psiZ_bands_by_window.csv                   psi_Z bands, both shocks",
  "  short_sample_results.rds                       all objects"
)
writeLines(readme, file.path(ss_out_dir, "README.md"))

cat("\n--- Short-sample summary ----------------------------------------------\n")
print(summary_dt[, .(label_plain, Td, p_ttc = round(p_ttc, 4), rho = round(rho, 4),
                     r2 = round(r2, 3), peak_pp = round(peak_pp, 4),
                     peak_h, band_ratio = round(band_ratio_vs_full, 2))])
cat("Outputs written to: ", ss_out_dir, "\n", sep = "")
cat("------------------------------------------------------------------------\n")

invisible(TRUE)
