# =============================================================================
#  code/dralacbn/07_dralacbn_short_sample_bma.R
#
#  DRALACBN -- OPERATIONAL short-default-sample demonstration.
#
#  Real-world setting: a practitioner observes only a short default series
#  (regulatory redefinitions, portfolio restructuring, new institution) but a
#  long macro-financial history. This script shows that the paper's pipeline
#  applies UNCHANGED in that situation: within each default window (1985-,
#  1995-, 2005-, 2015-2024), the entire credit block is re-estimated from
#  scratch -- Merton parameters, full BIC screen of candidate satellite
#  designs, model-averaged (BMA) mixture satellite -- with the same settings
#  as the paper's main specification. The macro-financial BVAR is unchanged
#  (full history): this is the asymmetric-sample architecture at work.
#
#  Rho policy (SS_RHO_POLICY, set before sourcing to override):
#    "recalibrate"  (default) p_ttc and rho re-estimated within each window.
#                   Naive practice: on short, calm windows the unit-variance
#                   calibration severely understates rho (systematic risk).
#    "full_sample"  rho fixed at its full-sample value; p_ttc stays
#                   window-specific. Recommended practice when T_d is short
#                   (mirrors the use of regulatory/long-run asset correlations).
#  Outputs are written to separate folders so both variants can coexist; when
#  both have been run, a cross-policy comparison figure is produced.
#
#  The fixed-design exercise (05) is the controlled counterpart: it isolates
#  the pure sample-length effect by freezing the satellite design.
#
#  Run from the repository root, ideally AFTER 05 (for the comparison figure):
#      source("code/dralacbn/07_dralacbn_short_sample_bma.R")
#
#  Outputs: output/applications/dralacbn/short_sample/bma/          (recalibrate)
#           output/applications/dralacbn/short_sample/bma_rho_fixed/ (full_sample)
#  WARNING: the BIC screen enumerates ~250k designs per window; expect a
#  runtime similar to ~4 information sets of the main pipeline.
# =============================================================================

source("code/00_setup.R")
source("code/shared/_helpers_var_merton.R")
source("code/shared/_model_averaging.R")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

# -----------------------------------------------------------------------------
# 0. Parameters and output tree
# -----------------------------------------------------------------------------

SS_INFO_SET    <- if (exists("DELINQUENCY_BASELINE_INFO_SET")) DELINQUENCY_BASELINE_INFO_SET else "real_side"
SS_START_YEARS <- c(1985, 1995, 2005, 2015)

# Rho policy: "recalibrate" (naive, default) or "full_sample" (recommended
# when T_d is short). Set SS_RHO_POLICY before sourcing to override.
if (!exists("SS_RHO_POLICY")) SS_RHO_POLICY <- "recalibrate"
SS_RHO_POLICY <- match.arg(SS_RHO_POLICY, c("recalibrate", "full_sample"))

p_lags       <- if (exists("P_LAGS")) P_LAGS else 2
nrep         <- if (exists("NREP")) NREP else 30000
M_target     <- if (exists("M_TARGET")) M_TARGET else 10000L
H            <- if (exists("HORIZON")) HORIZON else 12
seed         <- if (exists("SEED_ROBUSTNESS")) SEED_ROBUSTNESS else 123
seed_sat     <- if (exists("SEED_ZFACTOR")) SEED_ZFACTOR else 12345
impulse_name <- if (exists("IMPULSE_NAME")) IMPULSE_NAME else "log_GPRD"
NB_LAGS_SAT  <- if (exists("ALT_SPEC_NB_LAGS_SAT")) ALT_SPEC_NB_LAGS_SAT else 4
MAX_VARS_SAT <- if (exists("ALT_SPEC_MAX_VARIABLES_IN_MODEL")) ALT_SPEC_MAX_VARIABLES_IN_MODEL else 6
TOPN_AIC     <- if (exists("PRESELECT_BY_AIC_TOPN")) PRESELECT_BY_AIC_TOPN else 1000L

BMA_W   <- if (exists("BMA_WEIGHT")) BMA_WEIGHT else "bic"
BMA_R   <- if (exists("BMA_RULE")) BMA_RULE else "cum_weight"
BMA_CUM <- if (exists("BMA_CUM_THRESHOLD")) BMA_CUM_THRESHOLD else 0.95
BMA_TN  <- if (exists("BMA_TOP_N")) BMA_TOP_N else 10L
BMA_DM  <- if (exists("BMA_DELTA_MAX")) BMA_DELTA_MAX else 2

bma_out_dir <- file.path(
  "output", "applications", "dralacbn", "short_sample",
  if (SS_RHO_POLICY == "full_sample") "bma_rho_fixed" else "bma"
)
for (sub in c("", "tables", "figures", "irf")) {
  dir.create(file.path(bma_out_dir, sub), recursive = TRUE, showWarnings = FALSE)
}

COL_WINDOWS <- c("#6F1732", "#AB4A7D", "#D78FB4", "#6B6B6B")

theme_paper <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      axis.title = element_text(color = "#1A1A1A"),
      axis.text = element_text(color = "grey20"),
      strip.text = element_text(face = "bold", color = "#1A1A1A", size = base_size),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linetype = "dotted", linewidth = 0.25),
      plot.caption = element_blank()
    )
}

save_fig <- function(plot, stem, width, height, dpi = 300) {
  ggsave(file.path(bma_out_dir, "figures", paste0(stem, ".png")),
         plot, width = width, height = height, dpi = dpi)
  tryCatch(
    ggsave(file.path(bma_out_dir, "figures", paste0(stem, ".pdf")),
           plot, width = width, height = height, device = grDevices::pdf),
    error = function(e) tryCatch(
      ggsave(file.path(bma_out_dir, "figures", paste0(stem, ".pdf")),
             plot, width = width, height = height),
      error = function(e2) message("PDF export failed for ", stem)
    )
  )
}

cat("\n========== DRALACBN short-sample exercise, per-window BMA ==========\n")
cat("Rho policy: ", SS_RHO_POLICY, "\n", sep = "")

# -----------------------------------------------------------------------------
# 1. Data, BVAR kernel (unchanged), default proxy
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

# Full-sample Merton calibration (anchor for the "full_sample" rho policy).
z_full <- f_Z_estimation(delinq$value / 100)
RHO_FULL <- z_full$rho
if (SS_RHO_POLICY == "full_sample") {
  cat("Rho fixed at its full-sample value: ", sprintf("%.4f", RHO_FULL),
      " (p_ttc stays window-specific)\n", sep = "")
}

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
vars_sat <- setdiff(variables, impulse_name)

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
# 2. One window: full credit-block re-estimation (Merton + screen + BMA)
# -----------------------------------------------------------------------------

run_window_bma <- function(start_year) {
  cat("\n----- Window ", start_year, "-2024: full BIC screen + BMA -----\n", sep = "")

  d_w <- delinq[Date >= as.Date(paste0(start_year, "-01-01"))]
  z_out <- if (SS_RHO_POLICY == "full_sample") {
    f_Z_estimation(d_w$value / 100, rho_fixed = RHO_FULL)
  } else {
    f_Z_estimation(d_w$value / 100)
  }
  if (!is.finite(z_out$rho)) {
    stop("Rho calibration failed for window ", start_year, call. = FALSE)
  }

  Z_DT_w <- data.table(
    Date = as.Date(format(d_w$Date + 62, "%Y-%m-01")),
    Z = z_out$Z
  )
  sat_df_w <- prepare_lagged_satellite_data(
    z_input = Z_DT_w, DT_raw = DT_raw, vars_sat = vars_sat, nb_lags = NB_LAGS_SAT
  )

  screen <- select_satellite_fullsample(
    model_df = sat_df_w,
    max_vars = MAX_VARS_SAT,
    p_threshold = NULL,
    preselect_by_aic_topn = TOPN_AIC
  )

  bsat_w <- make_bma_satellite(
    candidates = screen$candidates,
    X_names = screen$X_names,
    model_df = sat_df_w,
    variables = variables,
    weight = BMA_W, rule = BMA_R,
    cum_threshold = BMA_CUM, top_n = BMA_TN, delta_max = BMA_DM,
    n_draws = M_target,
    seed = seed_sat
  )

  best_terms <- screen$best$vars
  cat("  Best model: ", paste(best_terms, collapse = " + "), "\n", sep = "")
  cat("  Models averaged: ", bsat_w$bma$n_models,
      " | p_ttc = ", sprintf("%.4f", z_out$p_ttc),
      " | rho = ", sprintf("%.4f", z_out$rho), "\n", sep = "")

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

  incl <- copy(bsat_w$bma$inclusion)
  incl[, window := start_year]
  models <- copy(bsat_w$bma$models)
  models[, window := start_year]

  list(
    start_year = start_year, Td = nrow(d_w),
    p_ttc = z_out$p_ttc, rho = z_out$rho,
    n_models = bsat_w$bma$n_models,
    best_model = paste(best_terms, collapse = " + "),
    inclusion = incl, models = models,
    psiZ = rbind(s1$psiZ, s2$psiZ),
    pd = rbind(s1$pd, s2$pd)
  )
}

runs <- lapply(SS_START_YEARS, run_window_bma)
names(runs) <- as.character(SS_START_YEARS)

# -----------------------------------------------------------------------------
# 3. Collect and export
# -----------------------------------------------------------------------------

pd_all <- rbindlist(lapply(runs, function(r) r$pd))
zz_all <- rbindlist(lapply(runs, function(r) r$psiZ))
incl_all <- rbindlist(lapply(runs, function(r) r$inclusion), fill = TRUE)
models_all <- rbindlist(lapply(runs, function(r) r$models), fill = TRUE)

meta <- rbindlist(lapply(runs, function(r) {
  data.table(window = r$start_year, Td = r$Td, p_ttc = r$p_ttc, rho = r$rho,
             n_models = r$n_models, best_model = r$best_model)
}))
meta[, label_plain := sprintf("%d-2024 (Td = %d)", window, Td)]
win_levels <- meta$label_plain

add_labels <- function(dt) {
  dt <- merge(dt, meta[, .(window, label_plain)], by = "window", sort = FALSE)
  dt[, window_label := factor(label_plain, levels = win_levels)]
  dt
}
pd_all <- add_labels(pd_all)
zz_all <- add_labels(zz_all)

num_cols <- c("lower90", "lower68", "median", "upper68", "upper90")
pd_all[, (num_cols) := lapply(.SD, function(x) 100 * x), .SDcols = num_cols]  # pp

fwrite(pd_all, file.path(bma_out_dir, "irf", "pd_girf_bands_by_window_bma_pp.csv"))
fwrite(zz_all, file.path(bma_out_dir, "irf", "psiZ_bands_by_window_bma.csv"))
fwrite(incl_all, file.path(bma_out_dir, "tables", "inclusion_probabilities_by_window.csv"))
fwrite(models_all, file.path(bma_out_dir, "tables", "top_models_by_window.csv"))
fwrite(meta, file.path(bma_out_dir, "tables", "bma_summary_by_window.csv"))

# -----------------------------------------------------------------------------
# 4. Figures
# -----------------------------------------------------------------------------

pd_1sd <- pd_all[shock == "one_sd"]

# (a) PD GIRF facets, common y-scale (mirror of the fixed-design main figure).
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
    caption = "One-s.d. GPR shock. Per-window BIC screen and model averaging (selection redone within each window). Posterior medians with 68% and 90% bands; common y-scale."
  ) +
  theme_paper() +
  theme(legend.position = "none")

save_fig(p_facets, "fig_short_sample_bma_pd_facets", width = 10.5, height = 3.6)

# (b) BMA vs fixed-design medians, if the 05 results are available.
fixed_path <- file.path("output", "applications", "dralacbn",
                        "short_sample", "short_sample_results.rds")
if (file.exists(fixed_path)) {
  fixed <- readRDS(fixed_path)
  pd_fixed <- as.data.table(fixed$pd_bands_pp)[shock == "one_sd"]
  pd_fixed <- pd_fixed[, .(window, horizon, median, satellite = "Fixed design (best-BIC)")]
  pd_bma_m <- pd_1sd[, .(window, horizon, median, satellite = "Per-window BMA")]
  both <- rbind(pd_fixed, pd_bma_m)
  both <- merge(both, meta[, .(window, label_plain)], by = "window", sort = FALSE)
  both[, window_label := factor(label_plain, levels = win_levels)]

  p_vs <- ggplot(both, aes(x = horizon, y = median,
                           color = satellite, linetype = satellite)) +
    geom_hline(yintercept = 0, color = "grey45", linewidth = 0.35, linetype = "dashed") +
    geom_line(linewidth = 0.95) +
    facet_wrap(~ window_label, nrow = 1) +
    scale_color_manual(values = c("Per-window BMA" = "#6F1732",
                                  "Fixed design (best-BIC)" = "#6B6B6B")) +
    scale_linetype_manual(values = c("Per-window BMA" = "solid",
                                     "Fixed design (best-BIC)" = "longdash")) +
    scale_x_continuous(breaks = seq(0, 12, by = 3), minor_breaks = NULL) +
    labs(
      x = "Horizon (quarters)",
      y = expression(Delta ~ "PD (percentage points)"),
      caption = "One-s.d. GPR shock, posterior medians. The gap between the two lines isolates the model-selection instability added by re-running the BIC screen on the short window."
    ) +
    theme_paper()

  save_fig(p_vs, "fig_short_sample_bma_vs_fixed", width = 10.5, height = 3.6)
} else {
  message("05 results not found; skipping the BMA-vs-fixed comparison figure.")
}

# (c) Selection instability: inclusion probabilities across windows.
incl_plot <- copy(incl_all)
incl_plot <- merge(incl_plot, meta[, .(window, label_plain)], by = "window", sort = FALSE)
incl_plot[, window_label := factor(label_plain, levels = win_levels)]
term_order <- incl_plot[window == min(SS_START_YEARS)][order(inclusion_prob)]$term
extra_terms <- setdiff(unique(incl_plot$term), term_order)
incl_plot[, term_f := factor(term, levels = c(extra_terms, term_order))]

p_incl <- ggplot(incl_plot, aes(x = window_label, y = term_f,
                                fill = inclusion_prob)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = formatC(inclusion_prob, format = "f", digits = 2)),
            size = 2.9, color = "grey15") +
  scale_fill_gradient(low = "#FBEFF5", high = "#AB4A7D",
                      limits = c(0, 1), name = "Inclusion prob.") +
  labs(
    x = "Default-sample window", y = NULL,
    caption = "Posterior inclusion probabilities (cumulative BIC weights) of each satellite term, per default window. Terms absent from a window's retained set have probability ~0."
  ) +
  theme_paper() +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 8))

n_terms <- length(unique(incl_plot$term))
save_fig(p_incl, "fig_short_sample_bma_inclusion",
         width = 8.0, height = max(3.5, 0.28 * n_terms + 1.5))

# (d) Cross-policy comparison (naive rho recalibration vs full-sample rho),
#     produced once both variants of this script have been run.
other_policy <- if (SS_RHO_POLICY == "full_sample") "bma" else "bma_rho_fixed"
other_rds <- file.path("output", "applications", "dralacbn", "short_sample",
                       other_policy, "short_sample_bma_results.rds")
if (file.exists(other_rds)) {
  other <- readRDS(other_rds)
  lab_this <- if (SS_RHO_POLICY == "full_sample") {
    "Rho fixed (full sample)"
  } else {
    "Rho recalibrated per window"
  }
  lab_other <- if (SS_RHO_POLICY == "full_sample") {
    "Rho recalibrated per window"
  } else {
    "Rho fixed (full sample)"
  }

  cmp_this <- pd_1sd[, .(window, horizon, median, lower68, upper68,
                         policy = lab_this)]
  cmp_other <- as.data.table(other$pd_bands_pp)[shock == "one_sd",
    .(window, horizon, median, lower68, upper68, policy = lab_other)]
  cmp <- rbind(cmp_this, cmp_other)
  cmp <- merge(cmp, meta[, .(window, label_plain)], by = "window", sort = FALSE)
  cmp[, window_label := factor(label_plain, levels = win_levels)]
  cmp[, policy := factor(policy, levels = c("Rho recalibrated per window",
                                            "Rho fixed (full sample)"))]

  p_pol <- ggplot(cmp, aes(x = horizon, y = median,
                           color = policy, fill = policy, linetype = policy)) +
    geom_ribbon(aes(ymin = lower68, ymax = upper68), alpha = 0.15, color = NA) +
    geom_hline(yintercept = 0, color = "grey45", linewidth = 0.35, linetype = "dashed") +
    geom_line(linewidth = 0.95) +
    facet_wrap(~ window_label, nrow = 1) +
    scale_color_manual(values = c("Rho recalibrated per window" = "#6B6B6B",
                                  "Rho fixed (full sample)" = "#6F1732")) +
    scale_fill_manual(values = c("Rho recalibrated per window" = "#6B6B6B",
                                 "Rho fixed (full sample)" = "#AB4A7D")) +
    scale_linetype_manual(values = c("Rho recalibrated per window" = "longdash",
                                     "Rho fixed (full sample)" = "solid")) +
    scale_x_continuous(breaks = seq(0, 12, by = 3), minor_breaks = NULL) +
    labs(
      x = "Horizon (quarters)",
      y = expression(Delta ~ "PD (percentage points)"),
      caption = "One-s.d. GPR shock; posterior medians with 68% bands. Recalibrating rho on short, calm windows compresses the PD scale and understates the systematic-risk response; fixing rho at its long-run value restores comparability."
    ) +
    theme_paper()

  save_fig(p_pol, "fig_short_sample_rho_policy_comparison",
           width = 10.5, height = 3.9)
  cat("Cross-policy comparison figure written (both rho policies found).\n")
}

# -----------------------------------------------------------------------------
# 5. Save everything + README
# -----------------------------------------------------------------------------

saveRDS(
  list(info_set = SS_INFO_SET, windows = SS_START_YEARS,
       rho_policy = SS_RHO_POLICY, rho_full_sample = RHO_FULL,
       bma_settings = list(weight = BMA_W, rule = BMA_R, cum = BMA_CUM),
       meta = meta, inclusion = incl_all, models = models_all,
       pd_bands_pp = pd_all, psiZ_bands = zz_all),
  file.path(bma_out_dir, "short_sample_bma_results.rds")
)

readme <- c(
  "# DRALACBN short-sample exercise -- operational per-window BMA",
  "",
  "Operational demonstration: a practitioner with a SHORT default series and a",
  "long macro-financial history applies the paper's pipeline UNCHANGED. Within",
  "each default window, the Merton parameters, the full BIC screen of candidate",
  "satellite designs AND the model-averaged satellite are re-estimated from",
  "scratch, with the main-specification settings. The BVAR is unchanged.",
  "",
  paste0("Rho policy of this run: ", SS_RHO_POLICY,
         " (rho_full_sample = ", sprintf("%.4f", RHO_FULL), ")."),
  "Run the script once with SS_RHO_POLICY <- \"recalibrate\" and once with",
  "SS_RHO_POLICY <- \"full_sample\" to obtain the cross-policy comparison figure",
  "(fig_short_sample_rho_policy_comparison): recalibrating rho on short, calm",
  "windows understates systematic risk; fixing rho at its long-run value is the",
  "recommended practice (mirrors regulatory long-run asset correlations).",
  "",
  "The fixed-design exercise (05) is the controlled counterpart: it isolates",
  "the pure sample-length effect (fig_short_sample_bma_vs_fixed shows the gap).",
  "",
  "Generated by code/dralacbn/07_dralacbn_short_sample_bma.R (standalone).",
  "",
  "Contents:",
  "  figures/fig_short_sample_bma_pd_facets           PD GIRF per window (BMA satellite)",
  "  figures/fig_short_sample_bma_vs_fixed            BMA vs fixed-design medians",
  "  figures/fig_short_sample_bma_inclusion           inclusion-probability heatmap",
  "  figures/fig_short_sample_rho_policy_comparison   naive vs recommended rho (if both runs exist)",
  "  tables/inclusion_probabilities_by_window.csv",
  "  tables/top_models_by_window.csv",
  "  tables/bma_summary_by_window.csv",
  "  irf/*.csv                                        PD (pp) and psi_Z bands, both shocks",
  "  short_sample_bma_results.rds                     all objects"
)
writeLines(readme, file.path(bma_out_dir, "README.md"))

cat("\n--- Per-window BMA summary ---------------------------------------------\n")
print(meta[, .(label_plain, n_models, p_ttc = round(p_ttc, 4),
               rho = round(rho, 4), best_model)])
cat("Outputs written to: ", bma_out_dir, "\n", sep = "")
cat("------------------------------------------------------------------------\n")

invisible(TRUE)
