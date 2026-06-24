# =============================================================================
#  code/dralacbn/04_dralacbn_direct_channel_paper_outputs.R
#
#  Publication-quality figures, tables and clean CSVs for the direct-channel
#  (control-function) exercise. Reads the results of
#  code/dralacbn/03_dralacbn_direct_channel.R and writes paper-ready outputs.
#
#  STANDALONE: modifies nothing in the existing pipeline; not in run_all.R.
#  Run from the repository root AFTER 03_dralacbn_direct_channel.R:
#      source("code/dralacbn/04_dralacbn_direct_channel_paper_outputs.R")
#
#  Outputs: output/applications/dralacbn/direct_channel/paper/
#      figures/  PNG (300 dpi) + PDF versions
#      tables/   LaTeX (booktabs) + clean CSVs
#      data/     rounded, self-describing CSVs (units: Delta PD in pp)
# =============================================================================

source("code/00_setup.R")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

res_path <- file.path("output", "applications", "dralacbn",
                      "direct_channel", "direct_channel_results.rds")
if (!file.exists(res_path)) {
  stop("Missing ", res_path,
       "\nRun code/dralacbn/03_dralacbn_direct_channel.R first.", call. = FALSE)
}
res <- readRDS(res_path)

paper_dir <- file.path("output", "applications", "dralacbn", "direct_channel", "paper")
for (sub in c("figures", "tables", "data")) {
  dir.create(file.path(paper_dir, sub), recursive = TRUE, showWarnings = FALSE)
}

# House style ------------------------------------------------------------------
COL_DARKRED <- "#6F1732"
COL_ROSE2   <- "#AB4A7D"
COL_ROSE    <- "#F3D6E3"
COL_GREY    <- "#6B6B6B"

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
  ggsave(file.path(paper_dir, "figures", paste0(stem, ".png")),
         plot, width = width, height = height, dpi = dpi)
  ok <- tryCatch({
    ggsave(file.path(paper_dir, "figures", paste0(stem, ".pdf")),
           plot, width = width, height = height, device = grDevices::pdf)
    TRUE
  }, error = function(e) FALSE)
  if (!ok) {
    tryCatch(
      ggsave(file.path(paper_dir, "figures", paste0(stem, ".pdf")),
             plot, width = width, height = height),
      error = function(e) message("PDF export failed for ", stem, ": ",
                                  conditionMessage(e))
    )
  }
}

fmt <- function(x, d = 4) formatC(x, format = "f", digits = d)

lam      <- res$lambda_draws
lam_hat  <- res$lambda_hat_ols
ls       <- res$lambda_summary
hac      <- res$hac_diagnostic
joint    <- res$joint_orthogonality_test
comp     <- as.data.table(res$coefficients_comparison)
pd_bands <- as.data.table(res$pd_bands)
zz_bands <- as.data.table(res$psiZ_bands)

# Variance decomposition of lambda: generated-regressor (between-VAR-draw)
# vs within-draw posterior noise. lambda_m = lambda_hat_m + conjugate noise.
v_total   <- stats::var(lam)
v_between <- stats::var(lam_hat)
v_within  <- stats::var(lam - lam_hat)
share_between <- 100 * v_between / v_total

# =============================================================================
# 1. Figure: posterior density of lambda, annotated
# =============================================================================

dens <- stats::density(lam)
dens_dt <- data.table(x = dens$x, y = dens$y)
ci90 <- c(ls$q05, ls$q95)

p_lambda <- ggplot(dens_dt, aes(x = x, y = y)) +
  geom_area(data = dens_dt[x >= ci90[1] & x <= ci90[2]],
            fill = COL_ROSE, alpha = 0.9) +
  geom_line(color = COL_DARKRED, linewidth = 0.9) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey30", linewidth = 0.45) +
  geom_vline(xintercept = ls$median, color = COL_DARKRED, linewidth = 0.6) +
  annotate("text", x = ls$median, y = max(dens_dt$y) * 1.045,
           label = paste0("median = ", fmt(ls$median, 3)),
           color = COL_DARKRED, size = 3.4, hjust = 1.05) +
  annotate("text", x = max(dens_dt$x), y = max(dens_dt$y) * 0.92,
           label = paste0("Pr(lambda < 0) = ", fmt(ls$prob_negative, 2)),
           color = "grey25", size = 3.4, hjust = 1) +
  labs(
    x = expression(lambda ~ "(impact response of" ~ Z ~ "to a one-s.d. GPR innovation)"),
    y = "Posterior density",
    caption = "Shaded area: 90% credible interval. Posterior pools 10,000 hierarchical draws (one conjugate satellite draw per BVAR draw)."
  ) +
  theme_paper()

save_fig(p_lambda, "fig_direct_channel_lambda_posterior", width = 6.8, height = 4.0)

# =============================================================================
# 2. Figure: Delta PD overlay, restricted vs augmented, both shocks (in pp)
# =============================================================================

pd <- copy(pd_bands)
num_cols <- c("lower90", "lower68", "median", "upper68", "upper90")
pd[, (num_cols) := lapply(.SD, function(x) 100 * x), .SDcols = num_cols]  # decimal -> pp

pd[, model_label := fifelse(model == "augmented",
                            "Augmented (direct channel)",
                            "Restricted (baseline)")]
pd[, model_label := factor(model_label,
                           levels = c("Restricted (baseline)",
                                      "Augmented (direct channel)"))]
pd[, shock_label := fifelse(shock == "one_sd",
                            "One-s.d. GPR shock", "2001Q3-sized GPR shock")]
pd[, shock_label := factor(shock_label,
                           levels = c("One-s.d. GPR shock",
                                      "2001Q3-sized GPR shock"))]

p_pd <- ggplot(pd, aes(x = horizon, y = median,
                       color = model_label, fill = model_label,
                       linetype = model_label)) +
  geom_ribbon(aes(ymin = lower68, ymax = upper68), alpha = 0.16, color = NA) +
  geom_hline(yintercept = 0, color = "grey45", linewidth = 0.35, linetype = "dashed") +
  geom_line(linewidth = 0.95) +
  facet_wrap(~ shock_label, nrow = 1, scales = "free_y") +
  scale_color_manual(values = c("Restricted (baseline)" = COL_GREY,
                                "Augmented (direct channel)" = COL_DARKRED)) +
  scale_fill_manual(values = c("Restricted (baseline)" = COL_GREY,
                               "Augmented (direct channel)" = COL_ROSE2)) +
  scale_linetype_manual(values = c("Restricted (baseline)" = "longdash",
                                   "Augmented (direct channel)" = "solid")) +
  scale_x_continuous(breaks = seq(0, 12, by = 2), minor_breaks = NULL) +
  labs(
    x = "Horizon (quarters)",
    y = expression(Delta ~ "PD (percentage points)"),
    caption = "Posterior medians; shaded bands are 68% pointwise credible intervals. The two responses differ only at impact (h = 0)."
  ) +
  theme_paper()

save_fig(p_pd, "fig_direct_channel_pd_girf", width = 9.2, height = 3.9)

# Single-panel version (one-s.d. shock), for a smaller figure in the text.
p_pd_1sd <- ggplot(pd[shock == "one_sd"],
                   aes(x = horizon, y = median,
                       color = model_label, fill = model_label,
                       linetype = model_label)) +
  geom_ribbon(aes(ymin = lower68, ymax = upper68), alpha = 0.16, color = NA) +
  geom_hline(yintercept = 0, color = "grey45", linewidth = 0.35, linetype = "dashed") +
  geom_line(linewidth = 0.95) +
  scale_color_manual(values = c("Restricted (baseline)" = COL_GREY,
                                "Augmented (direct channel)" = COL_DARKRED)) +
  scale_fill_manual(values = c("Restricted (baseline)" = COL_GREY,
                               "Augmented (direct channel)" = COL_ROSE2)) +
  scale_linetype_manual(values = c("Restricted (baseline)" = "longdash",
                                   "Augmented (direct channel)" = "solid")) +
  scale_x_continuous(breaks = seq(0, 12, by = 1), minor_breaks = NULL) +
  labs(x = "Horizon (quarters)", y = expression(Delta ~ "PD (percentage points)"),
       caption = "Posterior medians with 68% pointwise credible bands.") +
  theme_paper()

save_fig(p_pd_1sd, "fig_direct_channel_pd_girf_1sd", width = 6.8, height = 4.2)

# =============================================================================
# 3. Figure: psi_Z overlay (one-s.d. shock)
# =============================================================================

zz <- copy(zz_bands)
zz[, model_label := fifelse(model == "augmented",
                            "Augmented (direct channel)",
                            "Restricted (baseline)")]
zz[, model_label := factor(model_label,
                           levels = c("Restricted (baseline)",
                                      "Augmented (direct channel)"))]

p_z <- ggplot(zz[shock == "one_sd"],
              aes(x = horizon, y = median,
                  color = model_label, fill = model_label,
                  linetype = model_label)) +
  geom_ribbon(aes(ymin = lower68, ymax = upper68), alpha = 0.16, color = NA) +
  geom_hline(yintercept = 0, color = "grey45", linewidth = 0.35, linetype = "dashed") +
  geom_line(linewidth = 0.95) +
  scale_color_manual(values = c("Restricted (baseline)" = COL_GREY,
                                "Augmented (direct channel)" = COL_DARKRED)) +
  scale_fill_manual(values = c("Restricted (baseline)" = COL_GREY,
                               "Augmented (direct channel)" = COL_ROSE2)) +
  scale_linetype_manual(values = c("Restricted (baseline)" = "longdash",
                                   "Augmented (direct channel)" = "solid")) +
  scale_x_continuous(breaks = seq(0, 12, by = 1), minor_breaks = NULL) +
  labs(x = "Horizon (quarters)", y = expression(psi[Z](h) ~ "(s.d. of" ~ Z * ")"),
       caption = "Posterior medians with 68% pointwise credible bands. One-s.d. GPR shock.") +
  theme_paper()

save_fig(p_z, "fig_direct_channel_z_girf_1sd", width = 6.8, height = 4.2)

# =============================================================================
# 4. LaTeX table: restricted vs augmented satellite, with lambda block
# =============================================================================

cell <- function(m, s, d = 4) {
  if (!is.finite(m)) return("--")
  paste0(fmt(m, d), " (", fmt(s, d), ")")
}

var_pretty <- c(
  "(Intercept)"        = "Intercept",
  "log_inv_pc_lag0"    = "log(Investment p.c.)$_{t}$",
  "log_inv_pc_lag1"    = "log(Investment p.c.)$_{t-1}$",
  "log_inv_pc_lag2"    = "log(Investment p.c.)$_{t-2}$",
  "log_gdp_pc_lag2"    = "log(GDP p.c.)$_{t-2}$",
  "log_gdp_pc_lag4"    = "log(GDP p.c.)$_{t-4}$",
  "log_oil_real_lag0"  = "log(Oil real)$_{t}$",
  "log_oil_real_lag4"  = "log(Oil real)$_{t-4}$",
  "infl_yoy_pct_lag4"  = "Inflation YoY$_{t-4}$",
  "gpr_innov"          = "$e^{GPR}_{t}$ (direct channel, $\\lambda$)"
)
pretty_term <- function(x) {
  out <- unname(var_pretty[x])
  out[is.na(out)] <- gsub("_", "\\\\_", x[is.na(out)])
  out
}

tex_path <- file.path(paper_dir, "tables", "tab_direct_channel_satellite.tex")
con <- file(tex_path, open = "wt")
cat("\\begin{table}[!htbp]\n\\centering\n", file = con)
cat("\\caption{Relaxing the satellite exogeneity assumption: control-function estimates}\n", file = con)
cat("\\label{tab:direct_channel_satellite}\n", file = con)
cat("\\begin{tabular}{lcc}\n\\toprule\n", file = con)
cat(" & Restricted ($\\lambda = 0$) & Augmented \\\\\n\\midrule\n", file = con)
for (i in seq_len(nrow(comp))) {
  cat(pretty_term(comp$term[i]), " & ",
      cell(comp$mean_restricted[i], comp$sd_restricted[i]), " & ",
      cell(comp$mean_augmented[i], comp$sd_augmented[i]), " \\\\\n",
      sep = "", file = con)
}
cat("\\midrule\n", file = con)
cat(sprintf("90\\%% CI for $\\lambda$ & -- & [%s, %s] \\\\\n",
            fmt(ls$q05, 3), fmt(ls$q95, 3)), file = con)
cat(sprintf("$\\Pr(\\lambda < 0 \\mid \\text{data})$ & -- & %s \\\\\n",
            fmt(ls$prob_negative, 2)), file = con)
cat(sprintf("HAC $t$-statistic ($p$-value) & -- & %s (%s) \\\\\n",
            fmt(hac$t_value, 2), fmt(hac$p_value, 2)), file = con)
if (!is.null(joint) && is.finite(joint$F_stat[1])) {
  cat(sprintf("Joint orthogonality $F$-test ($p$-value) & -- & %s (%s) \\\\\n",
              fmt(joint$F_stat[1], 2), fmt(joint$p_value[1], 2)), file = con)
}
cat(sprintf("Generated-regressor share of $\\mathrm{Var}(\\lambda)$ & -- & %s\\%% \\\\\n",
            fmt(share_between, 1)), file = con)
cat(sprintf("Observations & %d & %d \\\\\n", ls$n_obs, ls$n_obs), file = con)
cat(sprintf("Posterior draws & %s & %s \\\\\n",
            format(ls$n_draws, big.mark = ","),
            format(ls$n_draws, big.mark = ",")), file = con)
cat("\\bottomrule\n\\end{tabular}\n", file = con)
cat("\\begin{tablenotes}\n\\small\n", file = con)
cat("\\item Notes: Posterior means with posterior standard deviations in parentheses. The augmented satellite adds the standardized reduced-form GPR innovation $e^{GPR}_t = u_{g,t}/\\sqrt{\\sigma_{gg}}$ to the best-BIC specification, relaxing the exclusion restriction $\\eta_t \\perp u_t$ ($\\lambda = 0$ in the baseline). Estimation is hierarchical: the innovation series and the conjugate satellite posterior are recomputed for each of the retained BVAR posterior draws, so generated-regressor uncertainty is fully propagated; its share of the posterior variance of $\\lambda$ is reported above. The HAC diagnostic fixes the innovation at its posterior-median series and uses Newey--West standard errors. The joint orthogonality $F$-test regresses the baseline ($\\lambda = 0$) satellite residual on the full vector of reduced-form VAR innovations (at posterior-median coefficients) and tests their joint nullity, i.e. $\\eta_t \\perp u_t$ for the entire innovation vector, not only the GPR component.\n", file = con)
cat("\\end{tablenotes}\n\\end{table}\n", file = con)
close(con)

# =============================================================================
# 5. LaTeX table: PD responses at impact and peak (pp), both shocks
# =============================================================================

peak_row <- function(shock_id, model_id) {
  b <- pd[shock == shock_id & model == model_id]
  i <- b[, which.max(abs(median))]
  list(impact = b[horizon == 0, median],
       impact_lo = b[horizon == 0, lower68], impact_hi = b[horizon == 0, upper68],
       peak = b$median[i], peak_h = b$horizon[i],
       peak_lo = b$lower68[i], peak_hi = b$upper68[i])
}

tex2 <- file.path(paper_dir, "tables", "tab_direct_channel_pd.tex")
con <- file(tex2, open = "wt")
cat("\\begin{table}[!htbp]\n\\centering\n", file = con)
cat("\\caption{Portfolio PD responses with and without the direct geopolitical channel}\n", file = con)
cat("\\label{tab:direct_channel_pd}\n", file = con)
cat("\\begin{tabular}{llcc}\n\\toprule\n", file = con)
cat("Shock & Model & Impact $\\Delta$PD ($h=0$) & Peak $\\Delta$PD \\\\\n\\midrule\n", file = con)
for (sh in c("one_sd", "2001Q3")) {
  sh_lab <- if (sh == "one_sd") "One s.d." else "2001Q3"
  for (mo in c("restricted", "augmented")) {
    mo_lab <- if (mo == "restricted") "Restricted" else "Augmented"
    r <- peak_row(sh, mo)
    cat(sprintf("%s & %s & %s [%s, %s] & %s (h=%d) [%s, %s] \\\\\n",
                sh_lab, mo_lab,
                fmt(r$impact, 3), fmt(r$impact_lo, 3), fmt(r$impact_hi, 3),
                fmt(r$peak, 3), as.integer(r$peak_h),
                fmt(r$peak_lo, 3), fmt(r$peak_hi, 3)),
        file = con)
  }
  if (sh == "one_sd") cat("\\midrule\n", file = con)
}
cat("\\bottomrule\n\\end{tabular}\n", file = con)
cat("\\begin{tablenotes}\n\\small\n", file = con)
cat("\\item Notes: Posterior medians in percentage points, with 68\\% credible intervals in brackets. The direct channel only affects the impact response: with serially independent satellite errors, the $\\lambda$ correction to the PD generalized impulse response applies at $h = 0$, leaving the propagation profile unchanged from $h \\geq 1$.\n", file = con)
cat("\\end{tablenotes}\n\\end{table}\n", file = con)
close(con)

# =============================================================================
# 6. Clean CSVs (units documented in column names)
# =============================================================================

pd_csv <- copy(pd)[, .(
  shock, model, horizon,
  median_pp = round(median, 5),
  lower68_pp = round(lower68, 5), upper68_pp = round(upper68, 5),
  lower90_pp = round(lower90, 5), upper90_pp = round(upper90, 5)
)]
fwrite(pd_csv, file.path(paper_dir, "data", "pd_girf_direct_channel_pp.csv"))

zz_csv <- copy(zz_bands)[, .(
  shock, model, horizon,
  median_sd = round(median, 5),
  lower68_sd = round(lower68, 5), upper68_sd = round(upper68, 5),
  lower90_sd = round(lower90, 5), upper90_sd = round(upper90, 5)
)]
fwrite(zz_csv, file.path(paper_dir, "data", "z_girf_direct_channel_sd.csv"))

lambda_csv <- data.table(
  statistic = c("posterior_mean", "posterior_sd", "q05", "q16", "median",
                "q84", "q95", "prob_lambda_negative", "hac_t", "hac_p",
                "var_share_generated_regressor_pct", "n_obs", "n_draws"),
  value = c(round(ls$mean, 5), round(ls$sd, 5), round(ls$q05, 5),
            round(ls$q16, 5), round(ls$median, 5), round(ls$q84, 5),
            round(ls$q95, 5), round(ls$prob_negative, 4),
            round(hac$t_value, 3), round(hac$p_value, 4),
            round(share_between, 2), ls$n_obs, ls$n_draws)
)
fwrite(lambda_csv, file.path(paper_dir, "data", "lambda_posterior_clean.csv"))

# =============================================================================
# 7. Console recap
# =============================================================================

cat("\n--- Paper outputs for the direct-channel exercise ---------------------\n")
cat("lambda: mean ", fmt(ls$mean, 4), ", 90% CI [", fmt(ls$q05, 4), ", ",
    fmt(ls$q95, 4), "], Pr(<0) = ", fmt(ls$prob_negative, 3), "\n", sep = "")
cat("Generated-regressor share of Var(lambda): ", fmt(share_between, 1), "%\n", sep = "")
cat("Figures: ", file.path(paper_dir, "figures"), "\n", sep = "")
cat("Tables:  ", file.path(paper_dir, "tables"), "\n", sep = "")
cat("Data:    ", file.path(paper_dir, "data"), "\n", sep = "")
cat("------------------------------------------------------------------------\n")

invisible(TRUE)
