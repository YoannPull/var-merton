# =============================================================================
#  code/dralacbn/09_dralacbn_short_sample_single_figures.R
#
#  One figure PER WINDOW (not side-by-side facets) for the short-sample
#  exercises, publication style. Reads the saved results of:
#    05  output/applications/dralacbn/short_sample/short_sample_results.rds
#    07  output/applications/dralacbn/short_sample/bma/short_sample_bma_results.rds
#        output/applications/dralacbn/short_sample/bma_rho_fixed/... (if present)
#  and writes single-panel Delta PD figures (PNG 300 dpi + PDF), with a
#  COMMON y-range across the windows of each exercise so that panels remain
#  directly comparable when placed on different pages or in different floats.
#
#  Cheap (reads saved bands only). STANDALONE; not in run_all.R.
#  Run from the repository root AFTER 05 (and 07 for the BMA variants):
#      source("code/dralacbn/09_dralacbn_short_sample_single_figures.R")
#
#  Outputs:
#    .../short_sample/figures/single/          (fixed-design exercise, 05)
#    .../short_sample/bma/figures/single/      (per-window BMA, 07)
#    .../short_sample/bma_rho_fixed/figures/single/  (if the variant was run)
# =============================================================================

source("code/00_setup.R")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

theme_paper <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      axis.title = element_text(color = "#1A1A1A"),
      axis.text = element_text(color = "grey20"),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linetype = "dotted", linewidth = 0.25),
      plot.caption = element_blank()
    )
}

save_fig <- function(plot, dir, stem, width = 6.2, height = 4.0, dpi = 300) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(dir, paste0(stem, ".png")),
         plot, width = width, height = height, dpi = dpi)
  tryCatch(
    ggsave(file.path(dir, paste0(stem, ".pdf")),
           plot, width = width, height = height, device = grDevices::pdf),
    error = function(e) tryCatch(
      ggsave(file.path(dir, paste0(stem, ".pdf")),
             plot, width = width, height = height),
      error = function(e2) message("PDF export failed for ", stem)
    )
  )
}

plot_single_window <- function(b, ylim_common, caption_txt) {
  ggplot(b, aes(x = horizon, y = median)) +
    geom_ribbon(aes(ymin = lower90, ymax = upper90), fill = "#F3D6E3", alpha = 0.65) +
    geom_ribbon(aes(ymin = lower68, ymax = upper68), fill = "#AB4A7D", alpha = 0.42) +
    geom_hline(yintercept = 0, color = "grey45", linewidth = 0.35, linetype = "dashed") +
    geom_line(color = "#6F1732", linewidth = 1.0) +
    scale_x_continuous(breaks = seq(0, 12, by = 1), minor_breaks = NULL) +
    coord_cartesian(ylim = ylim_common) +
    labs(x = "Horizon (quarters)",
         y = expression(Delta ~ "PD (percentage points)"),
         caption = caption_txt) +
    theme_paper()
}

# -----------------------------------------------------------------------------
# Generic driver: one exercise -> one figure per window x shock
# -----------------------------------------------------------------------------

make_single_figures <- function(rds_path, out_dir, exercise_label) {
  if (!file.exists(rds_path)) {
    message("Skipping (not found): ", rds_path)
    return(invisible(FALSE))
  }
  res <- readRDS(rds_path)
  pd <- as.data.table(res$pd_bands_pp)
  if (!all(c("window", "shock", "horizon", "median") %in% names(pd))) {
    message("Skipping (unexpected format): ", rds_path)
    return(invisible(FALSE))
  }

  # Window labels: reuse label_plain if present, else build from window.
  if (!"label_plain" %in% names(pd)) {
    pd[, label_plain := sprintf("%d-2024", window)]
  }

  for (sh in unique(pd$shock)) {
    sub_all <- pd[shock == sh]
    # Common y-range across windows (per shock), padded by 4%.
    rng <- range(c(sub_all$lower90, sub_all$upper90), na.rm = TRUE)
    pad <- 0.04 * diff(rng)
    ylim_common <- c(rng[1] - pad, rng[2] + pad)

    sh_txt <- if (sh == "one_sd") "One-s.d. GPR shock" else "2001Q3-sized GPR shock"

    for (w in sort(unique(sub_all$window))) {
      b <- sub_all[window == w]
      lab <- b$label_plain[1]
      caption_txt <- paste0(
        lab, ". ", sh_txt, ", ", exercise_label,
        ". Posterior medians with 68% (dark) and 90% (light) credible bands. ",
        "Y-scale common to all windows of this exercise."
      )
      stem <- sprintf("pd_girf_window_%d_%s", w,
                      if (sh == "one_sd") "1sd" else "2001Q3")
      save_fig(plot_single_window(b, ylim_common, caption_txt), out_dir, stem)
      cat("  written: ", file.path(out_dir, stem), ".png/.pdf\n", sep = "")
    }
  }
  invisible(TRUE)
}

base <- file.path("output", "applications", "dralacbn", "short_sample")

cat("\n========== Single-panel short-sample figures ==========\n")

cat("\nFixed-design exercise (05):\n")
make_single_figures(
  rds_path = file.path(base, "short_sample_results.rds"),
  out_dir = file.path(base, "figures", "single"),
  exercise_label = "fixed best-BIC satellite design"
)

cat("\nPer-window BMA exercise (07, rho recalibrated):\n")
make_single_figures(
  rds_path = file.path(base, "bma", "short_sample_bma_results.rds"),
  out_dir = file.path(base, "bma", "figures", "single"),
  exercise_label = "per-window BIC screen and model averaging"
)

cat("\nPer-window BMA exercise (07, rho fixed at full-sample value):\n")
make_single_figures(
  rds_path = file.path(base, "bma_rho_fixed", "short_sample_bma_results.rds"),
  out_dir = file.path(base, "bma_rho_fixed", "figures", "single"),
  exercise_label = "per-window BMA, rho fixed at its full-sample value"
)

cat("\nDone.\n")
invisible(TRUE)
