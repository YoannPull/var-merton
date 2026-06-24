# =============================================================================
#  code/shared/collate_paper_outputs.R
#
#  Gathers the final article deliverables (figures, tables, supporting data)
#  from the deep pipeline output tree into a single flat, self-describing folder
#  output/paper/, with consistent names of the form
#     <type>_<application>_<information_set>[_<shock>].<ext>
#  so they are easy to locate and drop into the paper. An INDEX.csv lists
#  everything that was collated.
#
#  Run automatically at the end of run_all.R; can also be sourced on its own
#  after a full run (it only copies existing files).
# =============================================================================

suppressPackageStartupMessages(library(data.table))

if (!exists("DIR_PAPER")) source("config.R")

paper_fig  <- DIR_PAPER_FIG
paper_tab  <- DIR_PAPER_TAB
paper_data <- DIR_PAPER_DATA

# Start from a clean tree so stale files (renamed/removed information sets,
# previous layouts) do not linger.
for (d in c(paper_fig, paper_tab, paper_data)) {
  if (dir.exists(d)) unlink(d, recursive = TRUE)
}
for (d in c(DIR_PAPER, paper_fig, paper_tab, paper_data)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

.index <- list()

# Copy src -> dest. Files are sorted into per-application sub-folders
# (figures/<app>/, tables/<app>/, data/<app>/). Records the mapping; silently
# skips if the source is missing.
.collate <- function(src, dest_dir, dest_name, kind, app, spec, shock = "") {
  if (!file.exists(src)) return(invisible(FALSE))
  dest_sub <- file.path(dest_dir, app)
  dir.create(dest_sub, recursive = TRUE, showWarnings = FALSE)
  dest <- file.path(dest_sub, dest_name)
  ok <- file.copy(src, dest, overwrite = TRUE)
  if (ok) {
    .index[[length(.index) + 1L]] <<- data.table(
      kind = kind, application = app, information_set = spec, shock = shock,
      file = dest, source = src
    )
  }
  invisible(ok)
}

# Applications and their directory conventions (mirrors the shared engine).
# Only the DRALACBN application feeds the current paper; the EBA application has
# been moved to code/extra/ (run via run_extra.R) and is no longer collated.
apps <- list(
  list(
    name = "dralacbn", label = "DRALACBN",
    main = if (exists("DIR_DRALACBN_MAIN")) DIR_DRALACBN_MAIN
           else file.path("output", "applications", "dralacbn", "main_specification"),
    robust = if (exists("DIR_ROBUST_DRALACBN_INFO")) DIR_ROBUST_DRALACBN_INFO
             else file.path("output", "robustness", "dralacbn", "information_set"),
    baseline = if (exists("DELINQUENCY_BASELINE_INFO_SET")) DELINQUENCY_BASELINE_INFO_SET else "real_side",
    specs = if (exists("DRALACBN_INFO_SETS_TO_COMPARE")) DRALACBN_INFO_SETS_TO_COMPARE
            else c("real_side", "real_side_lighter", "financial",
                   "monetary", "uncertainty", "caldara_style")
  )
)

ma_dir <- if (exists("DIR_MODEL_AVERAGING")) {
  DIR_MODEL_AVERAGING
} else {
  file.path("output", "model_averaging")
}

shock_tags <- c(one_sd = "1sd", `2001Q3` = "2001Q3")

for (app in apps) {
  for (spec in app$specs) {
    is_main <- identical(spec, app$baseline)
    d <- if (is_main) app$main else file.path(app$robust, "by_information_set", spec)
    a <- app$name
    tag <- paste0(a, "_", spec)

    # ---- Figures: mean PD GIRF, stressed-PD (VaR) GIRF, overlay (vector PDF) ----
    for (sh in names(shock_tags)) {
      st <- shock_tags[[sh]]
      .collate(file.path(d, "figures", paste0("pd_girf_", st, ".pdf")),
               paper_fig, paste0("pd_mean_", tag, "_", st, ".pdf"),
               "fig_pd_mean", a, spec, st)
      .collate(file.path(d, "irf", "PD_VaR", paste0("pd_mean_vs_var_", sh, ".pdf")),
               paper_fig, paste0("pd_mean_vs_var_", tag, "_", st, ".pdf"),
               "fig_pd_mean_vs_var", a, spec, st)
      for (atag in c("0p99", "0p999")) {
        .collate(file.path(d, "irf", "PD_VaR",
                           paste0("pd_var_girf_", sh, "_alpha", atag, ".pdf")),
                 paper_fig,
                 paste0("pd_var", atag, "_", tag, "_", st, ".pdf"),
                 "fig_pd_var", a, spec, st)
      }
      # Z-factor and macro GIRF figures (context).
      .collate(file.path(d, "figures", paste0("z_girf_", st, ".pdf")),
               paper_fig, paste0("z_girf_", tag, "_", st, ".pdf"),
               "fig_z_girf", a, spec, st)
    }

    # ---- Tables: main spec (BMA), horseshoe, HAC, OLS satellite ----
    .collate(file.path(ma_dir, a, paste0(tag, "_averaged_regression.tex")),
             paper_tab, paste0("bma_", tag, ".tex"), "tab_bma_main", a, spec)
    .collate(file.path(ma_dir, a, paste0(tag, "_averaged_regression.csv")),
             paper_data, paste0("bma_", tag, ".csv"), "tab_bma_main_csv", a, spec)
    .collate(file.path(d, "satellite", "horseshoe",
                       paste0(tag, "_horseshoe_regression.tex")),
             paper_tab, paste0("horseshoe_", tag, ".tex"), "tab_horseshoe", a, spec)
    .collate(file.path(d, "satellite", "satellite_hac_newey_west.tex"),
             paper_tab, paste0("hac_", tag, ".tex"), "tab_hac", a, spec)
    .collate(file.path(d, "satellite", "satellite_ols_table.tex"),
             paper_tab, paste0("satellite_ols_", tag, ".tex"), "tab_satellite_ols", a, spec)

    # ---- Supporting data (bands, peaks) ----
    for (sh in names(shock_tags)) {
      st <- shock_tags[[sh]]
      .collate(file.path(d, "irf", "PD", paste0("pd_girf_bands_", st, ".csv")),
               paper_data, paste0("pd_mean_bands_", tag, "_", st, ".csv"),
               "data_pd_mean_bands", a, spec, st)
    }
    .collate(file.path(d, "irf", "PD_VaR", "pd_var_response_bands.csv"),
             paper_data, paste0("pd_var_bands_", tag, ".csv"),
             "data_pd_var_bands", a, spec)
    .collate(file.path(d, "irf", "PD_VaR", "pd_var_peaks.csv"),
             paper_data, paste0("pd_var_peaks_", tag, ".csv"),
             "data_pd_var_peaks", a, spec)
  }

  # ---- Application-level cross-information-set comparison figures/tables ----
  for (sh_tag in c("1sd", "2001Q3")) {
    .collate(file.path(app$robust, paste0("pd_girf_overlay_", sh_tag, ".pdf")),
             paper_fig, paste0("pd_overlay_infosets_", app$name, "_", sh_tag, ".pdf"),
             "fig_pd_overlay_infosets", app$name, "all", sh_tag)
  }
  .collate(file.path(app$robust, "information_set_comparison.tex"),
           paper_tab, paste0("information_set_comparison_", app$name, ".tex"),
           "tab_infoset_comparison", app$name, "all")
  .collate(file.path(app$robust, "information_set_comparison.csv"),
           paper_data, paste0("information_set_comparison_", app$name, ".csv"),
           "data_infoset_comparison", app$name, "all")
}

if (length(.index) > 0L) {
  idx <- rbindlist(.index, use.names = TRUE, fill = TRUE)
  setorder(idx, application, kind, information_set, shock)
  fwrite(idx, file.path(DIR_PAPER, "INDEX.csv"))
  cat("Collated ", nrow(idx), " article deliverables into ", DIR_PAPER, "\n",
      sep = "")
  cat("  figures: ", paper_fig, "\n", sep = "")
  cat("  tables : ", paper_tab, "\n", sep = "")
  cat("  data   : ", paper_data, "\n", sep = "")
  cat("  index  : ", file.path(DIR_PAPER, "INDEX.csv"), "\n", sep = "")
} else {
  cat("collate_paper_outputs.R: no deliverables found. Run the pipeline first.\n")
}
