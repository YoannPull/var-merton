# =============================================================================
#  code/shared/replot_var_irf_gpr_first.R
#
#  Re-renders ALL saved macro-financial IRF facet figures with the GPR panel
#  first, WITHOUT re-estimating anything. Reads the macro_girf_*.rds objects
#  saved by the pipeline and overwrites the corresponding figures (PNG).
#
#  Requires the facet-ordering patch in code/shared/_helpers_var_merton.R
#  (facets ordered by VAR variable ordering, impulse variable first).
#
#  Covered locations (skipped silently if absent):
#    output/applications/<app>/main_specification/
#    output/robustness/<app>/information_set/by_information_set/<spec>/
#  for app in {dralacbn, coralacbn, eba}.
#
#  Run from the repository root (a few seconds):
#      source("code/shared/replot_var_irf_gpr_first.R")
#
#  NB: the collated paper figures (output/paper/, figures/main/) are copies;
#  re-run code/shared/collate_paper_outputs.R or re-copy manually afterwards.
# =============================================================================

source("code/00_setup.R")
source("code/shared/_helpers_var_merton.R")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

var_labels <- c(
  "log_GPRD"               = "Geopolitical Risk (log)",
  "vix"                    = "VIX",
  "log_sp500_real"         = "Real S&P 500 (log)",
  "log_oil_real"           = "Real Oil Price (log)",
  "log_private_pc"         = "Private Employment (p.c., log)",
  "log_gdp_pc"             = "Real GDP (p.c., log)",
  "log_inv_pc"             = "Investment (p.c., log)",
  "infl_yoy_pct"           = "Inflation YoY (%)",
  "gs2"                    = "2-Year Treasury Yield",
  "nfci"                   = "NFCI",
  "epu"                    = "EPU",
  "t10Y3M"                 = "10Y-3M Term Spread",
  "t10Y2Y"                 = "10Y-2Y Term Spread",
  "gdp_yoy_pct"            = "GDP Growth YoY (%)",
  "gdp_annualized_pct"     = "GDP Growth Ann. QoQ (%)",
  "infl_annualized_pct"    = "Inflation Ann. QoQ (%)",
  "log_payems"             = "Nonfarm Payrolls (log)",
  "unrate"                 = "Unemployment Rate (%)"
)

replot_dir <- function(info_set_dir) {
  done <- 0L
  for (sh in c("1sd", "2001Q3")) {
    rds <- file.path(info_set_dir, "irf", "VAR",
                     paste0("macro_girf_", sh, ".rds"))
    if (!file.exists(rds)) next
    obj <- readRDS(rds)
    if (is.null(obj$bands)) next
    ggsave_both(
      plot_bands_facets(obj$bands, ylab_txt = "VAR GIRF",
                        var_labels = var_labels),
      file.path(info_set_dir, "figures"),
      paste0("var_irf_", sh),
      width = 8.5, height = 5.6
    )
    done <- done + 1L
  }
  done
}

apps <- c("dralacbn", "coralacbn", "eba")
total <- 0L

for (app in apps) {
  main_dir <- file.path("output", "applications", app, "main_specification")
  if (dir.exists(main_dir)) {
    n <- replot_dir(main_dir)
    if (n > 0) cat("Re-rendered ", n, " figure(s): ", main_dir, "\n", sep = "")
    total <- total + n
  }
  rob_root <- file.path("output", "robustness", app, "information_set",
                        "by_information_set")
  if (dir.exists(rob_root)) {
    for (d in list.dirs(rob_root, recursive = FALSE)) {
      n <- replot_dir(d)
      if (n > 0) cat("Re-rendered ", n, " figure(s): ", d, "\n", sep = "")
      total <- total + n
    }
  }
}

cat("\nDone: ", total, " figures re-rendered with the GPR panel first.\n", sep = "")
cat("Remember to refresh the collated/paper copies (collate_paper_outputs.R\n")
cat("or manual copy to figures/main/ and figures/robustness/).\n")

invisible(TRUE)
