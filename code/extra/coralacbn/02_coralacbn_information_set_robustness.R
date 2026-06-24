# ----------------------------------------------------------------------------
# Run from the ROOT of the repository.
source("code/00_setup.R")
source("code/shared/_helpers_var_merton.R")
# ----------------------------------------------------------------------------

# =============================================================================
# CORALACBN application — information-set robustness
#
# Mirror of code/dralacbn/02_dralacbn_information_set_robustness.R with the
# charge-off rate on all loans and leases, all U.S. commercial banks
# (FRED: CORALACBN) as the default proxy. Charge-offs are realized losses and
# are conceptually closer to a default rate than the delinquency proxy; they
# also lag delinquencies, so the PD response may peak slightly later.
#
# This wrapper prepares the CORALACBN-specific Z factor and Merton parameters,
# then delegates the full transmission exercise to the shared engine used by
# the EBA and DRALACBN applications. Outputs are strictly symmetric to the
# DRALACBN application:
#   output/applications/coralacbn/main_specification/
#   output/robustness/coralacbn/information_set/
#   output/model_averaging/coralacbn/
#
# STANDALONE: not part of run_all.R. Run from the repository root:
#   source("code/extra/coralacbn/02_coralacbn_information_set_robustness.R")
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

APP_NAME <- "coralacbn"
APP_LABEL <- "CORALACBN"

APP_MAIN_OUT_DIR <- file.path("output", "applications", "coralacbn",
                              "main_specification")
APP_ROBUST_OUT_DIR <- file.path("output", "robustness", "coralacbn",
                                "information_set")

APP_INFO_SETS_TO_COMPARE <- if (exists("DRALACBN_INFO_SETS_TO_COMPARE")) {
  DRALACBN_INFO_SETS_TO_COMPARE
} else {
  c("real_side", "real_side_lighter", "financial",
    "monetary", "uncertainty", "caldara_style")
}

APP_BASELINE_INFO_SET <- if (exists("DELINQUENCY_BASELINE_INFO_SET")) {
  DELINQUENCY_BASELINE_INFO_SET
} else {
  APP_INFO_SETS_TO_COMPARE[1]
}

# Load the charge-off proxy and reconstruct the systematic factor Z once.
chargeoff_path <- file.path(
  if (exists("DIR_RAW")) DIR_RAW else "data/raw",
  "default",
  "CORALACBN.csv"
)

if (!file.exists(chargeoff_path)) {
  stop(
    "Missing charge-off series: ", chargeoff_path,
    "\nDownload CORALACBN from FRED and save it there with a date column and value column.",
    call. = FALSE
  )
}

chargeoff_raw <- fread(chargeoff_path)
date_col <- intersect(c("DATE", "Date", "date", "observation_date"), names(chargeoff_raw))[1]
val_col <- intersect(c("CORALACBN", "value", "VALUE", "Value"), names(chargeoff_raw))[1]
if (is.na(date_col)) date_col <- names(chargeoff_raw)[1]
if (is.na(val_col)) val_col <- setdiff(names(chargeoff_raw), date_col)[1]

chargeoff <- data.table(
  Date = as.Date(chargeoff_raw[[date_col]]),
  value = as.numeric(chargeoff_raw[[val_col]])
)
chargeoff <- chargeoff[is.finite(value)]
setorder(chargeoff, Date)
if (exists("DATA_END_DATE")) chargeoff <- chargeoff[Date <= DATA_END_DATE]

z_out <- f_Z_estimation(chargeoff$value / 100)
APP_RHO <- z_out$rho
APP_P_TTC <- z_out$p_ttc
Z_raw <- z_out$Z

# Align CORALACBN dates to the quarter-start convention used by data_var_for_model.csv.
APP_Z_DT <- data.table(
  Date = as.Date(format(chargeoff$Date + 62, "%Y-%m-01")),
  Z = Z_raw
)

APP_Z_EXPORT <- data.table(
  Date = chargeoff$Date,
  chargeoff = chargeoff$value,
  Z = Z_raw,
  rho = APP_RHO,
  p_ttc = APP_P_TTC,
  mean_Z_raw = mean(Z_raw, na.rm = TRUE),
  var_Z_raw = var(Z_raw, na.rm = TRUE)
)
APP_Z_EXPORT_NAME <- "zfactor_coralacbn.csv"
APP_NOTE <- "For CORALACBN, Z is reconstructed from the FRED charge-off proxy (realized losses) and kept in its original Merton--Vasicek scale."

source("code/shared/_application_information_set_robustness.R")
