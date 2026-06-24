# ----------------------------------------------------------------------------
# Run from the ROOT of the repository.
source("code/00_setup.R")
source("code/shared/_helpers_var_merton.R")
# ----------------------------------------------------------------------------

# =============================================================================
# DRALACBN application — information-set robustness
#
# This wrapper prepares the DRALACBN-specific Z factor and Merton parameters,
# then delegates the full transmission exercise to the shared engine used also
# by the EBA application. Outputs are PNG-only for figures.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

APP_NAME <- "dralacbn"
APP_LABEL <- "DRALACBN"
APP_MAIN_OUT_DIR <- if (exists("DIR_DRALACBN_MAIN")) {
  DIR_DRALACBN_MAIN
} else if (exists("DIR_DRALACBN_BASELINE")) {
  DIR_DRALACBN_BASELINE
} else {
  file.path("output", "applications", "dralacbn", "main_specification")
}

APP_ROBUST_OUT_DIR <- if (exists("DIR_ROBUST_DRALACBN_INFO")) {
  DIR_ROBUST_DRALACBN_INFO
} else if (exists("DIR_DRALACBN_INFO_SET")) {
  DIR_DRALACBN_INFO_SET
} else {
  file.path("output", "robustness", "dralacbn", "information_set")
}

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

# Load the delinquency proxy and reconstruct the systematic factor Z once.
delinquency_path <- file.path(
  if (exists("DIR_RAW")) DIR_RAW else "data/raw",
  "default",
  "DRALACBN.csv"
)

if (!file.exists(delinquency_path)) {
  stop(
    "Missing delinquency series: ", delinquency_path,
    "\nDownload DRALACBN from FRED and save it there with a date column and value column.",
    call. = FALSE
  )
}

delinq_raw <- fread(delinquency_path)
date_col <- intersect(c("DATE", "Date", "date", "observation_date"), names(delinq_raw))[1]
val_col <- intersect(c("DRALACBN", "value", "VALUE", "Value"), names(delinq_raw))[1]
if (is.na(date_col)) date_col <- names(delinq_raw)[1]
if (is.na(val_col)) val_col <- setdiff(names(delinq_raw), date_col)[1]

delinq <- data.table(Date = as.Date(delinq_raw[[date_col]]), value = as.numeric(delinq_raw[[val_col]]))
delinq <- delinq[is.finite(value)]
setorder(delinq, Date)
if (exists("DATA_END_DATE")) delinq <- delinq[Date <= DATA_END_DATE]

z_out <- f_Z_estimation(delinq$value / 100)
APP_RHO <- z_out$rho
APP_P_TTC <- z_out$p_ttc
Z_raw <- z_out$Z

# Align DRALACBN dates to the quarter-start convention used by data_var_for_model.csv.
APP_Z_DT <- data.table(
  Date = as.Date(format(delinq$Date + 62, "%Y-%m-01")),
  Z = Z_raw
)

APP_Z_EXPORT <- data.table(
  Date = delinq$Date,
  delinquency = delinq$value,
  Z = Z_raw,
  rho = APP_RHO,
  p_ttc = APP_P_TTC,
  mean_Z_raw = mean(Z_raw, na.rm = TRUE),
  var_Z_raw = var(Z_raw, na.rm = TRUE)
)
APP_Z_EXPORT_NAME <- "zfactor_dralacbn.csv"
APP_NOTE <- "For DRALACBN, Z is reconstructed from the FRED delinquency proxy and kept in its original Merton--Vasicek scale."

source("code/shared/_application_information_set_robustness.R")
