# ----------------------------------------------------------------------------
# Run from the ROOT of the repository.
source("code/00_setup.R")
source("code/shared/_helpers_var_merton.R")
# ----------------------------------------------------------------------------

# =============================================================================
# EBA application — information-set robustness
#
# This wrapper prepares the EBA-specific Z factor and Merton parameters, then
# delegates the full transmission exercise to the shared engine used also by the
# DRALACBN application. Outputs are PNG-only for figures.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

APP_NAME <- "eba"
APP_LABEL <- "EBA"
APP_MAIN_OUT_DIR <- if (exists("DIR_EBA_APP")) {
  DIR_EBA_APP
} else {
  file.path("output", "applications", "eba", "main_specification")
}

APP_ROBUST_OUT_DIR <- if (exists("DIR_ROBUST_EBA_INFO")) {
  DIR_ROBUST_EBA_INFO
} else {
  file.path("output", "robustness", "eba", "information_set")
}

APP_INFO_SETS_TO_COMPARE <- if (exists("EBA_INFO_SETS_TO_COMPARE")) {
  EBA_INFO_SETS_TO_COMPARE
} else if (exists("DRALACBN_INFO_SETS_TO_COMPARE")) {
  # Default to the same grid so that EBA and DRALACBN are perfectly symmetric.
  DRALACBN_INFO_SETS_TO_COMPARE
} else {
  c("real_side", "real_side_lighter", "financial",
    "monetary", "uncertainty", "caldara_style")
}

APP_BASELINE_INFO_SET <- if (exists("EBA_BASELINE_INFO_SET")) {
  EBA_BASELINE_INFO_SET
} else if (exists("DELINQUENCY_BASELINE_INFO_SET") && DELINQUENCY_BASELINE_INFO_SET %in% APP_INFO_SETS_TO_COMPARE) {
  DELINQUENCY_BASELINE_INFO_SET
} else {
  APP_INFO_SETS_TO_COMPARE[1]
}

z_factor_path <- if (exists("PATH_Z_FACTOR_USA")) {
  PATH_Z_FACTOR_USA
} else {
  file.path(if (exists("DIR_ZFACTOR")) DIR_ZFACTOR else "output/applications/eba/main_specification/zfactor", "z_factor_usa.csv")
}

merton_params_path <- if (exists("PATH_MERTON_PARAMS_USA")) {
  PATH_MERTON_PARAMS_USA
} else {
  file.path(if (exists("DIR_ZFACTOR")) DIR_ZFACTOR else "output/applications/eba/main_specification/zfactor", "merton_params_usa.csv")
}

merton_params_by_rho_path <- if (exists("PATH_MERTON_PARAMS_BY_RHO_USA")) {
  PATH_MERTON_PARAMS_BY_RHO_USA
} else {
  file.path(if (exists("DIR_ZFACTOR")) DIR_ZFACTOR else "output/applications/eba/main_specification/zfactor", "merton_params_by_rho_usa.csv")
}

if (!file.exists(z_factor_path)) {
  stop("Missing EBA Z-factor file: ", z_factor_path, "\nRun code/extra/eba/03_estimate_zfactor.R first.", call. = FALSE)
}

read_merton_params <- function(path_main, path_fallback = NULL) {
  if (file.exists(path_main)) {
    mp <- fread(path_main)
    source_path <- path_main
  } else if (!is.null(path_fallback) && file.exists(path_fallback)) {
    mp <- fread(path_fallback)
    source_path <- path_fallback
  } else {
    stop("Missing EBA Merton parameter file: ", path_main, call. = FALSE)
  }
  if (!"rho_value" %in% names(mp) && "rho" %in% names(mp)) setnames(mp, "rho", "rho_value")
  if ("model_tag" %in% names(mp) && any(mp$model_tag == "rho_estimated")) mp <- mp[model_tag == "rho_estimated"]
  if (!all(c("p_ttc", "rho_value") %in% names(mp))) {
    stop("Merton parameter file must contain p_ttc and rho/rho_value: ", source_path, call. = FALSE)
  }
  list(p_ttc = as.numeric(mp$p_ttc[1]), rho = as.numeric(mp$rho_value[1]), source = source_path)
}

z_eba <- fread(z_factor_path)
if (!all(c("Date", "Z") %in% names(z_eba))) {
  stop("EBA Z-factor file must contain Date and Z columns: ", z_factor_path, call. = FALSE)
}
z_eba[, Date := as.Date(Date)]

mp <- read_merton_params(merton_params_path, merton_params_by_rho_path)
APP_Z_DT <- z_eba[, .(Date, Z)]
APP_P_TTC <- mp$p_ttc
APP_RHO <- mp$rho
APP_Z_EXPORT <- z_eba
APP_Z_EXPORT_NAME <- "zfactor_eba.csv"
APP_NOTE <- paste0("For EBA, Z and Merton--Vasicek parameters are read from ", z_factor_path, "; Merton parameters from ", mp$source, ".")

source("code/shared/_application_information_set_robustness.R")
