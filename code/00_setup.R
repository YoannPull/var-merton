# =============================================================================
#  code/00_setup.R — Common initialization
#
#  Sourced at the top of every analysis script and by run_all.R.
#  Loads the configuration, the shared functions and the few packages used
#  almost everywhere. Script-specific packages stay loaded by each script.
# =============================================================================

# Configuration (paths + parameters).
source("config.R")

# Shared functions, if present (some installations keep helpers in R/functions.R).
if (file.exists("R/functions.R")) {
  source("R/functions.R")
}

# Shared VAR / BVAR engine (single source of truth). Sourced here so that every
# script starting with source("code/00_setup.R") gets the same definitions of
# simulate_bvar_niw(), build_companion(), max_root(), is_stable(),
# extract_c_A_list(), compute_ma_coefficients() and forecast_baseline_path().
# Scripts must NOT redefine these functions locally.
if (file.exists("code/shared/_bvar_niw_utils.R")) {
  source("code/shared/_bvar_niw_utils.R")
} else {
  stop("code/00_setup.R: missing code/shared/_bvar_niw_utils.R (shared VAR/BVAR engine).")
}


# Shared Bayesian satellite helpers.
if (file.exists("code/shared/_bayesian_satellite_utils.R")) {
  source("code/shared/_bayesian_satellite_utils.R")
} else {
  stop("code/00_setup.R: missing code/shared/_bayesian_satellite_utils.R.")
}

# Base packages present in almost every script.
suppressPackageStartupMessages({
  library(data.table)
  library(lubridate)
  library(zoo)
})

options(stringsAsFactors = FALSE)

# Reproducibility across R versions: keep the pre-R-3.6 sampling convention
# explicit when available. This matters for posterior-draw subsampling and any
# other stochastic selection based on sample().
try(
  RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection"),
  silent = TRUE
)

invisible(TRUE)
