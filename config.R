# =============================================================================
#  config.R — Single configuration entry point for the VAR-MERTON package
#
#  Sourced by run_all.R and by code/00_setup.R. This file is the single source
#  of truth for paths, seeds, VAR/GIRF settings, information sets, Merton-
#  Vasicek parameters, and figure colours.
#
#  All paths are relative to the repository root. Start R from the root folder
#  before running source("run_all.R") or Rscript run_all.R.
# =============================================================================

# --- Check that we are at the repository root --------------------------------
if (!dir.exists("data/raw")) {
  stop(
    "config.R: 'data/raw' not found.\n",
    "Start R from the VAR-MERTON repository root (the folder that contains run_all.R).\n",
    "Current working directory: ", getwd()
  )
}

# =============================================================================
# 1. Directory tree
# =============================================================================

DIR_RAW       <- "data/raw"
DIR_PROCESSED <- "data/processed"

# Backward-compatible aliases used by data-construction scripts.
DATA_RAW_DIR       <- DIR_RAW
DATA_PROCESSED_DIR <- DIR_PROCESSED

DIR_OUTPUT  <- "output"

# Clear output layout used by the paper/reproduction pipeline.
#
# Main applications:
#   output/applications/eba/main_specification/
#   output/applications/dralacbn/main_specification/
#
# Robustness exercises:
#   output/robustness/eba/information_set/
#   output/robustness/dralacbn/information_set/
#
# Shared cross-application objects:
#   output/shared/var_kernels/
#   output/shared/var_information_sets/
DIR_APPLICATIONS <- file.path(DIR_OUTPUT, "applications")
DIR_SHARED_OUTPUT <- file.path(DIR_OUTPUT, "shared")
DIR_ROBUSTNESS <- file.path(DIR_OUTPUT, "robustness")

# --- EBA main application ----------------------------------------------------
DIR_EBA_APP <- file.path(DIR_APPLICATIONS, "eba", "main_specification")
DIR_FIGURES <- file.path(DIR_EBA_APP, "figures")
DIR_TABLES  <- file.path(DIR_EBA_APP, "tables")
DIR_MODELS  <- file.path(DIR_EBA_APP, "models")
DIR_ZFACTOR <- file.path(DIR_EBA_APP, "zfactor")
DIR_BAYESIAN_SATELLITE <- file.path(DIR_EBA_APP, "bayesian_satellite")
DIR_DIAGNOSTICS <- file.path(DIR_EBA_APP, "diagnostics")

DIR_IRF        <- file.path(DIR_EBA_APP, "irf")
DIR_IRF_VAR         <- file.path(DIR_IRF, "VAR")
DIR_IRF_Z           <- file.path(DIR_IRF, "Z")
DIR_IRF_PD          <- file.path(DIR_IRF, "PD")
DIR_IRF_DIAGNOSTICS <- file.path(DIR_IRF, "diagnostics")
DIR_SCENARIO        <- file.path(DIR_EBA_APP, "scenario")

# --- DRALACBN main application ----------------------------------------------
DIR_DRALACBN_APP <- file.path(DIR_APPLICATIONS, "dralacbn")
DIR_DRALACBN_MAIN <- file.path(DIR_DRALACBN_APP, "main_specification")
DIR_DRALACBN_BASELINE <- DIR_DRALACBN_MAIN
DIR_DELINQUENCY_PROXY <- DIR_DRALACBN_MAIN

# --- Model-averaging output root --------------------------------------------
# Dedicated, clearly-labelled folder for the satellite model-averaging (BMA)
# diagnostics of each application. One sub-folder per application.
DIR_MODEL_AVERAGING          <- file.path(DIR_OUTPUT, "model_averaging")
DIR_MODEL_AVERAGING_EBA      <- file.path(DIR_MODEL_AVERAGING, "eba")
DIR_MODEL_AVERAGING_DRALACBN <- file.path(DIR_MODEL_AVERAGING, "dralacbn")

# --- Article deliverables (collated figures and tables) ----------------------
# A single, flat, self-describing folder gathering the final figures and tables
# for the paper, built at the end of the pipeline by collate_paper_outputs.R.
DIR_PAPER      <- file.path(DIR_OUTPUT, "paper")
DIR_PAPER_FIG  <- file.path(DIR_PAPER, "figures")
DIR_PAPER_TAB  <- file.path(DIR_PAPER, "tables")
DIR_PAPER_DATA <- file.path(DIR_PAPER, "data")

# --- Robustness output roots -------------------------------------------------
DIR_ROBUST_EBA      <- file.path(DIR_ROBUSTNESS, "eba")
DIR_ROBUST_EBA_INFO <- file.path(DIR_ROBUST_EBA, "information_set")
DIR_ROBUST_DRALACBN      <- file.path(DIR_ROBUSTNESS, "dralacbn")
DIR_ROBUST_DRALACBN_INFO <- file.path(DIR_ROBUST_DRALACBN, "information_set")

# --- Shared output roots -----------------------------------------------------
DIR_VAR_KERNELS     <- file.path(DIR_SHARED_OUTPUT, "var_kernels")
DIR_VAR_INFO_SETS   <- file.path(DIR_SHARED_OUTPUT, "var_information_sets")

# Backward-compatible aliases used by older helper scripts.
DIR_VAR_SPECS       <- DIR_VAR_INFO_SETS
DIR_ROBUST_ALT_SPEC <- DIR_ROBUST_EBA_INFO
DIR_DRALACBN        <- DIR_ROBUST_DRALACBN_INFO
DIR_DRALACBN_INFO_SET <- DIR_ROBUST_DRALACBN_INFO

# Figure sub-folders for the EBA main application.
DIR_FIG_MAIN        <- file.path(DIR_FIGURES, "main")
DIR_FIG_SCENARIO    <- file.path(DIR_FIGURES, "scenario")
DIR_FIG_ROBUSTNESS  <- file.path(DIR_FIGURES, "robustness")
DIR_FIG_BY_RHO      <- file.path(DIR_FIG_ROBUSTNESS, "by_rho")
DIR_FIG_COMPARE_RHO <- file.path(DIR_FIG_ROBUSTNESS, "comparison")

DIR_LOGS <- "logs"

for (d in c(
  DIR_OUTPUT, DIR_APPLICATIONS, DIR_SHARED_OUTPUT,
  DIR_EBA_APP, DIR_DRALACBN_APP, DIR_DRALACBN_MAIN,
  DIR_PROCESSED, DIR_FIGURES, DIR_TABLES, DIR_MODELS, DIR_ZFACTOR,
  DIR_BAYESIAN_SATELLITE, DIR_DIAGNOSTICS,
  DIR_IRF, DIR_IRF_VAR, DIR_IRF_Z, DIR_IRF_PD, DIR_IRF_DIAGNOSTICS, DIR_SCENARIO,
  DIR_ROBUSTNESS, DIR_ROBUST_EBA, DIR_ROBUST_EBA_INFO,
  DIR_ROBUST_DRALACBN, DIR_ROBUST_DRALACBN_INFO,
  DIR_MODEL_AVERAGING, DIR_MODEL_AVERAGING_EBA, DIR_MODEL_AVERAGING_DRALACBN,
  DIR_PAPER, DIR_PAPER_FIG, DIR_PAPER_TAB, DIR_PAPER_DATA,
  DIR_VAR_KERNELS, DIR_VAR_INFO_SETS,
  DIR_FIG_MAIN, DIR_FIG_SCENARIO, DIR_FIG_ROBUSTNESS,
  DIR_FIG_BY_RHO, DIR_FIG_COMPARE_RHO, DIR_LOGS
)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# =============================================================================
# 2. Pipeline pivot files
# =============================================================================

PATH_DATA_VAR <- file.path(DIR_PROCESSED, "data_var_for_model.csv")

PATH_BEST_MODEL               <- file.path(DIR_MODELS,  "best_model.rds") # legacy OLS-selected design, not the main estimator
PATH_BAYESIAN_SATELLITE_MAIN  <- file.path(DIR_BAYESIAN_SATELLITE, "bayesian_satellite_posterior_draws.rds")
PATH_Z_FACTOR_USA             <- file.path(DIR_ZFACTOR, "z_factor_usa.csv")
PATH_MERTON_PARAMS_USA        <- file.path(DIR_ZFACTOR, "merton_params_usa.csv")
PATH_MERTON_PARAMS_BY_RHO_USA <- file.path(DIR_ZFACTOR, "merton_params_by_rho_usa.csv")

PATH_E1_MEDIAN           <- file.path(DIR_IRF, "e1_median.csv")
PATH_GIRF_BANDS_VAR      <- file.path(DIR_IRF_VAR, "girf_bands.csv")
PATH_PSI_Z_BANDS         <- file.path(DIR_IRF_Z,  "psiZ_bands.csv")
PATH_PD_GIRF_BANDS       <- file.path(DIR_IRF_PD, "girf_pd_bands.csv")
PATH_PSI_Z_BANDS_ALL_RHO <- file.path(DIR_IRF, "psiZ_bands_all_rho.csv")
PATH_PD_BANDS_ALL_RHO    <- file.path(DIR_IRF, "girf_pd_bands_all_rho.csv")

PATH_RAW_RISK_EBA <- file.path(DIR_RAW, "data_risk_EBA.csv")

# =============================================================================
# 3. Global parameters
# =============================================================================

SEED            <- 123
SEED_ZFACTOR    <- 12345
SEED_ROBUSTNESS <- 123

ALLOW_DOWNLOADS <- FALSE

# End of the estimation sample. Observations after this date are dropped in
# 02_build_dataset.R (macro-financial block) and wherever default-rate series
# are loaded (credit block). Rationale: the trailing quarters of a fresh data
# download are typically incomplete (e.g., partial within-quarter averages of
# the daily GPR index), and the paper states a fixed sample end.
# IMPORTANT: after changing this date, delete output/shared/var_kernels/*.rds
# so that the BVAR kernels are re-estimated on the new sample.
DATA_END_DATE <- as.Date("2024-12-31")

# Z satellite selection for the main EBA application.
NB_LAGS                <- 4
MAX_VARIABLES_IN_MODEL <- 4
P_VALUE_THRESHOLD      <- 0.05

# BVAR + GIRF settings.
P_LAGS       <- 2
NREP         <- 30000
M_TARGET     <- 10000L
N_DRAWS_BAYESIAN_SATELLITE <- M_TARGET
INCLUDE_SATELLITE_RESIDUAL_UNCERTAINTY <- TRUE
HORIZON      <- 12
IMPULSE_NAME <- "log_GPRD"
IMPULSE_IX   <- 1

# Merton-Vasicek baseline calibration.
PD_UNCOND <- 0.0079
RHO_ASSET <- 0.024

# =============================================================================
# 4. Main EBA information set
# =============================================================================
# The main EBA pipeline uses VAR_BASELINE. The DRALACBN baseline is controlled
# separately by DELINQUENCY_BASELINE_INFO_SET.
#
# The activity variable keeps the legacy column name `log_hours_pc` for code
# compatibility. In the paper, label it according to the actual construction in
# 02_build_dataset.R: private hours per capita if hours/pop16 is used, or private
# employment per capita if priv_emp/pop16 is used.

VAR_BASELINE_ACTIVITY_VAR   <- "log_private_pc"
VAR_BASELINE_ACTIVITY_LABEL <- "Private activity per capita"

VAR_BASELINE <- c(
  "log_GPRD",
  "log_inv_pc",
  "log_gdp_pc",
  "log_private_pc",
  "log_oil_real",
  "infl_yoy_pct"
)

# =============================================================================
# 5. Retained VAR information sets for robustness exercises
# =============================================================================
# Each vector lists variables IN ADDITION to the impulse variable log_GPRD.
# The impulse variable is added first by the estimation helpers.

VAR_SPECS <- list(
  # Baseline real-economy channel.
  real_side = c(
    "log_inv_pc",
    "log_gdp_pc",
    "log_private_pc",
    "log_oil_real",
    "infl_yoy_pct"
  ),
  # Leaner real-economy channel: a single activity measure (investment, the
  # identified driver) to avoid stacking collinear activity variables.
  real_side_lighter = c(
    "log_inv_pc",
    "log_oil_real",
    "infl_yoy_pct"
  ),
  # Financial-markets channel: volatility, equity, financial conditions, rate.
  financial = c(
    "vix",
    "log_sp500_real",
    "nfci",
    "gs2"
  ),
  # Monetary / rates channel: policy rate, term-spread slope, inflation, output.
  monetary = c(
    "gs2",
    "t10Y3M",
    "infl_yoy_pct",
    "log_gdp_pc"
  ),
  # Discriminant validity: does GPR matter beyond generic uncertainty (EPU/VIX)?
  uncertainty = c(
    "epu",
    "vix",
    "log_gdp_pc",
    "log_sp500_real"
  ),
  # Caldara-Iacoviello-style anchor for comparison with the literature.
  caldara_style = c(
    "vix",
    "log_inv_pc",
    "log_private_pc",
    "log_sp500_real",
    "log_oil_real",
    "gs2",
    "nfci"
  )
)

INFO_SET_LABELS <- c(
  real_side         = "Real-side",
  real_side_lighter = "Real-side (lighter)",
  financial         = "Financial",
  monetary          = "Monetary",
  uncertainty       = "Uncertainty",
  caldara_style     = "Caldara-style"
)

# Main information set and comparison grid. The EBA and DRALACBN applications
# use the same names so that main-specification and robustness outputs are
# strictly symmetric.
EBA_BASELINE_INFO_SET <- "real_side"
EBA_INFO_SETS_TO_COMPARE <- c(
  "real_side",
  "real_side_lighter",
  "financial",
  "monetary",
  "uncertainty",
  "caldara_style"
)

# Backward-compatible alias for older EBA robustness scripts.
EBA_ALT_INFORMATION_SET <- "real_side"

# DRALACBN application: baseline and comparison set.
DELINQUENCY_BASELINE_INFO_SET <- "real_side"
DELINQUENCY_BASELINE_SPEC <- DELINQUENCY_BASELINE_INFO_SET  # backward-compatible alias

DRALACBN_INFO_SETS_TO_COMPARE <- c(
  "real_side",
  "real_side_lighter",
  "financial",
  "monetary",
  "uncertainty",
  "caldara_style"
)
DRALACBN_SPECS_TO_COMPARE <- DRALACBN_INFO_SETS_TO_COMPARE  # backward-compatible alias

CALDARA_STYLE_INFO_SET <- "caldara_style"
CALDARA_CLOSE_SPEC <- CALDARA_STYLE_INFO_SET  # backward-compatible alias

# Robustness satellite-selection parameters.
ALT_SPEC_NB_LAGS_SAT            <- 4
ALT_SPEC_MAX_VARIABLES_IN_MODEL <- 6
ALT_SPEC_P_VALUE_THRESHOLD      <- 0.10
PRESELECT_BY_AIC_TOPN           <- 1000L

# =============================================================================
# 5b. Satellite model averaging (all VAR information sets)
# =============================================================================
# To reduce model-selection uncertainty, the satellite regression of EVERY VAR
# information set (the main specification AND the robustness specifications) is
# replaced by an Akaike-weight average over its best specifications. The averaged
# (mixture) satellite feeds the GIRF -> Z -> PD pipeline, so the posterior of
# Delta PD integrates over which regressors enter the satellite. Set to FALSE to
# fall back to the single OOS/AIC-selected satellite.
USE_MODEL_AVERAGING_MAIN_SPEC <- TRUE

# Weighting scheme for the average: "aic" (Akaike weights, recommended),
# "bic" (Schwarz weights, more parsimonious) or "equal".
BMA_WEIGHT <- "bic"

# Rule defining the set of top models to average:
#   "cum_weight" keep best models until cumulative weight >= BMA_CUM_THRESHOLD
#   "top_n"      keep the BMA_TOP_N best models
#   "delta"      keep all models with delta-IC <= BMA_DELTA_MAX
BMA_RULE          <- "cum_weight"
BMA_CUM_THRESHOLD <- 0.95
BMA_TOP_N         <- 10L
BMA_DELTA_MAX     <- 2

# Out-of-sample evaluation. The EBA corporate default-rate sample is short, so
# its regression is estimated on the full sample with no OOS scoring. DRALACBN
# keeps its OOS diagnostics.
EBA_DROP_OOS <- TRUE

# Ridge/Lasso shrinkage robustness of the satellite (Z) coefficients. Applied in
# the shared engine to every information set of both applications. Ridge (alpha=0)
# probes coefficient stability under collinearity; Lasso (alpha=1) probes the
# robustness of the selected set. Elastic net is an optional compromise.
RUN_SATELLITE_REGULARIZATION <- TRUE
RUN_RIDGE        <- TRUE
RUN_LASSO        <- TRUE
RUN_ELASTIC_NET  <- FALSE
USE_LAMBDA       <- "lambda.1se"

# Horseshoe satellite robustness. A single Bayesian shrinkage regression on the
# full lagged design, propagated through the GIRF -> Z -> PD pipeline and
# compared with the model-averaged main specification. Outputs in
# <info_set_dir>/satellite/horseshoe/. The main spec remains the BMA satellite.
RUN_HORSESHOE_ROBUSTNESS <- TRUE
HORSESHOE_BURNIN         <- 1000L

# Value-at-Risk (stressed / downturn PD) response. Evaluates the tail of the
# systematic-factor distribution (Vasicek loss VaR) per posterior draw, so the
# reported bands are the credible bands AROUND the VaR. alpha indexes risk;
# the bands index estimation uncertainty. Outputs in <info_set_dir>/irf/PD_VaR/.
RUN_PD_VAR    <- TRUE
PD_VAR_ALPHAS <- c(0.99, 0.999)

# =============================================================================
# 6. Figure colours
# =============================================================================

SQUARE_BURGUNDY  <- "#7A1E3A"
SQUARE_BURGUNDY2 <- "#6F1732"
SQUARE_ROSE      <- "#F3D6E3"
SQUARE_ROSE2     <- "#AB4A7D"
SQUARE_DARK      <- "#1A1A1A"
SQUARE_GREY      <- "#6B6B6B"

invisible(TRUE)
