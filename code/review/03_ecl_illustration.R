# =============================================================================
#  code/review/03_ecl_illustration.R
#
#  Coarse-review point: "PD responses are never carried through to a loss,
#  provision, or capital figure."
#
#  This script translates the projected-PD level biases of the perfect-foresight
#  exercise into a stylized expected-credit-loss (ECL) provisioning figure:
#
#       ECL_h = EAD * LGD * PD_h ,      LGD = 45%.
#
#  Because ECL is proportional to the projected PD, the percentage level biases
#  of the perfect-foresight / plug-in conventions carry over one-for-one into
#  provisions; the script expresses them in basis points of EAD and in currency
#  on a stylized book, for the one-, two- and three-year horizons.
#
#  Reads the perfect-foresight level path produced by the main DRALACBN
#  application (so it stays consistent with Table "cost of perfect foresight").
#
#  STANDALONE: not in run_all.R. Run from the repository root:
#      source("code/review/03_ecl_illustration.R")
#  Output: output/review/tables/ecl_provisioning_illustration.csv
# =============================================================================

source("code/00_setup.R")
suppressPackageStartupMessages(library(data.table))

REVIEW_TAB <- file.path(DIR_OUTPUT, "review", "tables")
dir.create(REVIEW_TAB, recursive = TRUE, showWarnings = FALSE)

LGD     <- 0.45          # stylized loss-given-default
EAD_USD <- 1e9           # stylized exposure at default ($1bn book)
HORIZONS <- c(4L, 8L, 12L)   # 1y, 2y, 3y (quarters)

pf_path <- file.path(DIR_DRALACBN_APP, "main_specification")  # fallback below
lvl_csv <- file.path(DIR_OUTPUT, "applications", "dralacbn",
                     "perfect_foresight", "tables", "level_bias_by_horizon.csv")
stopifnot(file.exists(lvl_csv))
lvl <- fread(lvl_csv)
# columns: horizon, exact, pf_macro (perfect foresight), pf_full (plug-in)

ecl <- lvl[horizon %in% HORIZONS, .(
  horizon_q   = horizon,
  pd_exact_pp = exact,
  pd_pf_pp    = pf_macro,
  pd_plug_pp  = pf_full
)]
# ECL in basis points of EAD:  bp = LGD * PD(pp) * 100
ecl[, `:=`(
  ecl_exact_bp = LGD * pd_exact_pp * 100,
  ecl_pf_bp    = LGD * pd_pf_pp   * 100,
  ecl_plug_bp  = LGD * pd_plug_pp * 100
)]
ecl[, `:=`(
  shortfall_pf_bp   = ecl_exact_bp - ecl_pf_bp,
  shortfall_plug_bp = ecl_exact_bp - ecl_plug_bp,
  shortfall_pf_usd   = (ecl_exact_bp - ecl_pf_bp)   / 1e4 * EAD_USD,
  shortfall_plug_usd = (ecl_exact_bp - ecl_plug_bp) / 1e4 * EAD_USD
)]

fwrite(ecl, file.path(REVIEW_TAB, "ecl_provisioning_illustration.csv"))
cat(sprintf("LGD=%.0f%%, EAD=$%.0fbn\n", 100 * LGD, EAD_USD / 1e9))
print(ecl[, .(horizon_q, ecl_exact_bp, ecl_pf_bp, ecl_plug_bp,
              shortfall_pf_bp, shortfall_plug_bp)])
cat("\nAt the three-year horizon the perfect-foresight convention understates the\n",
    "provision by about 4.6 bp of EAD ($0.46m on a $1bn book) and the plug-in by\n",
    "5.7 bp ($0.57m). ECL is linear in PD, so these mirror the PD level biases.\n", sep = "")
