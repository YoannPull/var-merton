# =============================================================================
#  dependencies.R — Required packages
#
#  RECOMMENDED (this repository ships a renv.lock): restore the exact package
#  versions with
#       install.packages("renv"); renv::restore()
#  Run this script only if you are NOT using renv.
#
#  Usage:  Rscript dependencies.R   (or source("dependencies.R") in R)
# =============================================================================

pkgs <- c(
  # Data / dates
  "data.table", "dplyr", "lubridate", "zoo", "readxl",
  # Market data
  "tidyquant",
  # Econometrics / time series
  "tseries", "lmtest", "sandwich", "MASS", "car", "strucchange",
  # Regularization
  "glmnet",
  # BVAR / linear algebra
  "mniw", "matrixStats",
  # Parallelism (base R) + optional BLAS thread control
  "parallel", "RhpcBLASctl",
  # Graphics
  "ggplot2", "ggrepel", "scales"
)

to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install)) {
  message("Installing: ", paste(to_install, collapse = ", "))
  install.packages(to_install)
} else {
  message("All required packages are already installed.")
}

invisible(lapply(pkgs, function(p)
  if (!requireNamespace(p, quietly = TRUE))
    warning("Package not available: ", p)))
