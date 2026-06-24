# ----------------------------------------------------------------------------
# Run from the ROOT of the repository.
source("code/00_setup.R")
source("code/shared/_helpers_var_merton.R")
# ----------------------------------------------------------------------------

# =============================================================================
# Documentation — Macro-financial variable definitions and transformations
#
# Emits the reference table documenting every macro-financial variable code used
# across the VAR information sets, its definition/source, and its transformation.
#
# This script is documentation-only: it does not estimate any model and does not
# depend on the application-specific default proxy (EBA or DRALACBN).
#
# Outputs:
#   output/shared/variable_definitions/variable_definitions.csv
#   output/shared/variable_definitions/variable_definitions.tex
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------------------------------------------------------
# Output directory
# -----------------------------------------------------------------------------

out_dir <- "output/shared/variable_definitions"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Variable dictionary
# -----------------------------------------------------------------------------

def_dt <- rbindlist(
  list(
    # -------------------------------------------------------------------------
    # Baseline VAR variables
    # -------------------------------------------------------------------------
    data.table(
      group = "Baseline VAR",
      code = "log_GPRD",
      definition = "Log geopolitical risk index (Caldara and Iacoviello, 2022).",
      transform = "Log level; quarterly average."
    ),
    data.table(
      group = "Baseline VAR",
      code = "vix",
      definition = "CBOE Volatility Index (FRED: VIXCLS).",
      transform = "Level; quarterly average (index)."
    ),
    data.table(
      group = "Baseline VAR",
      code = "log_sp500_real",
      definition = "Real S\\&P 500 index (Yahoo Finance: adjusted close, deflated by CPIAUCSL).",
      transform = "Log real level."
    ),
    data.table(
      group = "Baseline VAR",
      code = "log_oil_real",
      definition = "Real WTI spot oil price (FRED: WTISPLC, deflated by CPIAUCSL).",
      transform = "Log real level."
    ),
    data.table(
      group = "Baseline VAR",
      code = "log_private_pc",
      definition = "Private employment per capita (FRED: USPRIV divided by CNP16OV). The label retains the historical name \\texttt{log\\_hours\\_pc} for backward compatibility.",
      transform = "Log; per capita."
    ),
    data.table(
      group = "Baseline VAR",
      code = "log_gdp_pc",
      definition = "Real GDP per capita (FRED: GDPC1 divided by CNP16OV).",
      transform = "Log; per capita."
    ),
    
    # -------------------------------------------------------------------------
    # Robustness variables
    # -------------------------------------------------------------------------
    data.table(
      group = "Robustness variables",
      code = "nfci",
      definition = "National Financial Conditions Index (Chicago Fed / FRED: NFCI).",
      transform = "Level; values above zero indicate tighter-than-average financial conditions."
    ),
    data.table(
      group = "Robustness variables",
      code = "gs2",
      definition = "U.S. Treasury 2-year constant maturity yield (FRED: GS2).",
      transform = "Level; quarterly average (\\%)."
    ),
    data.table(
      group = "Robustness variables",
      code = "t10Y2Y",
      definition = "10-year minus 2-year Treasury spread (FRED: T10Y2Y).",
      transform = "Level; quarterly average (percentage points)."
    ),
    data.table(
      group = "Robustness variables",
      code = "t10Y3M",
      definition = "10-year minus 3-month Treasury spread (FRED: T10Y3M).",
      transform = "Level; quarterly average (percentage points)."
    ),
    data.table(
      group = "Robustness variables",
      code = "infl_yoy_pct",
      definition = "CPI year-over-year inflation (FRED: CPIAUCSL).",
      transform = "$100[\\log(CPI_t)-\\log(CPI_{t-4})]$."
    ),
    data.table(
      group = "Robustness variables",
      code = "infl_annualized_pct",
      definition = "CPI quarter-on-quarter inflation, annualized (FRED: CPIAUCSL).",
      transform = "$400[\\log(CPI_t)-\\log(CPI_{t-1})]$."
    ),
    data.table(
      group = "Robustness variables",
      code = "gdp_yoy_pct",
      definition = "Real GDP year-over-year growth (FRED: GDPC1).",
      transform = "$100[\\log(GDP_t)-\\log(GDP_{t-4})]$."
    ),
    data.table(
      group = "Robustness variables",
      code = "gdp_annualized_pct",
      definition = "Real GDP quarter-on-quarter growth, annualized (FRED: GDPC1).",
      transform = "$400[\\log(GDP_t)-\\log(GDP_{t-1})]$."
    ),
    data.table(
      group = "Robustness variables",
      code = "log_inv_pc",
      definition = "Real private fixed investment per capita (FRED: GPDIC1 divided by CNP16OV).",
      transform = "Log; per capita."
    ),
    data.table(
      group = "Robustness variables",
      code = "epu",
      definition = "Economic Policy Uncertainty index (Baker et al., 2016).",
      transform = "Level; quarterly average (index)."
    ),
    data.table(
      group = "Robustness variables",
      code = "log_payems",
      definition = "Nonfarm payroll employment (FRED: PAYEMS).",
      transform = "Log level."
    ),
    data.table(
      group = "Robustness variables",
      code = "unrate",
      definition = "Civilian unemployment rate (FRED: UNRATE).",
      transform = "Level; quarterly average (\\%)."
    )
  ),
  use.names = TRUE,
  fill = TRUE
)

# -----------------------------------------------------------------------------
# Basic validation
# -----------------------------------------------------------------------------

required_cols <- c("group", "code", "definition", "transform")
missing_cols <- setdiff(required_cols, names(def_dt))

if (length(missing_cols) > 0) {
  stop(
    "variable_definitions.R: missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

if (anyDuplicated(def_dt$code)) {
  duplicated_codes <- unique(def_dt$code[duplicated(def_dt$code)])
  stop(
    "variable_definitions.R: duplicated variable code(s): ",
    paste(duplicated_codes, collapse = ", ")
  )
}

setcolorder(def_dt, required_cols)
setorder(def_dt, group, code)

# -----------------------------------------------------------------------------
# Write CSV
# -----------------------------------------------------------------------------

csv_path <- file.path(out_dir, "variable_definitions.csv")
fwrite(def_dt, csv_path)

# -----------------------------------------------------------------------------
# Write LaTeX
# -----------------------------------------------------------------------------

tex_path <- file.path(out_dir, "variable_definitions.tex")

if (!exists("write_variable_definitions_tex", mode = "function")) {
  stop(
    "variable_definitions.R: function write_variable_definitions_tex() not found. ",
    "Please check code/shared/_helpers_var_merton.R."
  )
}

write_variable_definitions_tex(
  def_dt,
  tex_path,
  caption = "Macro-financial variables, definitions and transformations",
  label = "tab:variable_definitions"
)

# -----------------------------------------------------------------------------
# Completion message
# -----------------------------------------------------------------------------

cat("\nVariable-definition outputs written to: ", out_dir, "\n", sep = "")
cat("  - ", csv_path, "\n", sep = "")
cat("  - ", tex_path, "\n", sep = "")
