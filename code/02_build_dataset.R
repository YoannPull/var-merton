# ----------------------------------------------------------------------------
# VAR-MERTON package — run from the ROOT of the repository
#
# Purpose:
#   Build the quarterly macroeconomic dataset used in the VAR block.
#
# Main outputs:
#   - data/processed/data_raw.csv
#   - data/processed/data_macro_full.csv
#   - data/processed/data_var_for_model.csv
#
# Notes:
#   - The S&P500 series is loaded from a frozen snapshot if available.
#   - For final replication, the snapshot should be committed to data/raw/.
#   - Quarterly dates are represented as the first month of the quarter:
#       Q1 -> YYYY-01-01
#       Q2 -> YYYY-04-01
#       Q3 -> YYYY-07-01
#       Q4 -> YYYY-10-01
#     The end-of-quarter shift, when needed for matching EBA risk data, is done
#     later in the Z-factor / satellite script.
# ----------------------------------------------------------------------------

source("code/00_setup.R")

# ─────────────────────────────────────────────────────────────────────────────
# Load packages
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(tidyquant)
  library(data.table)
  library(readxl)
  library(lubridate)
  library(zoo)
})

# ─────────────────────────────────────────────────────────────────────────────
# Paths
# ─────────────────────────────────────────────────────────────────────────────

data_dir <- if (exists("DATA_RAW_DIR")) DATA_RAW_DIR else "data/raw"
processed_dir <- if (exists("DATA_PROCESSED_DIR")) DATA_PROCESSED_DIR else "data/processed"

dir.create(processed_dir, showWarnings = FALSE, recursive = TRUE)

# ─────────────────────────────────────────────────────────────────────────────
# Utility functions
# ─────────────────────────────────────────────────────────────────────────────

assert_file_exists <- function(path) {
  if (!file.exists(path)) {
    stop("Missing required file: ", path, call. = FALSE)
  }
}

load_csv_dt <- function(filename) {
  path <- file.path(data_dir, filename)
  assert_file_exists(path)
  
  # FRED sometimes uses "." for missing values.
  dt <- fread(path, na.strings = c("", "NA", "."))
  
  if (!("observation_date" %in% names(dt))) {
    stop(
      "File ", filename, " must contain a column named 'observation_date'.",
      call. = FALSE
    )
  }
  
  dt[, DATE := as.Date(observation_date)]
  dt[, observation_date := NULL]
  
  if (anyNA(dt$DATE)) {
    stop("Date conversion produced NA values in ", filename, ".", call. = FALSE)
  }
  
  # Convert value columns to numeric if fread imported them as character.
  value_cols <- setdiff(names(dt), "DATE")
  for (col in value_cols) {
    if (is.character(dt[[col]])) {
      dt[, (col) := gsub(",", ".", get(col))]
      dt[, (col) := suppressWarnings(as.numeric(get(col)))]
    }
  }
  
  return(dt)
}

convert_qtr_to_decimal <- function(qtr) {
  qtr <- as.yearqtr(qtr)
  year <- as.integer(format(qtr, "%Y"))
  qnum <- as.integer(cycle(qtr))
  year + (qnum - 1L) * 0.25
}

decimal_quarter_to_date <- function(x) {
  year <- floor(x)
  qnum <- round((x - year) / 0.25) + 1L
  month <- 3L * (qnum - 1L) + 1L
  as.Date(sprintf("%04d-%02d-01", year, month))
}

quarterly_mean <- function(dt, date_col, value_col, out_col) {
  stopifnot(date_col %in% names(dt), value_col %in% names(dt))
  
  tmp <- copy(dt)
  tmp[, quarter := as.yearqtr(get(date_col))]
  tmp <- tmp[!is.na(get(value_col))]
  
  out <- tmp[, .(value = mean(get(value_col), na.rm = TRUE)), by = quarter]
  setnames(out, "value", out_col)
  out[, quarter := convert_qtr_to_decimal(quarter)]
  setorder(out, quarter)
  
  out
}

# ─────────────────────────────────────────────────────────────────────────────
# Monthly / daily / weekly data -> quarterly
# ─────────────────────────────────────────────────────────────────────────────

# ----------------------------------------------------------------------------
# S&P 500
# ----------------------------------------------------------------------------
# Reproducibility note:
# Yahoo Finance data can change over time. For final replication, the frozen
# snapshot data/raw/sp500_GSPC_snapshot.csv should be committed and reused.
# ----------------------------------------------------------------------------

sp500_snapshot <- file.path(data_dir, "sp500_GSPC_snapshot.csv")

if (file.exists(sp500_snapshot)) {
  dt_sp500 <- fread(sp500_snapshot)
  
  if (!("date" %in% names(dt_sp500))) {
    stop("S&P500 snapshot must contain a column named 'date'.", call. = FALSE)
  }
  if (!("adjusted" %in% names(dt_sp500))) {
    stop("S&P500 snapshot must contain a column named 'adjusted'.", call. = FALSE)
  }
  
  dt_sp500[, date := as.Date(date)]
  message("S&P500 loaded from frozen snapshot: ", sp500_snapshot)
  
} else {
  allow_downloads <- exists("ALLOW_DOWNLOADS") && isTRUE(ALLOW_DOWNLOADS)
  
  if (!allow_downloads) {
    stop(
      "Missing frozen S&P500 snapshot: ", sp500_snapshot, "\n",
      "Set ALLOW_DOWNLOADS <- TRUE for development, or add the snapshot for replication.",
      call. = FALSE
    )
  }
  
  message(
    "S&P500 snapshot not found. Downloading once from Yahoo Finance and ",
    "creating a frozen snapshot. For final replication, commit this file."
  )
  
  dt_sp500 <- tq_get("^GSPC", get = "stock.prices", from = "1900-01-01")
  dt_sp500 <- as.data.table(dt_sp500)
  
  if (nrow(dt_sp500) == 0L) {
    stop("Yahoo Finance returned an empty S&P500 dataset.", call. = FALSE)
  }
  
  fwrite(dt_sp500, sp500_snapshot)
  message("S&P500 frozen into: ", sp500_snapshot)
}

dt_sp500[, quarter := as.yearqtr(date)]
dt_sp500_quarterly <- dt_sp500[
  !is.na(adjusted),
  .(sp500 = mean(adjusted, na.rm = TRUE)),
  by = quarter
]
dt_sp500_quarterly[, quarter := convert_qtr_to_decimal(quarter)]
setorder(dt_sp500_quarterly, quarter)

# ----------------------------------------------------------------------------
# VIX
# ----------------------------------------------------------------------------

dt_vix <- load_csv_dt("VIXCLS.csv")
dt_vix_quarterly <- quarterly_mean(dt_vix, "DATE", "VIXCLS", "vix")

# ----------------------------------------------------------------------------
# Nonfarm payrolls
# ----------------------------------------------------------------------------

dt_payems <- load_csv_dt("PAYEMS.csv")
dt_payems_quarterly <- quarterly_mean(dt_payems, "DATE", "PAYEMS", "payems")

# ----------------------------------------------------------------------------
# Nonfarm business hours (HOANBS)
#
# NOTE: HOANBS is retained in the processed data as column `hours` for possible
# robustness checks. It is not used to build the baseline activity variable.
# The baseline activity variable is based on private employment (USPRIV), below.
# ----------------------------------------------------------------------------

dt_hours <- load_csv_dt("HOANBS.csv")
dt_hours_quarterly <- quarterly_mean(dt_hours, "DATE", "HOANBS", "hours")

# ----------------------------------------------------------------------------
# Private employment (USPRIV) — baseline real-activity variable
#
# Baseline activity variable for the paper. The processed column is named
# `log_hours_pc` for backward compatibility with the existing pipeline, but the
# economic content is private employment per capita:
#   100 * log(USPRIV / CNP16OV).
# In the paper/tables, describe it as private employment per capita, not hours.
# ----------------------------------------------------------------------------

dt_uspriv <- load_csv_dt("USPRIV.csv")
dt_uspriv_quarterly <- quarterly_mean(dt_uspriv, "DATE", "USPRIV", "priv_emp")

# ----------------------------------------------------------------------------
# VIX composite with Caldara historical VIX/SPVXO
# ----------------------------------------------------------------------------

vix_caldara_path <- file.path(data_dir, "vix_caldara.csv")
assert_file_exists(vix_caldara_path)

dt_vix_caldara <- fread(vix_caldara_path)

if (!all(c("quarter", "SPVXO") %in% names(dt_vix_caldara))) {
  stop("vix_caldara.csv must contain columns 'quarter' and 'SPVXO'.",
       call. = FALSE)
}

dt_vix_caldara[, quarter := as.numeric(quarter)]

dt_merged_vix <- merge(
  dt_vix_quarterly,
  dt_vix_caldara[, .(quarter, SPVXO)],
  by = "quarter",
  all = TRUE
)

all_quarters_vix <- seq(
  from = min(dt_merged_vix$quarter, na.rm = TRUE),
  to   = max(dt_merged_vix$quarter, na.rm = TRUE),
  by   = 0.25
)

dt_complete_vix <- data.table(quarter = all_quarters_vix)

dt_merged_vix <- merge(
  dt_complete_vix,
  dt_merged_vix,
  by = "quarter",
  all.x = TRUE
)

setorder(dt_merged_vix, quarter)

dt_merged_vix[, vix_composite := fifelse(!is.na(vix), vix, SPVXO)]

dt_vix_quarterly <- dt_merged_vix[, .(
  quarter,
  vix = vix_composite
)]

# ----------------------------------------------------------------------------
# CPI
# ----------------------------------------------------------------------------

dt_cpi <- load_csv_dt("CPIAUCSL.csv")
dt_cpi_quarterly <- quarterly_mean(dt_cpi, "DATE", "CPIAUCSL", "cpi")

# ----------------------------------------------------------------------------
# Population 16+
# ----------------------------------------------------------------------------

dt_pop16 <- load_csv_dt("CNP16OV.csv")
dt_pop16_quarterly <- quarterly_mean(dt_pop16, "DATE", "CNP16OV", "pop16")

# ----------------------------------------------------------------------------
# WTI oil price
# ----------------------------------------------------------------------------

dt_wti <- load_csv_dt("WTISPLC.csv")
dt_wti_quarterly <- quarterly_mean(dt_wti, "DATE", "WTISPLC", "wti")

# ----------------------------------------------------------------------------
# 2-year Treasury rate
# ----------------------------------------------------------------------------

dt_gs2 <- load_csv_dt("GS2.csv")
dt_gs2_quarterly <- quarterly_mean(dt_gs2, "DATE", "GS2", "gs2")

# ----------------------------------------------------------------------------
# 10Y-2Y Treasury spread
# ----------------------------------------------------------------------------

dt_T10Y2Y <- load_csv_dt("T10Y2Y.csv")
dt_t10Y2Y_quarterly <- quarterly_mean(dt_T10Y2Y, "DATE", "T10Y2Y", "t10Y2Y")

# ----------------------------------------------------------------------------
# 10Y-3M Treasury spread
# ----------------------------------------------------------------------------

dt_T10Y3M <- load_csv_dt("T10Y3M.csv")
dt_t10Y3M_quarterly <- quarterly_mean(dt_T10Y3M, "DATE", "T10Y3M", "t10Y3M")

# ----------------------------------------------------------------------------
# NFCI
# ----------------------------------------------------------------------------

dt_nfci <- load_csv_dt("NFCI.csv")
dt_nfci_quarterly <- quarterly_mean(dt_nfci, "DATE", "NFCI", "nfci")

# ----------------------------------------------------------------------------
# EPU index
# ----------------------------------------------------------------------------

epu_path <- file.path(data_dir, "EPU_US.xlsx")
assert_file_exists(epu_path)

dt_epu <- as.data.table(read_excel(epu_path, sheet = 1))

required_epu_cols <- c("Year", "Month", "News_Based_Policy_Uncert_Index")
if (!all(required_epu_cols %in% names(dt_epu))) {
  stop(
    "EPU_US.xlsx must contain columns: ",
    paste(required_epu_cols, collapse = ", "),
    call. = FALSE
  )
}

dt_epu <- dt_epu[
  !is.na(Year) &
    !is.na(Month) &
    !is.na(News_Based_Policy_Uncert_Index)
]

dt_epu[, Year := as.integer(Year)]
dt_epu[, Month := as.integer(Month)]

dt_epu[, date := as.Date(as.yearmon(paste(Year, Month), format = "%Y %m"))]
dt_epu[, quarter := as.yearqtr(date)]

dt_epu_quarterly <- dt_epu[
  ,
  .(epu = mean(News_Based_Policy_Uncert_Index, na.rm = TRUE)),
  by = quarter
]
dt_epu_quarterly[, quarter := convert_qtr_to_decimal(quarter)]
setorder(dt_epu_quarterly, quarter)

# ----------------------------------------------------------------------------
# GPR daily data
# ----------------------------------------------------------------------------

gpr_path <- file.path(data_dir, "data_gpr_daily.csv")
assert_file_exists(gpr_path)

dt_gpr <- fread(gpr_path)

cols_to_drop_gpr <- intersect(c("date", "event", "N10D"), names(dt_gpr))
if (length(cols_to_drop_gpr) > 0L) {
  dt_gpr[, (cols_to_drop_gpr) := NULL]
}

if (!("DAY" %in% names(dt_gpr))) {
  stop("data_gpr_daily.csv must contain a column named 'DAY'.", call. = FALSE)
}

cols_to_skip <- c("DAY")
for (col in setdiff(names(dt_gpr), cols_to_skip)) {
  if (is.character(dt_gpr[[col]])) {
    dt_gpr[, (col) := gsub(",", ".", get(col))]
  }
  dt_gpr[, (col) := suppressWarnings(as.numeric(get(col)))]
}

setnames(dt_gpr, old = "DAY", new = "date")
dt_gpr[, date := as.Date(as.character(date), format = "%Y%m%d")]

if (anyNA(dt_gpr$date)) {
  stop("Date conversion produced NA values in data_gpr_daily.csv.", call. = FALSE)
}

required_gpr_cols <- c("GPRD", "GPRD_ACT", "GPRD_THREAT")
if (!all(required_gpr_cols %in% names(dt_gpr))) {
  stop(
    "data_gpr_daily.csv must contain columns: ",
    paste(required_gpr_cols, collapse = ", "),
    call. = FALSE
  )
}

dt_gpr[, quarter := as.yearqtr(date)]

dt_gpr_quarterly <- dt_gpr[
  ,
  .(
    GPRD        = mean(GPRD, na.rm = TRUE),
    GPRD_ACT    = mean(GPRD_ACT, na.rm = TRUE),
    GPRD_THREAT = mean(GPRD_THREAT, na.rm = TRUE)
  ),
  by = quarter
]
dt_gpr_quarterly[, quarter := convert_qtr_to_decimal(quarter)]
setorder(dt_gpr_quarterly, quarter)

# ----------------------------------------------------------------------------
# Unemployment rate
# ----------------------------------------------------------------------------

dt_unemp <- load_csv_dt("UNRATE.csv")
dt_unemp_quarterly <- quarterly_mean(dt_unemp, "DATE", "UNRATE", "unrate")

# ─────────────────────────────────────────────────────────────────────────────
# Already-quarterly data
# ─────────────────────────────────────────────────────────────────────────────

# ----------------------------------------------------------------------------
# Real GDP
# ----------------------------------------------------------------------------

dt_gdp <- load_csv_dt("GDPC1.csv")
dt_gdp[, quarter := as.yearqtr(DATE)]
dt_gdp_quarterly <- dt_gdp[, .(gdp = GDPC1), by = quarter]
dt_gdp_quarterly[, quarter := convert_qtr_to_decimal(quarter)]
setorder(dt_gdp_quarterly, quarter)

# ----------------------------------------------------------------------------
# Real fixed investment
# ----------------------------------------------------------------------------

dt_inv <- load_csv_dt("GPDIC1.csv")
dt_inv[, quarter := as.yearqtr(DATE)]
dt_inv_quarterly <- dt_inv[, .(inv = GPDIC1), by = quarter]
dt_inv_quarterly[, quarter := convert_qtr_to_decimal(quarter)]
setorder(dt_inv_quarterly, quarter)

# ─────────────────────────────────────────────────────────────────────────────
# Merge quarterly data
# ─────────────────────────────────────────────────────────────────────────────

list_tables <- list(
  dt_sp500_quarterly,
  dt_vix_quarterly,
  dt_cpi_quarterly,
  dt_pop16_quarterly,
  dt_hours_quarterly,
  dt_uspriv_quarterly,
  dt_wti_quarterly,
  dt_gs2_quarterly,
  dt_nfci_quarterly,
  dt_epu_quarterly,
  dt_gdp_quarterly,
  dt_inv_quarterly,
  dt_unemp_quarterly,
  dt_gpr_quarterly,
  dt_t10Y2Y_quarterly,
  dt_t10Y3M_quarterly,
  dt_payems_quarterly
)

data_macro <- Reduce(
  function(x, y) merge(x, y, by = "quarter", all = TRUE),
  list_tables
)

setorder(data_macro, quarter)

# ─────────────────────────────────────────────────────────────────────────────
# Save untransformed dataset for comparison / diagnostics
# ─────────────────────────────────────────────────────────────────────────────

data_no_transform <- copy(data_macro)

setnames(data_no_transform, old = "quarter", new = "Date_quarter")
data_no_transform[, Date := decimal_quarter_to_date(Date_quarter)]
setcolorder(data_no_transform, c("Date_quarter", "Date"))

data_no_transform <- data_no_transform[
  Date <= as.Date("2022-01-01") &
    Date >= as.Date("1990-01-01")
]

fwrite(
  data_no_transform,
  file = file.path(processed_dir, "data_raw.csv")
)

# ─────────────────────────────────────────────────────────────────────────────
# Transformations following Caldara & Iacoviello-style VAR variables
# ─────────────────────────────────────────────────────────────────────────────

data_macro[, log_inv_pc := 100 * log(inv / pop16)]
data_macro[, log_gdp_pc := 100 * log(gdp / pop16)]
data_macro[, log_gdp    := 100 * log(gdp)]

data_macro[, log_gdp_yoy_pct :=
             100 * (log(gdp) - shift(log(gdp), 4L, type = "lag"))]

data_macro[, log_gdp_annualized_pct :=
             400 * (log(gdp) - shift(log(gdp), 1L, type = "lag"))]

data_macro[, gdp_yoy_pct :=
             100 * (log(gdp) - shift(log(gdp), 4L, type = "lag"))]

data_macro[, gdp_annualized_pct :=
             400 * (log(gdp) - shift(log(gdp), 1L, type = "lag"))]

data_macro[, log_private_pc := 100 * log(priv_emp / pop16)]

data_macro[, log_sp500_real := 100 * log(sp500 / cpi)]
data_macro[, log_oil_real   := 100 * log(wti / cpi)]

data_macro[, infl_annualized_pct :=
             400 * (log(cpi) - shift(log(cpi), 1L, type = "lag"))]

data_macro[, infl_yoy_pct :=
             100 * (log(cpi) - shift(log(cpi), 4L, type = "lag"))]

data_macro[, log_payems := 100 * log(payems)]

# GPR transformations
if (any(data_macro$GPRD <= 0, na.rm = TRUE)) {
  stop("GPRD contains non-positive values; log_GPRD cannot be computed safely.",
       call. = FALSE)
}
if (any(data_macro$GPRD_ACT <= 0, na.rm = TRUE)) {
  stop("GPRD_ACT contains non-positive values; log_GPRD_ACT cannot be computed safely.",
       call. = FALSE)
}
if (any(data_macro$GPRD_THREAT <= 0, na.rm = TRUE)) {
  stop("GPRD_THREAT contains non-positive values; log_GPRD_THREAT cannot be computed safely.",
       call. = FALSE)
}

data_macro[, log_GPRD       := 100 * log(GPRD)]
data_macro[, logdiff_GPRD   := 100 * (log(GPRD) - shift(log(GPRD), 1L, type = "lag"))]
data_macro[, diff_GPRD      := 100 * (GPRD - shift(GPRD, 1L, type = "lag"))]

data_macro[, log_GPRD_ACT     := 100 * log(GPRD_ACT)]
data_macro[, logdiff_GPRD_ACT := 100 * (log(GPRD_ACT) - shift(log(GPRD_ACT), 1L, type = "lag"))]
data_macro[, diff_GPRD_ACT    := 100 * (GPRD_ACT - shift(GPRD_ACT, 1L, type = "lag"))]

data_macro[, log_GPRD_THREAT     := 100 * log(GPRD_THREAT)]
data_macro[, logdiff_GPRD_THREAT := 100 * (log(GPRD_THREAT) - shift(log(GPRD_THREAT), 1L, type = "lag"))]
data_macro[, diff_GPRD_THREAT    := 100 * (GPRD_THREAT - shift(GPRD_THREAT, 1L, type = "lag"))]

# Add date columns
setnames(data_macro, old = "quarter", new = "Date_quarter")
data_macro[, Date := decimal_quarter_to_date(Date_quarter)]
setcolorder(data_macro, c("Date_quarter", "Date"))

# ─────────────────────────────────────────────────────────────────────────────
# Build final dataset for VAR estimation
# ─────────────────────────────────────────────────────────────────────────────

required_final_cols <- c(
  "Date_quarter",
  "Date",
  "log_inv_pc",
  "log_gdp_pc",
  "log_gdp",
  "log_gdp_yoy_pct",
  "log_gdp_annualized_pct",
  "gdp_yoy_pct",
  "gdp_annualized_pct",
  "log_private_pc",
  "log_sp500_real",
  "log_oil_real",
  "infl_annualized_pct",
  "infl_yoy_pct",
  "log_payems",
  "vix",
  "gs2",
  "t10Y2Y",
  "t10Y3M",
  "nfci",
  "epu",
  "cpi",
  "sp500",
  "gdp",
  "pop16",
  "inv",
  "hours",
  "priv_emp",
  "wti",
  "unrate",
  "payems",
  "logdiff_GPRD",
  "diff_GPRD",
  "GPRD",
  "log_GPRD",
  "GPRD_ACT",
  "log_GPRD_ACT",
  "logdiff_GPRD_ACT",
  "diff_GPRD_ACT",
  "GPRD_THREAT",
  "log_GPRD_THREAT",
  "logdiff_GPRD_THREAT",
  "diff_GPRD_THREAT"
)

missing_final_cols <- setdiff(required_final_cols, names(data_macro))
if (length(missing_final_cols) > 0L) {
  stop(
    "Missing columns in data_macro: ",
    paste(missing_final_cols, collapse = ", "),
    call. = FALSE
  )
}

data_var <- data_macro[, ..required_final_cols]

# ─────────────────────────────────────────────────────────────────────────────
# Keep the longest NA-free span for the main baseline VAR
# ─────────────────────────────────────────────────────────────────────────────

# Important:
# We do NOT define the main VAR sample using every auxiliary variable.
# Otherwise unused variables such as EPU or spreads could shorten the baseline
# sample. Robustness scripts should perform their own complete-case filtering.
cols_to_check <- c(
  "log_GPRD",
  "vix",
  "log_sp500_real",
  "log_oil_real",
  "log_private_pc",
  "log_gdp_pc"
)

missing_check_cols <- setdiff(cols_to_check, names(data_var))
if (length(missing_check_cols) > 0L) {
  stop(
    "Missing baseline complete-case columns: ",
    paste(missing_check_cols, collapse = ", "),
    call. = FALSE
  )
}

data_var[, complete_row := !Reduce(`|`, lapply(.SD, is.na)), .SDcols = cols_to_check]

data_var[, group := rleid(complete_row)]

valid_groups <- data_var[
  complete_row == TRUE,
  .N,
  by = group
][order(-N)]

if (nrow(valid_groups) == 0L) {
  stop("No complete NA-free span found in data_var.", call. = FALSE)
}

max_valid_group <- valid_groups[1, group]

data_var_clean <- data_var[
  group == max_valid_group &
    complete_row == TRUE
]

data_var_clean[, c("complete_row", "group") := NULL]

if (nrow(data_var_clean) == 0L) {
  stop("data_var_clean is empty after complete-span selection.", call. = FALSE)
}

# ─────────────────────────────────────────────────────────────────────────────
# Diagnostics / sanity checks
# ─────────────────────────────────────────────────────────────────────────────

message("Final VAR dataset span:")
message("  Start: ", min(data_var_clean$Date, na.rm = TRUE))
message("  End:   ", max(data_var_clean$Date, na.rm = TRUE))
message("  Rows:  ", nrow(data_var_clean))
message("  Cols:  ", ncol(data_var_clean))

if (anyDuplicated(data_var_clean$Date) > 0L) {
  stop("Duplicate Date values detected in data_var_clean.", call. = FALSE)
}

if (!all(diff(data_var_clean$Date_quarter) > 0)) {
  stop("Date_quarter is not strictly increasing in data_var_clean.",
       call. = FALSE)
}

# ─────────────────────────────────────────────────────────────────────────────
# Sample truncation
# ─────────────────────────────────────────────────────────────────────────────
# The estimation sample ends at DATA_END_DATE (config.R). Later observations
# are dropped here so that every downstream block uses the same cutoff; the
# trailing quarters of a fresh download are typically incomplete (e.g.,
# partial within-quarter averages of the daily GPR index).

if (exists("DATA_END_DATE")) {
  data_macro <- data_macro[as.Date(Date) <= DATA_END_DATE]
  data_var_clean <- data_var_clean[as.Date(Date) <= DATA_END_DATE]
  message("Sample truncated at DATA_END_DATE = ", format(DATA_END_DATE),
          " (last macro observation: ",
          format(max(as.Date(data_var_clean$Date))), ")")
}

# ─────────────────────────────────────────────────────────────────────────────
# Save outputs
# ─────────────────────────────────────────────────────────────────────────────

fwrite(
  data_macro,
  file = file.path(processed_dir, "data_macro_full.csv")
)

fwrite(
  data_var_clean,
  file = file.path(processed_dir, "data_var_for_model.csv")
)

message("Saved:")
message("  - ", file.path(processed_dir, "data_raw.csv"))
message("  - ", file.path(processed_dir, "data_macro_full.csv"))
message("  - ", file.path(processed_dir, "data_var_for_model.csv"))