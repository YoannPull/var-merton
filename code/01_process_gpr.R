# ----------------------------------------------------------------------------
# VAR-MERTON package — run from the ROOT of the repository (see README.md / run_all.R).
# Loads config.R (paths + parameters) and R/functions.R (shared functions).
source("code/00_setup.R")
# ----------------------------------------------------------------------------

# ==================== LOADING & GPR STATIONARITY ==================== #

library(data.table)
library(lubridate)
library(tseries)
library(zoo)


# Convert yearqtr to decimal
convert_qtr_to_decimal <- function(qtr) {
  year <- floor(as.numeric(format(qtr, "%Y")))
  qnum <- as.numeric(cycle(qtr))  # returns 1 to 4
  return(year + (qnum - 1) * 0.25)
}

# 1. Load the GPR file (comma decimal separator)
df_gpr_daily <- read.csv2("data/raw/data_gpr_daily.csv", dec = ",",
                          stringsAsFactors = FALSE)
df_gpr_daily$date <- NULL

# 2. Convert numeric columns
cols_to_skip <- c("DAY", "event")
df_gpr_daily[] <- lapply(names(df_gpr_daily), function(col) {
  x <- df_gpr_daily[[col]]
  if (col %in% cols_to_skip) return(x)
  if (is.character(x)) x <- gsub(",", ".", x)
  suppressWarnings(as.numeric(x))
})
names(df_gpr_daily)[names(df_gpr_daily) == "DAY"] <- "Date"
df_gpr_daily$Date <- as.Date(as.character(df_gpr_daily$Date), format = "%Y%m%d")

# 3. Add quarters
setDT(df_gpr_daily)
df_gpr_daily[, `:=`(
  Year = year(Date),
  Quarter = quarter(Date)
)]

df_gpr_daily[, `:=`(
  Date_quarter = as.Date(paste0(Year, "-", (Quarter - 1) * 3 + 1, "-01"))
)]

df_gpr_daily <- df_gpr_daily[Date >= as.Date("1990-01-01")]

# 4. Compute quarterly means
cols_num <- setdiff(names(df_gpr_daily), c("Date", "event", "Year", 
                                           "Quarter", "Date_quarter"))
gpr_quarterly <- df_gpr_daily[, lapply(.SD, mean, na.rm = TRUE),
                              by = Date_quarter, .SDcols = cols_num]


# 5. List of GPR columns
gpr_columns <- grep("(?i)^gpr", names(gpr_quarterly), value = TRUE)

# 6. Stationarity test function
test_stationarity <- function(x) {
  ts_clean <- na.omit(x)
  adf <- tryCatch(adf.test(ts_clean, k = 2, alternative = "stationary"), error = function(e) NULL)
  kpss <- tryCatch(kpss.test(ts_clean, null = "Level"), error = function(e) NULL)
  
  list(
    ADF_Stat = if (!is.null(adf)) round(adf$statistic, 3) else NA,
    ADF_p = if (!is.null(adf)) round(adf$p.value, 3) else NA,
    KPSS_Stat = if (!is.null(kpss)) round(kpss$statistic, 3) else NA,
    KPSS_p = if (!is.null(kpss)) round(kpss$p.value, 3) else NA,
    Conclusion = if (!is.null(adf) && !is.null(kpss)) {
      if (adf$p.value < 0.05 && kpss$p.value >= 0.05) "Stationnaire"
      else if (adf$p.value >= 0.05 && kpss$p.value < 0.05) "Non stationnaire"
      else "Ambigu"
    } else {
      "Erreur"
    }
  )
}

# 7. Apply the tests to all GPR columns
results_list <- lapply(gpr_columns, function(col) {
  res <- test_stationarity(gpr_quarterly[[col]])
  c(Variable = col, res)
})
results_dt <- rbindlist(lapply(results_list, as.data.table), fill = TRUE)

# 8. Affichage
cat("=== Summary of stationarity tests on the GPR indices ===\n\n")
print(results_dt)


# ==================== GPRD STATIONARITY (diff & log-diff) ==================== #


# --- 2. Transformations GPR ---
gpr_quarterly[, log_GPRD := 100 * log(GPRD)]
gpr_quarterly[, logdiff_GPRD := 100 * (log_GPRD - shift(log_GPRD))]
gpr_quarterly[, diff_GPRD := 100 * (GPRD - shift(GPRD))]

gpr_quarterly[, log_GPRD_ACT := 100 * log(GPRD_ACT)]
gpr_quarterly[, logdiff_GPRD_ACT := 100 * (log_GPRD_ACT - shift(log_GPRD_ACT))]
gpr_quarterly[, diff_GPRD_ACT := 100 * (GPRD_ACT - shift(GPRD_ACT))]

gpr_quarterly[, log_GPRD_THREAT := 100 * log(GPRD_THREAT)]
gpr_quarterly[, logdiff_GPRD_THREAT := 100 * (log_GPRD_THREAT - shift(log_GPRD_THREAT))]
gpr_quarterly[, diff_GPRD_THREAT := 100 * (GPRD_THREAT - shift(GPRD_THREAT))]


# --- 3. Stationarity test function ---
test_stationarity <- function(ts_vector, name) {
  ts_clean <- na.omit(ts_vector)
  adf <- tryCatch(adf.test(ts_clean, k = 2, alternative = "stationary"), error = function(e) NULL)
  kpss <- tryCatch(kpss.test(ts_clean, null = "Level"), error = function(e) NULL)
  
  adf_stat <- if (!is.null(adf)) round(adf$statistic, 3) else NA
  adf_p <- if (!is.null(adf)) round(adf$p.value, 3) else NA
  kpss_stat <- if (!is.null(kpss)) round(kpss$statistic, 3) else NA
  kpss_p <- if (!is.null(kpss)) round(kpss$p.value, 3) else NA
  
  conclusion <- if (!is.null(adf) && !is.null(kpss)) {
    if (adf$p.value < 0.05 && kpss$p.value >= 0.05) "Stationnaire"
    else if (adf$p.value >= 0.05 && kpss$p.value < 0.05) "Non stationnaire"
    else "Ambigu"
  } else {
    "Erreur"
  }
  
  list(
    Variable = name,
    ADF_Stat = adf_stat,
    ADF_p = adf_p,
    KPSS_Stat = kpss_stat,
    KPSS_p = kpss_p,
    Conclusion = conclusion
  )
}

# --- 4. Apply to the GPR series ---
gpr_vars <- c("log_GPRD","logdiff_GPRD", "diff_GPRD", "GPRD","GPRD_ACT", "log_GPRD_ACT",
              "logdiff_GPRD_ACT","diff_GPRD_ACT","GPRD_THREAT",
              "log_GPRD_THREAT","logdiff_GPRD_THREAT","diff_GPRD_THREAT")
gpr_results <- lapply(gpr_vars, 
                      function(v) test_stationarity(gpr_quarterly[[v]], v))
gpr_results_dt <- rbindlist(gpr_results)

# --- 5. Print the summary ---
cat("\n=== Summary of stationarity tests - GPR indices ===\n\n")
print(gpr_results_dt)

# --- 6. Add quarterly date columns ---
gpr_quarterly$Date <- gpr_quarterly$Date_quarter
gpr_quarterly[, `:=`(
  Annee = year(Date),
  Mois = month(Date)
)]
gpr_quarterly[, `:=`(
  Date_quarter = Annee + (floor((Mois - 1) / 3) + 1) * 0.25
)]
gpr_quarterly[, Date_quarter := as.yearqtr(Date_quarter) - 0.25]
gpr_quarterly[, Date_quarter := convert_qtr_to_decimal(Date_quarter)]
# --- 7. Select the final columns for modeling ---
gpr_quarterly <- gpr_quarterly[, .(Date, Date_quarter, log_GPRD, logdiff_GPRD,
                                   diff_GPRD, GPRD,
                                   GPRD_ACT, log_GPRD_ACT,logdiff_GPRD_ACT,
                                   diff_GPRD_ACT,GPRD_THREAT,
                                   log_GPRD_THREAT,logdiff_GPRD_THREAT,
                                   diff_GPRD_THREAT
)]



fwrite(gpr_quarterly, file = "data/processed/data_gpr_quarterly.csv")