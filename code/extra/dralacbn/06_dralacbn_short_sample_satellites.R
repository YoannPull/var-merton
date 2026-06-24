# =============================================================================
#  code/extra/dralacbn/06_dralacbn_short_sample_satellites.R
#
#  Satellite regression outputs for each default window of the short-sample
#  exercise (05). For every window (1985-, 1995-, 2005-, 2015-2024):
#    - OLS estimates of the fixed best-BIC design with Newey-West HAC SEs,
#    - conjugate Bayesian posterior summaries (Jeffreys prior),
#    - window-specific Merton parameters (p_ttc, rho).
#  Plus a combined side-by-side LaTeX table (coefficient stability across
#  default samples) and combined CSVs.
#
#  Fast: no GIRF computation, regressions only (a few seconds).
#
#  STANDALONE: modifies nothing in the existing pipeline; not in run_all.R.
#  Run from the repository root:
#      source("code/extra/dralacbn/06_dralacbn_short_sample_satellites.R")
#
#  Outputs: output/applications/dralacbn/short_sample/satellites/
# =============================================================================

source("code/00_setup.R")
source("code/shared/_helpers_var_merton.R")

suppressPackageStartupMessages({
  library(data.table)
})

SS_INFO_SET    <- if (exists("DELINQUENCY_BASELINE_INFO_SET")) DELINQUENCY_BASELINE_INFO_SET else "real_side"
SS_START_YEARS <- c(1985, 1995, 2005, 2015)

NB_LAGS_SAT <- if (exists("ALT_SPEC_NB_LAGS_SAT")) ALT_SPEC_NB_LAGS_SAT else 4
seed_sat    <- if (exists("SEED_ZFACTOR")) SEED_ZFACTOR else 12345
n_draws_sat <- if (exists("M_TARGET")) M_TARGET else 10000L
impulse_name <- if (exists("IMPULSE_NAME")) IMPULSE_NAME else "log_GPRD"

out_dir <- file.path("output", "applications", "dralacbn",
                     "short_sample", "satellites")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- Default proxy ------------------------------------------------------------
delinquency_path <- file.path(
  if (exists("DIR_RAW")) DIR_RAW else "data/raw", "default", "DRALACBN.csv"
)
delinq_raw <- fread(delinquency_path)
date_col <- intersect(c("DATE", "Date", "date", "observation_date"), names(delinq_raw))[1]
val_col  <- intersect(c("DRALACBN", "value", "VALUE", "Value"), names(delinq_raw))[1]
if (is.na(date_col)) date_col <- names(delinq_raw)[1]
if (is.na(val_col))  val_col  <- setdiff(names(delinq_raw), date_col)[1]
delinq <- data.table(Date = as.Date(delinq_raw[[date_col]]),
                     value = as.numeric(delinq_raw[[val_col]]))
delinq <- delinq[is.finite(value)]
setorder(delinq, Date)
if (exists("DATA_END_DATE")) delinq <- delinq[Date <= DATA_END_DATE]

DT_raw <- safe_fread(if (exists("PATH_DATA_VAR")) PATH_DATA_VAR else "data/processed/data_var_for_model.csv")
DT_raw[, Date := as.Date(Date)]

# --- Fixed design from the full-sample best-BIC satellite ---------------------
best_sat_path <- file.path("output", "applications", "dralacbn",
                           "main_specification", "satellite", "best_satellite.rds")
if (!file.exists(best_sat_path)) {
  stop("Missing ", best_sat_path, call. = FALSE)
}
best_sat <- readRDS(best_sat_path)
sat_formula <- formula(best_sat)
sat_terms <- setdiff(names(coef(best_sat)), "(Intercept)")
all_terms <- c("(Intercept)", sat_terms)

vars_sat <- setdiff(
  c(impulse_name, VAR_SPECS[[SS_INFO_SET]]),
  impulse_name
)

have_hac <- requireNamespace("sandwich", quietly = TRUE) &&
  requireNamespace("lmtest", quietly = TRUE)
if (!have_hac) {
  message("sandwich/lmtest unavailable: classical OLS SEs will be reported.")
}

# --- One window ---------------------------------------------------------------
run_window <- function(start_year) {
  d_w <- delinq[Date >= as.Date(paste0(start_year, "-01-01"))]
  z_out <- f_Z_estimation(d_w$value / 100)

  Z_DT_w <- data.table(
    Date = as.Date(format(d_w$Date + 62, "%Y-%m-01")),
    Z = z_out$Z
  )
  sat_df_w <- prepare_lagged_satellite_data(
    z_input = Z_DT_w, DT_raw = DT_raw, vars_sat = vars_sat, nb_lags = NB_LAGS_SAT
  )

  fit <- lm(sat_formula, data = sat_df_w)
  sm <- summary(fit)

  ct <- if (have_hac) {
    lmtest::coeftest(fit, vcov. = sandwich::NeweyWest(fit, lag = NB_LAGS_SAT))
  } else {
    sm$coefficients
  }
  ct_m <- as.matrix(ct)

  bay <- fit_bayesian_lm_jeffreys(
    model_df = as.data.frame(sat_df_w),
    formula_obj = sat_formula,
    n_draws = n_draws_sat,
    seed = seed_sat
  )

  ols_dt <- data.table(
    window = start_year,
    term = rownames(ct_m),
    estimate = ct_m[, 1],
    se = ct_m[, 2],
    t_value = ct_m[, 3],
    p_value = ct_m[, ncol(ct_m)],
    se_type = if (have_hac) "Newey-West HAC" else "classical"
  )

  bay_dt <- as.data.table(
    summarise_draw_matrix(bay$beta_draws, ols = bay$beta_ols)
  )
  bay_dt[, window := start_year]

  fwrite(ols_dt, file.path(out_dir, sprintf("satellite_ols_hac_%d.csv", start_year)))
  fwrite(bay_dt, file.path(out_dir, sprintf("satellite_bayes_summary_%d.csv", start_year)))

  list(
    window = start_year, Td = nrow(d_w), n_sat = stats::nobs(fit),
    p_ttc = z_out$p_ttc, rho = z_out$rho,
    r2 = sm$r.squared, adj_r2 = sm$adj.r.squared, sigma = sm$sigma,
    ols = ols_dt, bayes = bay_dt
  )
}

cat("\n========== Short-sample satellite regressions ==========\n")
runs <- lapply(SS_START_YEARS, run_window)
names(runs) <- as.character(SS_START_YEARS)

ols_all <- rbindlist(lapply(runs, function(r) r$ols))
bay_all <- rbindlist(lapply(runs, function(r) r$bayes))
fwrite(ols_all, file.path(out_dir, "satellite_ols_hac_all_windows.csv"))
fwrite(bay_all, file.path(out_dir, "satellite_bayes_summary_all_windows.csv"))

meta <- rbindlist(lapply(runs, function(r) {
  data.table(window = r$window, Td = r$Td, n_sat = r$n_sat,
             p_ttc = r$p_ttc, rho = r$rho,
             r2 = r$r2, adj_r2 = r$adj_r2, sigma = r$sigma)
}))
fwrite(meta, file.path(out_dir, "satellite_fit_by_window.csv"))

# --- Combined LaTeX table: coefficients side by side ---------------------------
fmt <- function(x, d = 4) formatC(x, format = "f", digits = d)
stars <- function(p) {
  ifelse(p < 0.01, "$^{***}$", ifelse(p < 0.05, "$^{**}$", ifelse(p < 0.10, "$^{*}$", "")))
}

var_pretty <- c(
  "(Intercept)"        = "Intercept",
  "log_inv_pc_lag0"    = "log(Investment p.c.)$_{t}$",
  "log_inv_pc_lag2"    = "log(Investment p.c.)$_{t-2}$",
  "log_gdp_pc_lag2"    = "log(GDP p.c.)$_{t-2}$",
  "log_gdp_pc_lag4"    = "log(GDP p.c.)$_{t-4}$",
  "log_oil_real_lag0"  = "log(Oil real)$_{t}$",
  "log_oil_real_lag4"  = "log(Oil real)$_{t-4}$",
  "infl_yoy_pct_lag4"  = "Inflation YoY$_{t-4}$"
)
pretty_term <- function(x) {
  out <- unname(var_pretty[x])
  out[is.na(out)] <- gsub("_", "\\\\_", x[is.na(out)])
  out
}

win_head <- paste(sprintf("%d--2024", SS_START_YEARS), collapse = " & ")

tex_path <- file.path(out_dir, "satellite_coefficients_by_window.tex")
con <- file(tex_path, open = "wt")
cat("\\begin{table}[!htbp]\n\\centering\n", file = con)
cat("\\caption{Satellite coefficients across default-sample windows (fixed best-BIC design)}\n", file = con)
cat("\\label{tab:dralacbn_satellite_by_window}\n", file = con)
cat("\\begin{tabular}{l", paste(rep("c", length(SS_START_YEARS)), collapse = ""),
    "}\n\\toprule\n", sep = "", file = con)
cat(" & ", win_head, " \\\\\n\\midrule\n", sep = "", file = con)

for (tm in all_terms) {
  cells <- vapply(SS_START_YEARS, function(w) {
    r <- ols_all[window == w & term == tm]
    if (nrow(r) == 0L) return("--")
    paste0(fmt(r$estimate), stars(r$p_value))
  }, character(1))
  se_cells <- vapply(SS_START_YEARS, function(w) {
    r <- ols_all[window == w & term == tm]
    if (nrow(r) == 0L) return("")
    paste0("(", fmt(r$se), ")")
  }, character(1))
  cat(pretty_term(tm), " & ", paste(cells, collapse = " & "), " \\\\\n",
      sep = "", file = con)
  cat(" & ", paste(se_cells, collapse = " & "), " \\\\\n", sep = "", file = con)
}

cat("\\midrule\n", file = con)
cat("Observations & ",
    paste(meta$n_sat, collapse = " & "), " \\\\\n", sep = "", file = con)
cat("$R^2$ & ",
    paste(fmt(meta$r2, 3), collapse = " & "), " \\\\\n", sep = "", file = con)
cat("$\\hat p$ (TTC, \\%) & ",
    paste(fmt(100 * meta$p_ttc, 2), collapse = " & "), " \\\\\n", sep = "", file = con)
cat("$\\hat\\rho$ & ",
    paste(fmt(meta$rho, 4), collapse = " & "), " \\\\\n", sep = "", file = con)
cat("\\bottomrule\n\\end{tabular}\n", file = con)
cat("\\begin{tablenotes}\n\\small\n", file = con)
cat("\\item Notes: OLS estimates of the satellite equation on each truncated default sample, with Newey--West HAC standard errors in parentheses (", NB_LAGS_SAT, " lags). The design is fixed at the full-sample best-BIC specification; the dependent variable is the systematic factor reconstructed within each window with the window-specific Merton parameters reported in the bottom panel. Significance: $^{*}p<0.10$, $^{**}p<0.05$, $^{***}p<0.01$.\n", sep = "", file = con)
cat("\\end{tablenotes}\n\\end{table}\n", file = con)
close(con)

cat("\n--- Coefficient stability across windows (OLS, HAC) --------------------\n")
wide <- dcast(ols_all, term ~ window, value.var = "estimate")
print(wide)
cat("\nMerton parameters by window:\n")
print(meta[, .(window, Td, p_ttc = round(p_ttc, 4), rho = round(rho, 4),
               r2 = round(r2, 3))])
cat("\nOutputs written to: ", out_dir, "\n", sep = "")

invisible(TRUE)
