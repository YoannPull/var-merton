# =============================================================================
#  code/shared/_application_information_set_robustness.R
#
#  Shared information-set robustness engine for the EBA and DRALACBN
#  applications. Application-specific wrappers must define the APP_* objects
#  listed below, then source this file.
#
#  Required objects supplied by wrapper:
#    APP_NAME                 character scalar, e.g. "eba" or "dralacbn"
#    APP_LABEL                character scalar used in messages/tables
#    APP_MAIN_OUT_DIR         main-specification output root
#    APP_ROBUST_OUT_DIR       robustness output root for alternative specifications
#    APP_Z_DT                 data.table/data.frame with Date and Z columns
#    APP_P_TTC                through-the-cycle default probability
#    APP_RHO                  asset correlation used in Merton--Vasicek mapping
#    APP_INFO_SETS_TO_COMPARE character vector of VAR_SPECS names
#    APP_BASELINE_INFO_SET    one element of APP_INFO_SETS_TO_COMPARE
#    APP_Z_EXPORT             optional data.table to write at the output root
#    APP_Z_EXPORT_NAME        optional CSV file name for APP_Z_EXPORT
#    APP_NOTE                 optional character note for metadata
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(lubridate)
})

required_app_objects <- c(
  "APP_NAME", "APP_LABEL", "APP_Z_DT", "APP_P_TTC", "APP_RHO",
  "APP_INFO_SETS_TO_COMPARE", "APP_BASELINE_INFO_SET"
)

missing_app_objects <- required_app_objects[
  !vapply(required_app_objects, exists, logical(1))
]

if (length(missing_app_objects) > 0L) {
  stop(
    "Shared robustness engine missing required APP_* objects: ",
    paste(missing_app_objects, collapse = ", "),
    call. = FALSE
  )
}

set.seed(if (exists("SEED_ROBUSTNESS")) SEED_ROBUSTNESS else 123)

# -----------------------------------------------------------------------------
# 0. Paths and generic parameters
# -----------------------------------------------------------------------------

data_var_path <- if (exists("PATH_DATA_VAR")) {
  PATH_DATA_VAR
} else {
  "data/processed/data_var_for_model.csv"
}

kernel_dir <- if (exists("DIR_VAR_KERNELS")) {
  DIR_VAR_KERNELS
} else {
  file.path("output", "shared", "var_kernels")
}

main_out_dir <- if (exists("APP_MAIN_OUT_DIR")) {
  APP_MAIN_OUT_DIR
} else if (exists("APP_OUT_DIR")) {
  file.path(APP_OUT_DIR, "main_specification")
} else {
  file.path("output", "applications", APP_NAME, "main_specification")
}

robust_out_dir <- if (exists("APP_ROBUST_OUT_DIR")) {
  APP_ROBUST_OUT_DIR
} else if (exists("APP_OUT_DIR")) {
  APP_OUT_DIR
} else {
  file.path("output", "robustness", APP_NAME, "information_set")
}

dir.create(main_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(robust_out_dir, recursive = TRUE, showWarnings = FALSE)

p_lags <- if (exists("P_LAGS")) P_LAGS else 2
nrep <- if (exists("NREP")) NREP else 30000
M_target <- if (exists("M_TARGET")) M_TARGET else 10000L
H <- if (exists("HORIZON")) HORIZON else 12
seed <- if (exists("SEED_ROBUSTNESS")) SEED_ROBUSTNESS else 123

impulse_name <- if (exists("IMPULSE_NAME")) IMPULSE_NAME else "log_GPRD"

NB_LAGS_SAT <- if (exists("ALT_SPEC_NB_LAGS_SAT")) {
  ALT_SPEC_NB_LAGS_SAT
} else if (exists("NB_LAGS")) {
  NB_LAGS
} else {
  4
}

MAX_VARIABLES_IN_MODEL <- if (exists("ALT_SPEC_MAX_VARIABLES_IN_MODEL")) {
  ALT_SPEC_MAX_VARIABLES_IN_MODEL
} else {
  6
}

PRESELECT_BY_AIC_TOPN <- if (exists("PRESELECT_BY_AIC_TOPN")) {
  PRESELECT_BY_AIC_TOPN
} else {
  1000L
}

n_draws_sat <- if (exists("N_DRAWS_SATELLITE")) {
  N_DRAWS_SATELLITE
} else if (exists("N_DRAWS_Z_SATELLITE")) {
  N_DRAWS_Z_SATELLITE
} else {
  M_target
}

seed_sat <- if (exists("SEED_ZFACTOR")) SEED_ZFACTOR else 12345

# -----------------------------------------------------------------------------
# Satellite model averaging (main specification only) and OOS policy.
# -----------------------------------------------------------------------------
source("code/shared/_model_averaging.R")
source("code/shared/_satellite_regularization.R")
source("code/shared/_horseshoe_satellite.R")
source("code/shared/_pd_var.R")

APP_RUN_PD_VAR <- if (exists("RUN_PD_VAR")) isTRUE(RUN_PD_VAR) else TRUE
APP_PD_VAR_ALPHAS <- if (exists("PD_VAR_ALPHAS")) PD_VAR_ALPHAS else c(0.99, 0.999)

APP_RUN_SATELLITE_REGULARIZATION <- if (exists("RUN_SATELLITE_REGULARIZATION")) {
  isTRUE(RUN_SATELLITE_REGULARIZATION)
} else {
  TRUE
}

APP_RUN_HORSESHOE <- if (exists("RUN_HORSESHOE_ROBUSTNESS")) {
  isTRUE(RUN_HORSESHOE_ROBUSTNESS)
} else {
  TRUE
}

APP_USE_MODEL_AVERAGING <- if (exists("USE_MODEL_AVERAGING_MAIN_SPEC")) {
  isTRUE(USE_MODEL_AVERAGING_MAIN_SPEC)
} else {
  FALSE
}

APP_BMA_WEIGHT        <- if (exists("BMA_WEIGHT")) BMA_WEIGHT else "aic"
APP_BMA_RULE          <- if (exists("BMA_RULE")) BMA_RULE else "cum_weight"
APP_BMA_CUM_THRESHOLD <- if (exists("BMA_CUM_THRESHOLD")) BMA_CUM_THRESHOLD else 0.95
APP_BMA_TOP_N         <- if (exists("BMA_TOP_N")) BMA_TOP_N else 10L
APP_BMA_DELTA_MAX     <- if (exists("BMA_DELTA_MAX")) BMA_DELTA_MAX else 2

# Out-of-sample evaluation is dropped for short-sample applications (EBA).
APP_DROP_OOS <- if (identical(APP_NAME, "eba") && exists("EBA_DROP_OOS")) {
  isTRUE(EBA_DROP_OOS)
} else {
  FALSE
}

model_avg_dir <- if (exists("DIR_MODEL_AVERAGING")) {
  file.path(DIR_MODEL_AVERAGING, APP_NAME)
} else {
  file.path("output", "model_averaging", APP_NAME)
}
dir.create(model_avg_dir, recursive = TRUE, showWarnings = FALSE)

options(robustness.core_fraction = getOption("robustness.core_fraction", 0.6))

if (!exists("VAR_SPECS")) {
  stop(APP_LABEL, " robustness: VAR_SPECS not found. Source config.R first.",
       call. = FALSE)
}

missing_info_sets <- setdiff(APP_INFO_SETS_TO_COMPARE, names(VAR_SPECS))

if (length(missing_info_sets) > 0L) {
  stop(
    APP_LABEL, " robustness: missing information sets in VAR_SPECS: ",
    paste(missing_info_sets, collapse = ", "),
    call. = FALSE
  )
}

selected_info_sets <- APP_INFO_SETS_TO_COMPARE[
  APP_INFO_SETS_TO_COMPARE %in% names(VAR_SPECS)
]

if (length(selected_info_sets) < 1L) {
  stop(APP_LABEL, " robustness: no valid information set selected.",
       call. = FALSE)
}

if (!(APP_BASELINE_INFO_SET %in% selected_info_sets)) {
  stop(
    APP_LABEL, " robustness: APP_BASELINE_INFO_SET = '",
    APP_BASELINE_INFO_SET,
    "' is not in APP_INFO_SETS_TO_COMPARE.",
    call. = FALSE
  )
}

INFO_SETS_TO_COMPARE <- VAR_SPECS[selected_info_sets]

INFO_SET_LABELS_LOCAL <- if (exists("INFO_SET_LABELS")) {
  INFO_SET_LABELS
} else {
  c(
    real_side      = "Real-side",
    rates_activity = "Rates/activity",
    caldara_style  = "Caldara-style",
    baseline       = "Baseline"
  )
}

info_set_labels <- selected_info_sets
names(info_set_labels) <- selected_info_sets

matched_labels <- INFO_SET_LABELS_LOCAL[selected_info_sets]
info_set_labels[!is.na(matched_labels)] <- matched_labels[!is.na(matched_labels)]

# -----------------------------------------------------------------------------
# Publication-friendly variable labels
# -----------------------------------------------------------------------------
# Notes:
#   - These labels are used in figures, especially facet titles.
#   - They are intentionally shorter than the full variable definitions.
#   - Full definitions should remain in variable_definitions.R / tables.
#
# Key variable:
#   log_private_pc = 100 * log(USPRIV / CNP16OV)
#   where USPRIV is private employment and CNP16OV is the civilian
#   noninstitutional population aged 16 and over.
# -----------------------------------------------------------------------------

var_labels <- if (exists("VAR_LABELS")) {
  VAR_LABELS
} else {
  c(
    "log_GPRD"               = "Geopolitical Risk (log)",
    "vix"                    = "VIX",
    "log_sp500_real"         = "Real S&P 500 (log)",
    "log_oil_real"           = "Real Oil Price (log)",
    "log_private_pc"         = "Private Employment (p.c., log)",
    "log_gdp_pc"             = "Real GDP (p.c., log)",
    "log_inv_pc"             = "Investment (p.c., log)",
    "infl_yoy_pct"           = "Inflation YoY (%)",
    "gs2"                    = "2-Year Treasury Yield",
    "nfci"                   = "NFCI",
    "epu"                    = "EPU",
    "t10Y3M"                 = "10Y-3M Term Spread",
    "t10Y2Y"                 = "10Y-2Y Term Spread",
    "gdp_yoy_pct"            = "GDP Growth YoY (%)",
    "gdp_annualized_pct"     = "GDP Growth Ann. QoQ (%)",
    "infl_annualized_pct"    = "Inflation Ann. QoQ (%)",
    "log_payems"             = "Nonfarm Payrolls (log)",
    "unrate"                 = "Unemployment Rate (%)"
  )
}

# Enforce short labels even when VAR_LABELS is supplied externally.
# This avoids truncated facet labels in robustness figures.
label_overrides <- c(
  "log_private_pc"      = "Private Employment (p.c., log)",
  "log_gdp_pc"          = "Real GDP (p.c., log)",
  "log_inv_pc"          = "Investment (p.c., log)",
  "gdp_annualized_pct"  = "GDP Growth Ann. QoQ (%)",
  "infl_annualized_pct" = "Inflation Ann. QoQ (%)"
)

common_override_names <- intersect(names(label_overrides), names(var_labels))
var_labels[common_override_names] <- label_overrides[common_override_names]

# Helper for wrapping long facet labels without requiring extra packages.
wrap_facet_label <- function(x, width = 22) {
  vapply(
    as.character(x),
    function(s) paste(strwrap(s, width = width), collapse = "\n"),
    character(1)
  )
}

APP_Z_DT <- as.data.table(copy(APP_Z_DT))

if (!all(c("Date", "Z") %in% names(APP_Z_DT))) {
  stop("APP_Z_DT must contain Date and Z columns.", call. = FALSE)
}

APP_Z_DT[, Date := as.Date(Date)]
APP_Z_DT <- APP_Z_DT[is.finite(Z)]
setorder(APP_Z_DT, Date)

APP_P_TTC <- as.numeric(APP_P_TTC)[1]
APP_RHO <- as.numeric(APP_RHO)[1]

stopifnot(is.finite(APP_P_TTC), APP_P_TTC > 0, APP_P_TTC < 1)
stopifnot(is.finite(APP_RHO), APP_RHO > 0, APP_RHO < 1)

if (exists("APP_Z_EXPORT")) {
  export_name <- if (exists("APP_Z_EXPORT_NAME")) {
    APP_Z_EXPORT_NAME
  } else {
    paste0("zfactor_", APP_NAME, ".csv")
  }
  fwrite(as.data.table(APP_Z_EXPORT), file.path(main_out_dir, export_name))
}

DT_raw <- safe_fread(data_var_path)

if (!"Date" %in% names(DT_raw)) {
  stop("data_var_for_model.csv must contain Date.", call. = FALSE)
}

DT_raw[, Date := as.Date(Date)]

cat("\n========== ", APP_LABEL, " information-set robustness ==========\n",
    sep = "")
cat("Compared information sets: ", paste(selected_info_sets, collapse = ", "),
    "\n", sep = "")
cat("Baseline information set: ", APP_BASELINE_INFO_SET, "\n", sep = "")
cat("Merton parameters: p_ttc = ", sprintf("%.8f", APP_P_TTC),
    ", rho = ", sprintf("%.8f", APP_RHO), "\n", sep = "")
cat("Satellite selection: OOS/AIC screen, max variables = ",
    MAX_VARIABLES_IN_MODEL, ", lags = ", NB_LAGS_SAT, "\n", sep = "")

# -----------------------------------------------------------------------------
# 1. Local helpers
# -----------------------------------------------------------------------------

peak_info <- function(bands) {
  i <- which.max(abs(bands$median))
  list(peak = bands$median[i], horizon = bands$horizon[i])
}

save_png_plot <- function(plot, dir, stem, width, height, dpi = 300) {
  # ggsave_both writes a PNG and a vector PDF twin (preferred for the paper).
  ggsave_both(plot, dir, stem, width = width, height = height, dpi = dpi)
}

plot_overlay <- function(dt_all, ylab_txt, baseline = APP_BASELINE_INFO_SET) {
  dt_all <- as.data.table(copy(dt_all))
  
  dt_all[
    ,
    information_set := factor(
      as.character(information_set),
      levels = selected_info_sets
    )
  ]
  
  dt_all[
    ,
    info_set_label := factor(
      as.character(info_set_label),
      levels = unname(info_set_labels[selected_info_sets])
    )
  ]
  
  base_band <- dt_all[as.character(information_set) == baseline]
  
  baseline_label <- if (baseline %in% names(info_set_labels)) {
    info_set_labels[[baseline]]
  } else {
    baseline
  }
  
  ggplot(dt_all, aes(x = horizon, y = median, color = information_set)) +
    geom_ribbon(
      data = base_band,
      aes(x = horizon, ymin = lower68, ymax = upper68),
      inherit.aes = FALSE,
      fill = "#AB4A7D",
      alpha = 0.16
    ) +
    geom_hline(
      yintercept = 0,
      color = "grey45",
      linewidth = 0.35,
      linetype = "dashed"
    ) +
    geom_line(linewidth = 1.0) +
    scale_color_discrete(
      breaks = selected_info_sets,
      labels = unname(info_set_labels[selected_info_sets]),
      name = "Information set"
    ) +
    scale_x_continuous(
      breaks = function(x) seq(max(0, ceiling(x[1])), floor(x[2]), by = 1),
      minor_breaks = NULL
    ) +
    labs(
      x = "Horizon (quarters)",
      y = ylab_txt,
      caption = paste0(
        "Shaded band: 68% posterior interval of the ",
        baseline_label,
        " baseline."
      )
    ) +
    theme_robustness() +
    theme(plot.caption = element_blank())
}

plot_pd_bands_all_information_sets <- function(dt_all) {
  dt_all <- as.data.table(copy(dt_all))
  
  dt_all[
    ,
    information_set := factor(
      as.character(information_set),
      levels = selected_info_sets
    )
  ]
  
  dt_all[
    ,
    info_set_label := factor(
      as.character(info_set_label),
      levels = unname(info_set_labels[selected_info_sets])
    )
  ]
  
  dt_all[
    ,
    shock_label := fifelse(
      shock == "one_sd",
      "One-s.d. GPR shock",
      "2001Q3-sized GPR shock"
    )
  ]
  
  dt_all[
    ,
    shock_label := factor(
      shock_label,
      levels = c("One-s.d. GPR shock", "2001Q3-sized GPR shock")
    )
  ]
  
  ggplot(dt_all, aes(x = horizon, y = median)) +
    geom_ribbon(aes(ymin = lower90, ymax = upper90),
                fill = "grey70", alpha = 0.35) +
    geom_ribbon(aes(ymin = lower68, ymax = upper68),
                fill = "grey45", alpha = 0.35) +
    geom_hline(
      yintercept = 0,
      color = "grey35",
      linewidth = 0.35,
      linetype = "dashed"
    ) +
    geom_line(color = "#6F1732", linewidth = 0.9) +
    facet_grid(shock_label ~ info_set_label, scales = "free_y") +
    scale_x_continuous(
      breaks = function(x) seq(max(0, ceiling(x[1])), floor(x[2]), by = 2),
      minor_breaks = NULL
    ) +
    labs(
      x = "Horizon (quarters)",
      y = expression(Delta ~ PD),
      caption = "Dark band: 68% posterior interval. Light band: 90% posterior interval."
    ) +
    theme_robustness() +
    theme(
      strip.text = element_text(face = "bold"),
      plot.caption = element_blank(),
      legend.position = "none"
    )
}

plot_macro_bands_all_information_sets <- function(dt_all) {
  dt_all <- as.data.table(copy(dt_all))
  
  dt_all[
    ,
    information_set := factor(
      as.character(information_set),
      levels = selected_info_sets
    )
  ]
  
  dt_all[
    ,
    info_set_label := factor(
      as.character(info_set_label),
      levels = unname(info_set_labels[selected_info_sets])
    )
  ]
  
  dt_all[
    ,
    shock_label := fifelse(
      shock == "one_sd",
      "One-s.d. GPR shock",
      "2001Q3-sized GPR shock"
    )
  ]
  
  dt_all[
    ,
    shock_label := factor(
      shock_label,
      levels = c("One-s.d. GPR shock", "2001Q3-sized GPR shock")
    )
  ]
  
  dt_all[, variable_label := as.character(variable)]
  dt_all[, variable_label := unname(var_labels[variable_label])]
  dt_all[is.na(variable_label), variable_label := as.character(variable)]

  # Order facet columns by variable ordering (impulse variable first),
  # instead of ggplot's default alphabetical ordering.
  dt_all[, variable_label := factor(variable_label, levels = unique(variable_label))]

  ggplot(dt_all, aes(x = horizon, y = median)) +
    geom_hline(
      yintercept = 0,
      color = "grey45",
      linewidth = 0.30,
      linetype = "dashed"
    ) +
    geom_line(color = "#6F1732", linewidth = 0.75) +
    facet_grid(
      shock_label + info_set_label ~ variable_label,
      scales = "free_y",
      labeller = labeller(
        variable_label = function(x) wrap_facet_label(x, width = 20),
        info_set_label = function(x) wrap_facet_label(x, width = 18),
        shock_label = function(x) wrap_facet_label(x, width = 24)
      )
    ) +
    scale_x_continuous(
      breaks = function(x) seq(max(0, ceiling(x[1])), floor(x[2]), by = 3),
      minor_breaks = NULL
    ) +
    labs(x = "Horizon (quarters)", y = "VAR GIRF") +
    theme_robustness() +
    theme(
      strip.text = element_text(face = "bold", size = 8),
      legend.position = "none"
    )
}

write_information_set_comparison_tex <- function(dt, path) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  
  caption_app <- gsub("_", "\\\\_", APP_LABEL)
  
  cat("\\begin{table}[!htbp]\n\\centering\n", file = con)
  cat("\\caption{", caption_app,
      " credit-risk responses across VAR information sets}\n",
      sep = "", file = con)
  cat("\\label{tab:", APP_NAME, "_information_set_robustness}\n",
      sep = "", file = con)
  cat("\\begin{tabular}{lrrrrr}\n\\toprule\n", file = con)
  cat(
    "Information set & $k$ & Sat. $R^2$ & OOS RMSE & Peak $\\Delta$PD (1 s.d.) & Peak $\\Delta$PD (2001Q3) \\\\\n",
    file = con
  )
  cat("\\midrule\n", file = con)
  
  for (i in seq_len(nrow(dt))) {
    cat(
      dt$info_set_label[i], " & ",
      dt$k[i], " & ",
      fmt_num(dt$r2[i], 3), " & ",
      fmt_num(dt$oos_rmse[i], 4), " & ",
      sprintf("%.3e (h=%d)",
              dt$peak_dpd_1sd[i],
              as.integer(dt$peak_h_1sd[i])), " & ",
      sprintf("%.3e (h=%d)",
              dt$peak_dpd_2001[i],
              as.integer(dt$peak_h_2001[i])), " \\\\\n",
      sep = "",
      file = con
    )
  }
  
  cat("\\bottomrule\n\\end{tabular}\n", file = con)
  cat("\\begin{tablenotes}\n\\small\n", file = con)
  cat(
    "\\item Notes: The Merton--Vasicek $Z$, asset correlation $\\rho$ and TTC PD are fixed within the application. Across information sets, only VAR transmission and the selected satellite bridge change. Peak $\\Delta$PD is the largest-magnitude posterior-median PD response and its horizon.\n",
    file = con
  )
  
  if (exists("APP_NOTE")) {
    cat("\\item ", APP_NOTE, "\n", sep = "", file = con)
  }
  
  cat("\\end{tablenotes}\n\\end{table}\n", file = con)
}

# -----------------------------------------------------------------------------
# 2. Transmission for one VAR information set
# -----------------------------------------------------------------------------

run_transmission_for_information_set <- function(info_set_name, vars_extra) {
  info_set_label <- if (info_set_name %in% names(info_set_labels)) {
    info_set_labels[[info_set_name]]
  } else {
    info_set_name
  }
  
  cat("\n----- ", APP_LABEL, " | information set: ", info_set_name,
      " (", info_set_label, ") -----\n", sep = "")
  
  info_set_dir <- if (identical(info_set_name, APP_BASELINE_INFO_SET)) {
    main_out_dir
  } else {
    file.path(robust_out_dir, "by_information_set", info_set_name)
  }
  
  for (sub in c(
    "satellite",
    "figures",
    file.path("irf", "VAR"),
    file.path("irf", "Z"),
    file.path("irf", "PD"),
    file.path("irf", "diagnostics")
  )) {
    dir.create(file.path(info_set_dir, sub),
               recursive = TRUE,
               showWarnings = FALSE)
  }
  
  kern <- get_or_estimate_kernel(
    spec_name = info_set_name,
    vars_extra = vars_extra,
    DT_raw = DT_raw,
    p_lags = p_lags,
    nrep = nrep,
    H = H,
    seed = seed,
    M_target = M_target,
    kernel_dir = kernel_dir,
    impulse_name = impulse_name
  )
  
  variables <- kern$variables
  impulse_ix <- match(impulse_name, variables)
  
  if (is.na(impulse_ix)) {
    stop("Impulse variable not found in ", info_set_name, " VAR.",
         call. = FALSE)
  }
  
  vars_sat <- setdiff(variables, impulse_name)
  
  sat_df <- prepare_lagged_satellite_data(
    z_input = APP_Z_DT,
    DT_raw = DT_raw,
    vars_sat = vars_sat,
    nb_lags = NB_LAGS_SAT
  )
  
  # Model averaging (BMA) is applied to every VAR information set (the main
  # specification and the robustness specifications), so each one is summarised
  # by an Akaike-weight average rather than a single selected satellite.
  use_bma_here <- APP_USE_MODEL_AVERAGING

  if (use_bma_here) {
    cat("    Satellite: full-sample model averaging (",
        APP_BMA_WEIGHT, " weights, rule = ", APP_BMA_RULE, ").\n", sep = "")

    sat_screen <- select_satellite_fullsample(
      model_df = sat_df,
      max_vars = MAX_VARIABLES_IN_MODEL,
      p_threshold = NULL,
      preselect_by_aic_topn = PRESELECT_BY_AIC_TOPN
    )

    if (is.null(sat_screen$best) || is.null(sat_screen$best$model)) {
      stop("No valid satellite model selected for ", APP_LABEL, " / ",
           info_set_name, call. = FALSE)
    }

    # The single best model is retained for the OLS coefficient table / labels.
    sat_model <- sat_screen$best$model
    sat_summary <- summary(sat_model)

    # The mixture Bayesian satellite is what feeds the GIRF -> Z -> PD pipeline.
    bsat <- make_bma_satellite(
      candidates = sat_screen$candidates,
      X_names = sat_screen$X_names,
      model_df = sat_df,
      variables = variables,
      weight = APP_BMA_WEIGHT,
      rule = APP_BMA_RULE,
      cum_threshold = APP_BMA_CUM_THRESHOLD,
      top_n = APP_BMA_TOP_N,
      delta_max = APP_BMA_DELTA_MAX,
      n_draws = n_draws_sat,
      seed = seed_sat
    )

    write_bma_diagnostics(
      bsat = bsat,
      out_dir = model_avg_dir,
      app_name = APP_NAME,
      info_set_name = info_set_name,
      best_model = sat_model
    )

    write_bma_regression_table(
      bsat = bsat,
      out_dir = model_avg_dir,
      app_name = APP_NAME,
      info_set_name = info_set_name,
      info_set_label = info_set_label,
      var_labels = var_labels,
      n_obs = stats::nobs(sat_model),
      app_label = APP_LABEL
    )

    # Out-of-sample diagnostics: dropped for short-sample applications (EBA),
    # otherwise scored on the single best model for the comparison table.
    oos <- if (APP_DROP_OOS) {
      list(rmse = NA_real_, mae = NA_real_, r2 = NA_real_, preds = NULL)
    } else {
      rolling_oos_error(
        as.data.table(sat_df),
        formula(sat_model),
        initial = max(24L, ceiling(nrow(sat_df) / 2))
      )
    }
  } else {
    sat_sel <- select_satellite_oos(
      model_df = sat_df,
      max_vars = MAX_VARIABLES_IN_MODEL,
      p_threshold = NULL,
      preselect_by_aic_topn = PRESELECT_BY_AIC_TOPN
    )

    if (is.null(sat_sel$best) || is.null(sat_sel$best$model)) {
      stop("No valid satellite model selected for ", APP_LABEL, " / ",
           info_set_name, call. = FALSE)
    }

    sat_model <- sat_sel$best$model
    oos <- sat_sel$best$oos
    sat_summary <- summary(sat_model)

    bsat <- make_bayesian_satellite(
      sat_model = sat_model,
      model_df = sat_df,
      variables = variables,
      n_draws = n_draws_sat,
      seed = seed_sat
    )
  }

  saveRDS(sat_model, file.path(info_set_dir, "satellite", "best_satellite.rds"))
  
  sat_coefs <- as.data.table(
    as.data.frame(sat_summary$coefficients),
    keep.rownames = "term"
  )
  
  setnames(
    sat_coefs,
    old = c("Estimate", "Std. Error", "t value", "Pr(>|t|)"),
    new = c("estimate", "std_error", "t_value", "p_value"),
    skip_absent = TRUE
  )
  
  fwrite(
    sat_coefs,
    file.path(info_set_dir, "satellite", "satellite_ols_coefficients.csv")
  )
  
  write_satellite_table_tex(
    sat_model,
    file.path(info_set_dir, "satellite", "satellite_ols_table.tex"),
    caption = paste0(
      APP_LABEL,
      ": selected satellite under the ",
      info_set_label,
      " information set"
    ),
    label = paste0("tab:", APP_NAME, "_", info_set_name, "_satellite"),
    extra_rows = list(
      "Information set" = info_set_label,
      "Selection" = if (use_bma_here) {
        paste0("Akaike model averaging (", bsat$bma$n_models, " models)")
      } else {
        "OOS/AIC screen"
      }
    ),
    notes = "The dependent variable is the reconstructed systematic factor Z. Bayesian posterior draws are generated conditionally on this selected linear satellite."
  )

  # Robust SEs (HAC Newey-West + HC3) for the selected satellite. Frequentist
  # diagnostic alongside the model-averaged Bayesian satellite; every spec.
  write_satellite_robust_se(
    sat_model = sat_model,
    out_dir = file.path(info_set_dir, "satellite"),
    nb_lags = NB_LAGS_SAT,
    app_label = APP_LABEL, info_set_label = info_set_label,
    app_name = APP_NAME, info_set_name = info_set_name
  )

  # Ridge/Lasso shrinkage robustness of the satellite (Z) coefficients, mirroring
  # the EBA z-factor stage but applied to every information set of both apps.
  if (APP_RUN_SATELLITE_REGULARIZATION) {
    tryCatch(
      run_satellite_regularization(
        model_df = sat_df,
        out_dir = file.path(info_set_dir, "satellite", "regularization"),
        label = paste0(APP_LABEL, " / ", info_set_name),
        run_ridge = if (exists("RUN_RIDGE")) RUN_RIDGE else TRUE,
        run_lasso = if (exists("RUN_LASSO")) RUN_LASSO else TRUE,
        run_elastic_net = if (exists("RUN_ELASTIC_NET")) RUN_ELASTIC_NET else FALSE,
        use_lambda = if (exists("USE_LAMBDA")) USE_LAMBDA else "lambda.1se",
        seed = seed_sat
      ),
      error = function(e) {
        cat("    Ridge/Lasso robustness skipped for ", info_set_name, ": ",
            conditionMessage(e), "\n", sep = "")
      }
    )
  }

  # bsat (the Bayesian satellite feeding the GIRF -> Z -> PD pipeline) is built
  # above: a BMA mixture for the main specification, otherwise the single
  # OOS/AIC-selected satellite.
  saveRDS(
    bsat,
    file.path(info_set_dir, "satellite", "bayesian_satellite.rds")
  )
  
  fwrite(
    summarise_draw_matrix(
      cbind(`(Intercept)` = bsat$beta0_draws, bsat$beta_term_draws),
      ols = bsat$beta_ols
    ),
    file.path(
      info_set_dir,
      "satellite",
      "bayesian_satellite_coefficients_summary.csv"
    )
  )
  
  fwrite(
    summarise_vector_draws(bsat$sigma2_draws, name = "sigma2"),
    file.path(
      info_set_dir,
      "satellite",
      "bayesian_satellite_sigma_summary.csv"
    )
  )
  
  moments <- compute_factor_moment_draws_bayes(
    kernel = kern,
    bsat = bsat,
    impulse_idx = impulse_ix,
    p_lags = p_lags,
    H = H
  )
  
  macro_1sd <- compute_var_girf(
    kern$Psi_draws,
    kern$Sigma_kept,
    impulse_idx = impulse_ix,
    shock_scale = 1
  )
  
  psiZ_1sd <- inject_girf_into_Z_bayes(
    kern$Psi_draws,
    kern$Sigma_kept,
    bsat,
    impulse_ix,
    shock_scale = 1
  )
  
  pd_1sd <- compute_pd_girf(
    psiZ_1sd$draws,
    moments$mu_draws,
    moments$s2_draws,
    moments$s2_delta_draws,
    p = APP_P_TTC,
    rho = APP_RHO
  )
  
  decomp_1sd <- compute_sobol_decomposition(
    Psi_draws = kern$Psi_draws,
    Sigma_kept = kern$Sigma_kept,
    B_kept = kern$B_kept,
    variables = variables,
    DT_var = as.data.table(kern$DT),
    terms_df = bsat$terms,
    beta0_draws = bsat$beta0_draws,
    beta_term_draws = bsat$beta_term_draws,
    sigma2_draws = bsat$sigma2_draws,
    impulse_idx = impulse_ix,
    p_lags = p_lags,
    H = H,
    p_val = APP_P_TTC,
    rho_val = APP_RHO,
    shock_scale = 1,
    seed = seed
  )
  
  shock_2001 <- compute_structural_e1_median(
    kernel = kern,
    dates_raw = as.Date(DT_raw$Date),
    target_quarter = "2001Q3",
    impulse_idx = impulse_ix
  )
  
  macro_2001 <- compute_var_girf(
    kern$Psi_draws,
    kern$Sigma_kept,
    impulse_idx = impulse_ix,
    shock_scale = shock_2001$shock_scale
  )
  
  psiZ_2001 <- inject_girf_into_Z_bayes(
    kern$Psi_draws,
    kern$Sigma_kept,
    bsat,
    impulse_ix,
    shock_scale = shock_2001$shock_scale
  )
  
  pd_2001 <- compute_pd_girf(
    psiZ_2001$draws,
    moments$mu_draws,
    moments$s2_draws,
    moments$s2_delta_draws,
    p = APP_P_TTC,
    rho = APP_RHO
  )
  
  decomp_2001 <- compute_sobol_decomposition(
    Psi_draws = kern$Psi_draws,
    Sigma_kept = kern$Sigma_kept,
    B_kept = kern$B_kept,
    variables = variables,
    DT_var = as.data.table(kern$DT),
    terms_df = bsat$terms,
    beta0_draws = bsat$beta0_draws,
    beta_term_draws = bsat$beta_term_draws,
    sigma2_draws = bsat$sigma2_draws,
    impulse_idx = impulse_ix,
    p_lags = p_lags,
    H = H,
    p_val = APP_P_TTC,
    rho_val = APP_RHO,
    shock_scale = shock_2001$shock_scale,
    seed = seed
  )
  
  tag_bands <- function(dt, shock_lab) {
    dt <- as.data.table(copy(dt))
    dt[
      ,
      `:=`(
        information_set = info_set_name,
        info_set_label = info_set_label,
        shock = shock_lab
      )
    ]
    dt[]
  }
  
  macro_1sd$bands <- tag_bands(macro_1sd$bands, "one_sd")
  macro_2001$bands <- tag_bands(macro_2001$bands, "2001Q3")
  psiZ_1sd$bands <- tag_bands(psiZ_1sd$bands, "one_sd")
  psiZ_2001$bands <- tag_bands(psiZ_2001$bands, "2001Q3")
  pd_1sd$bands <- tag_bands(pd_1sd$bands, "one_sd")
  pd_2001$bands <- tag_bands(pd_2001$bands, "2001Q3")
  
  fwrite(
    macro_1sd$bands,
    file.path(info_set_dir, "irf", "VAR", "macro_girf_bands_1sd.csv")
  )
  
  fwrite(
    macro_2001$bands,
    file.path(info_set_dir, "irf", "VAR", "macro_girf_bands_2001Q3.csv")
  )
  
  fwrite(
    psiZ_1sd$bands,
    file.path(info_set_dir, "irf", "Z", "psiZ_bands_1sd.csv")
  )
  
  fwrite(
    psiZ_2001$bands,
    file.path(info_set_dir, "irf", "Z", "psiZ_bands_2001Q3.csv")
  )
  
  fwrite(
    pd_1sd$bands,
    file.path(info_set_dir, "irf", "PD", "pd_girf_bands_1sd.csv")
  )
  
  fwrite(
    pd_2001$bands,
    file.path(info_set_dir, "irf", "PD", "pd_girf_bands_2001Q3.csv")
  )
  
  fwrite(
    shock_2001$e1_dt,
    file.path(info_set_dir, "irf", "VAR", "standardized_gpr_innovation.csv")
  )
  
  fwrite(
    shock_2001$target,
    file.path(info_set_dir, "irf", "VAR", "shock_scale_2001Q3.csv")
  )
  
  saveRDS(
    macro_1sd,
    file.path(info_set_dir, "irf", "VAR", "macro_girf_1sd.rds")
  )
  
  saveRDS(
    macro_2001,
    file.path(info_set_dir, "irf", "VAR", "macro_girf_2001Q3.rds")
  )
  
  saveRDS(
    psiZ_1sd,
    file.path(info_set_dir, "irf", "Z", "psiZ_1sd.rds")
  )
  
  saveRDS(
    psiZ_2001,
    file.path(info_set_dir, "irf", "Z", "psiZ_2001Q3.rds")
  )
  
  saveRDS(
    pd_1sd,
    file.path(info_set_dir, "irf", "PD", "pd_girf_1sd.rds")
  )
  
  saveRDS(
    pd_2001,
    file.path(info_set_dir, "irf", "PD", "pd_girf_2001Q3.rds")
  )

  # Mirror the model-averaged PD (and psi_Z) GIRFs into the model-averaging
  # folder, so the averaged-satellite responses are self-contained there.
  if (use_bma_here) {
    ma_irf_dir <- file.path(model_avg_dir, "irf")
    dir.create(ma_irf_dir, recursive = TRUE, showWarnings = FALSE)
    ma_stem <- paste0(APP_NAME, "_", info_set_name, "_")

    fwrite(pd_1sd$bands,
           file.path(ma_irf_dir, paste0(ma_stem, "pd_girf_bands_1sd.csv")))
    fwrite(pd_2001$bands,
           file.path(ma_irf_dir, paste0(ma_stem, "pd_girf_bands_2001Q3.csv")))
    fwrite(psiZ_1sd$bands,
           file.path(ma_irf_dir, paste0(ma_stem, "psiZ_girf_bands_1sd.csv")))
    fwrite(psiZ_2001$bands,
           file.path(ma_irf_dir, paste0(ma_stem, "psiZ_girf_bands_2001Q3.csv")))

    saveRDS(pd_1sd,
            file.path(ma_irf_dir, paste0(ma_stem, "pd_girf_1sd.rds")))
    saveRDS(pd_2001,
            file.path(ma_irf_dir, paste0(ma_stem, "pd_girf_2001Q3.rds")))

    save_png_plot(
      plot_single_bands(pd_1sd$bands, ylab_txt = expression(Delta ~ PD)),
      ma_irf_dir, paste0(ma_stem, "pd_girf_1sd"),
      width = 6.5, height = 4.2
    )
    save_png_plot(
      plot_single_bands(pd_2001$bands, ylab_txt = expression(Delta ~ PD)),
      ma_irf_dir, paste0(ma_stem, "pd_girf_2001Q3"),
      width = 6.5, height = 4.2
    )
  }

  # Value-at-Risk (stressed / downturn PD) response (see _pd_var.R).
  if (APP_RUN_PD_VAR) {
    run_pd_var_robustness(
      out_dir = file.path(info_set_dir, "irf", "PD_VaR"),
      psiZ_1sd = psiZ_1sd, psiZ_2001 = psiZ_2001,
      mean_bands_1sd = pd_1sd$bands, mean_bands_2001 = pd_2001$bands,
      moments = moments,
      app_name = APP_NAME, info_set_name = info_set_name,
      info_set_label = info_set_label,
      p = APP_P_TTC, rho = APP_RHO, alphas = APP_PD_VAR_ALPHAS
    )
  }

  # Horseshoe satellite robustness (see _horseshoe_satellite.R).
  if (APP_RUN_HORSESHOE) {
    run_horseshoe_robustness(
      out_dir = file.path(info_set_dir, "satellite", "horseshoe"),
      sat_df = sat_df, variables = variables,
      n_draws = n_draws_sat,
      burnin = if (exists("HORSESHOE_BURNIN")) HORSESHOE_BURNIN else 1000L,
      seed = seed_sat,
      kern = kern, impulse_ix = impulse_ix, p_lags = p_lags, H = H,
      shock_scale_2001 = shock_2001$shock_scale,
      p = APP_P_TTC, rho = APP_RHO,
      app_name = APP_NAME, app_label = APP_LABEL,
      info_set_name = info_set_name, info_set_label = info_set_label,
      var_labels = var_labels, n_obs = stats::nobs(sat_model)
    )
  }

  pd_dec_1sd <- as.data.table(copy(decomp_1sd))[
    ,
    `:=`(
      information_set = info_set_name,
      info_set_label = info_set_label,
      shock = "one_sd"
    )
  ]
  
  pd_dec_2001 <- as.data.table(copy(decomp_2001))[
    ,
    `:=`(
      information_set = info_set_name,
      info_set_label = info_set_label,
      shock = "2001Q3"
    )
  ]
  
  fwrite(
    pd_dec_1sd,
    file.path(
      info_set_dir,
      "irf",
      "diagnostics",
      "pd_variance_decomposition_90pct_1sd.csv"
    )
  )
  
  fwrite(
    pd_dec_2001,
    file.path(
      info_set_dir,
      "irf",
      "diagnostics",
      "pd_variance_decomposition_90pct_2001Q3.csv"
    )
  )
  
  # Per-specification IRF figures, for both main and robustness specifications.
  save_png_plot(
    plot_bands_facets(
      macro_1sd$bands,
      ylab_txt = "VAR GIRF",
      var_labels = var_labels
    ),
    file.path(info_set_dir, "figures"),
    "var_irf_1sd",
    width = 8.5,
    height = 5.6
  )
  
  save_png_plot(
    plot_bands_facets(
      macro_2001$bands,
      ylab_txt = "VAR GIRF",
      var_labels = var_labels
    ),
    file.path(info_set_dir, "figures"),
    "var_irf_2001Q3",
    width = 8.5,
    height = 5.6
  )
  
  save_png_plot(
    plot_single_bands(psiZ_1sd$bands, ylab_txt = expression(psi[Z](h))),
    file.path(info_set_dir, "figures"),
    "z_girf_1sd",
    width = 6.5,
    height = 4.2
  )
  
  save_png_plot(
    plot_single_bands(psiZ_2001$bands, ylab_txt = expression(psi[Z](h))),
    file.path(info_set_dir, "figures"),
    "z_girf_2001Q3",
    width = 6.5,
    height = 4.2
  )
  
  save_png_plot(
    plot_single_bands(pd_1sd$bands, ylab_txt = expression(Delta ~ PD)),
    file.path(info_set_dir, "figures"),
    "pd_girf_1sd",
    width = 6.5,
    height = 4.2
  )
  
  save_png_plot(
    plot_single_bands(pd_2001$bands, ylab_txt = expression(Delta ~ PD)),
    file.path(info_set_dir, "figures"),
    "pd_girf_2001Q3",
    width = 6.5,
    height = 4.2
  )
  
  save_png_plot(
    plot_variance_decomposition_shares_robust(decomp_1sd),
    file.path(info_set_dir, "figures"),
    "pd_variance_decomposition_90pct_1sd",
    width = 6.5,
    height = 4.2
  )
  
  save_png_plot(
    plot_variance_decomposition_shares_robust(decomp_2001),
    file.path(info_set_dir, "figures"),
    "pd_variance_decomposition_90pct_2001Q3",
    width = 6.5,
    height = 4.2
  )
  
  pk1 <- peak_info(pd_1sd$bands)
  pk2 <- peak_info(pd_2001$bands)
  
  summary_row <- data.table(
    application = APP_NAME,
    information_set = info_set_name,
    info_set_label = info_set_label,
    is_main_specification = identical(info_set_name, APP_BASELINE_INFO_SET),
    k = length(variables),
    n_obs_satellite = nobs(sat_model),
    satellite_terms = paste(
      setdiff(names(coef(sat_model)), "(Intercept)"),
      collapse = " + "
    ),
    r2 = sat_summary$r.squared,
    adj_r2 = sat_summary$adj.r.squared,
    oos_rmse = oos$rmse,
    oos_r2 = oos$r2,
    shock_scale_2001 = shock_2001$shock_scale,
    peak_dpd_1sd = pk1$peak,
    peak_h_1sd = pk1$horizon,
    peak_dpd_2001 = pk2$peak,
    peak_h_2001 = pk2$horizon
  )
  
  list(
    information_set = info_set_name,
    info_set_label = info_set_label,
    is_main_specification = identical(info_set_name, APP_BASELINE_INFO_SET),
    merton = list(p_ttc = APP_P_TTC, rho = APP_RHO),
    satellite = list(model = sat_model, bayesian = bsat, oos = oos),
    bands = list(
      macro_1sd = macro_1sd$bands,
      macro_2001 = macro_2001$bands,
      psiZ_1sd = psiZ_1sd$bands,
      psiZ_2001 = psiZ_2001$bands,
      pd_1sd = pd_1sd$bands,
      pd_2001 = pd_2001$bands
    ),
    decomposition = list(
      pd_summary_1sd = pd_dec_1sd,
      pd_summary_2001 = pd_dec_2001
    ),
    summary = summary_row,
    shock_2001 = shock_2001
  )
}

# -----------------------------------------------------------------------------
# 3. Run all selected information sets
# -----------------------------------------------------------------------------

runs <- lapply(names(INFO_SETS_TO_COMPARE), function(s) {
  run_transmission_for_information_set(s, INFO_SETS_TO_COMPARE[[s]])
})

names(runs) <- names(INFO_SETS_TO_COMPARE)

collect_bands <- function(key) {
  rbindlist(
    lapply(runs, function(r) r$bands[[key]]),
    use.names = TRUE,
    fill = TRUE
  )
}

collect_decomp <- function(key) {
  rbindlist(
    lapply(runs, function(r) r$decomposition[[key]]),
    use.names = TRUE,
    fill = TRUE
  )
}

macro_all <- rbindlist(
  list(collect_bands("macro_1sd"), collect_bands("macro_2001")),
  use.names = TRUE,
  fill = TRUE
)

z_all <- rbindlist(
  list(collect_bands("psiZ_1sd"), collect_bands("psiZ_2001")),
  use.names = TRUE,
  fill = TRUE
)

pd_all <- rbindlist(
  list(collect_bands("pd_1sd"), collect_bands("pd_2001")),
  use.names = TRUE,
  fill = TRUE
)

pd_decomp_all <- rbindlist(
  list(collect_decomp("pd_summary_1sd"), collect_decomp("pd_summary_2001")),
  use.names = TRUE,
  fill = TRUE
)

fwrite(
  macro_all,
  file.path(robust_out_dir, "macro_girf_bands_all_information_sets_all_shocks.csv")
)

fwrite(
  z_all,
  file.path(robust_out_dir, "z_girf_bands_all_information_sets_all_shocks.csv")
)

fwrite(
  pd_all,
  file.path(robust_out_dir, "pd_girf_bands_all_information_sets_all_shocks.csv")
)

fwrite(
  pd_decomp_all,
  file.path(robust_out_dir, "pd_variance_decomposition_90pct_all_information_sets.csv")
)

save_png_plot(
  plot_overlay(pd_all[shock == "one_sd"], expression(Delta ~ PD)),
  robust_out_dir,
  "pd_girf_overlay_1sd",
  width = 7,
  height = 4.4
)

save_png_plot(
  plot_overlay(pd_all[shock == "2001Q3"], expression(Delta ~ PD)),
  robust_out_dir,
  "pd_girf_overlay_2001Q3",
  width = 7,
  height = 4.4
)

save_png_plot(
  plot_overlay(z_all[shock == "one_sd"], expression(psi[Z](h))),
  robust_out_dir,
  "z_girf_overlay_1sd",
  width = 7,
  height = 4.4
)

save_png_plot(
  plot_overlay(z_all[shock == "2001Q3"], expression(psi[Z](h))),
  robust_out_dir,
  "z_girf_overlay_2001Q3",
  width = 7,
  height = 4.4
)

save_png_plot(
  plot_pd_bands_all_information_sets(pd_all),
  robust_out_dir,
  "pd_girf_bands_all_information_sets_all_shocks",
  width = 11,
  height = 6.5
)

save_png_plot(
  plot_macro_bands_all_information_sets(macro_all),
  robust_out_dir,
  "macro_girf_medians_all_information_sets_all_shocks",
  width = 14.5,
  height = max(7.5, 2.6 * length(selected_info_sets))
)

comparison <- rbindlist(
  lapply(runs, function(r) r$summary),
  use.names = TRUE,
  fill = TRUE
)

comparison[, information_set := factor(information_set, levels = selected_info_sets)]
setorder(comparison, information_set)
comparison[, information_set := as.character(information_set)]

fwrite(comparison, file.path(robust_out_dir, "information_set_comparison.csv"))

write_information_set_comparison_tex(
  comparison,
  file.path(robust_out_dir, "information_set_comparison.tex")
)

saveRDS(
  list(
    application = APP_NAME,
    label = APP_LABEL,
    note = if (exists("APP_NOTE")) APP_NOTE else NA_character_,
    merton = list(p_ttc = APP_P_TTC, rho = APP_RHO),
    information_sets = INFO_SETS_TO_COMPARE,
    info_set_labels = info_set_labels,
    baseline = APP_BASELINE_INFO_SET,
    runs = runs,
    comparison = comparison,
    macro_all = macro_all,
    z_all = z_all,
    pd_all = pd_all,
    pd_decomp_all = pd_decomp_all
  ),
  file.path(robust_out_dir, paste0(APP_NAME, "_information_set_robustness.rds"))
)

cat("\n--- ", APP_LABEL, " information-set comparison ---\n", sep = "")

print(
  comparison[
    ,
    .(
      information_set,
      label = info_set_label,
      main = is_main_specification,
      k,
      r2 = round(r2, 3),
      oos_rmse = round(oos_rmse, 4),
      peak_dpd_2001 = signif(peak_dpd_2001, 3),
      peak_h_2001
    )
  ]
)

cat("\n", APP_LABEL, " main-specification outputs written to: ",
    main_out_dir, "\n", sep = "")
cat(APP_LABEL, " robustness outputs written to: ",
    robust_out_dir, "\n", sep = "")