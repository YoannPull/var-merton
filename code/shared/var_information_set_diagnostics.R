# ----------------------------------------------------------------------------
# Run from the ROOT of the repository.
source("code/00_setup.R")
source("code/shared/_helpers_var_merton.R")
# ----------------------------------------------------------------------------

# =============================================================================
# Robustness — VAR information-set stability diagnostics
#
# For each retained macro-financial information set defined in config.R:
#   - estimate a Bayesian VAR with a Normal-Inverse-Wishart posterior;
#   - retain dynamically stable posterior draws;
#   - report stable-share and companion-root diagnostics.
#
# Kernels are saved to a shared directory so that downstream robustness scripts
# can reuse them without re-estimating the BVARs.
#
# The retained information sets are defined in config.R through VAR_SPECS.
# For example:
#   - real_side
#   - rates_activity
#   - caldara_style
#
# Outputs:
#   output/robustness/var_kernels/<information_set>_var_kernel.rds
#   output/robustness/var_information_sets/var_stability_diagnostics.csv
#   output/robustness/var_information_sets/var_stability_diagnostics.tex
#   output/robustness/var_information_sets/var_information_sets_list.csv
#   output/robustness/var_information_sets/var_information_set_kernels_index.rds
# =============================================================================

# =============================================================================
# 0. Paths
# =============================================================================

data_var_path <- if (exists("PATH_DATA_VAR")) {
  PATH_DATA_VAR
} else {
  "data/processed/data_var_for_model.csv"
}

out_dir <- if (exists("DIR_VAR_SPECS")) {
  DIR_VAR_SPECS
} else {
  "output/robustness/var_information_sets"
}

kernel_dir <- if (exists("DIR_VAR_KERNELS")) {
  DIR_VAR_KERNELS
} else {
  "output/robustness/var_kernels"
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(kernel_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. BVAR parameters
# =============================================================================

p_lags <- if (exists("P_LAGS")) {
  P_LAGS
} else {
  2
}

nrep <- if (exists("NREP")) {
  NREP
} else {
  30000
}

M_target <- if (exists("M_TARGET")) {
  M_TARGET
} else {
  10000L
}

H <- if (exists("HORIZON")) {
  HORIZON
} else {
  12
}

seed <- if (exists("SEED_ROBUSTNESS")) {
  SEED_ROBUSTNESS
} else if (exists("SEED")) {
  SEED
} else {
  123
}

impulse_name <- if (exists("IMPULSE_NAME")) {
  IMPULSE_NAME
} else {
  "log_GPRD"
}

set.seed(seed)

# =============================================================================
# 2. VAR information sets
# =============================================================================
#
# The canonical information sets live in config.R.
# Each vector lists variables included IN ADDITION to the impulse variable.
# The impulse variable is added first by the BVAR estimation helper.
#
# Because the downstream analysis uses generalized impulse responses, the
# within-list order is not a structural identifying assumption. It is fixed only
# for reproducibility and plotting.

if (!exists("VAR_SPECS")) {
  stop(
    "var_information_set_diagnostics.R: VAR_SPECS not found. ",
    "Source config.R first."
  )
}

var_info_sets <- VAR_SPECS

if (length(var_info_sets) == 0L) {
  stop("var_information_set_diagnostics.R: VAR_SPECS is empty.")
}

# Optional human-readable labels. If config.R defines INFO_SET_LABELS, use it.
info_set_labels <- if (exists("INFO_SET_LABELS")) {
  INFO_SET_LABELS
} else {
  c(
    real_side = "Real-side",
    rates_activity = "Rates/activity",
    caldara_style = "Caldara-style"
  )
}

# Make sure every information set has a printable label.
for (nm in names(var_info_sets)) {
  if (!(nm %in% names(info_set_labels))) {
    info_set_labels[nm] <- nm
  }
}

cat("\n========== Robustness: VAR information-set diagnostics ==========\n")
cat("Information sets retained:\n")
for (nm in names(var_info_sets)) {
  cat(
    "  - ", nm,
    " (", info_set_labels[[nm]], "): ",
    paste(var_info_sets[[nm]], collapse = ", "),
    "\n",
    sep = ""
  )
}

cat("\nBVAR parameters:\n")
cat("  p_lags       = ", p_lags, "\n", sep = "")
cat("  nrep         = ", nrep, "\n", sep = "")
cat("  M_target     = ", M_target, "\n", sep = "")
cat("  H            = ", H, "\n", sep = "")
cat("  seed         = ", seed, "\n", sep = "")
cat("  impulse_name = ", impulse_name, "\n", sep = "")

# =============================================================================
# 3. Load macro data
# =============================================================================

DT_raw <- safe_fread(data_var_path)

if (!"Date" %in% names(DT_raw)) {
  stop("data_var_for_model.csv must contain Date.")
}

DT_raw[, Date := as.Date(Date)]

if (anyNA(DT_raw$Date)) {
  stop("Date conversion produced NA values in data_var_for_model.csv.")
}

# Quick check that the impulse is present.
if (!(impulse_name %in% names(DT_raw))) {
  stop(
    "Impulse variable '", impulse_name,
    "' not found in data_var_for_model.csv."
  )
}

# Check all information-set variables before starting estimation.
all_required_vars <- unique(c(impulse_name, unlist(var_info_sets, use.names = FALSE)))
missing_vars <- setdiff(all_required_vars, names(DT_raw))

if (length(missing_vars) > 0L) {
  stop(
    "The following variables required by VAR_SPECS are missing from ",
    "data_var_for_model.csv: ",
    paste(missing_vars, collapse = ", ")
  )
}

# =============================================================================
# 4. Estimate every information set
# =============================================================================

results <- list()
table_rows <- list()

for (info_set_name in names(var_info_sets)) {
  
  info_set_label <- info_set_labels[[info_set_name]]
  vars_extra <- var_info_sets[[info_set_name]]
  var_order <- c(impulse_name, vars_extra)
  
  cat("\n", strrep("-", 78), "\n", sep = "")
  cat("Estimating information set: ", info_set_name, " (", info_set_label, ")\n", sep = "")
  cat("Variables: ", paste(var_order, collapse = " -> "), "\n", sep = "")
  
  complete_n <- DT_raw[complete.cases(DT_raw[, ..var_order]), .N]
  
  cat("Complete-case observations: ", complete_n, "\n", sep = "")
  
  if (complete_n < 40L) {
    warning(
      "Information set '", info_set_name,
      "' has only ", complete_n,
      " complete-case observations. Check the input data."
    )
  }
  
  res_info_set <- estimate_bvar_kernel(
    DT_raw       = DT_raw,
    spec_name    = info_set_name,
    vars_extra   = vars_extra,
    p_lags       = p_lags,
    nrep         = nrep,
    H            = H,
    seed         = seed,
    M_target     = M_target,
    kernel_dir   = kernel_dir,
    impulse_name = impulse_name
  )
  
  results[[info_set_name]] <- res_info_set$kernel
  
  row_i <- as.data.table(res_info_set$table_row)
  
  # Add clean paper-facing labels and variable lists.
  row_i[, information_set := info_set_name]
  row_i[, information_set_label := info_set_label]
  row_i[, variables_extra := paste(vars_extra, collapse = ", ")]
  row_i[, variables_ordered := paste(var_order, collapse = ", ")]
  row_i[, complete_obs := complete_n]
  
  # If helper returns column `spec`, keep it for backward compatibility but
  # make sure it matches the new naming convention.
  if ("spec" %in% names(row_i)) {
    row_i[, spec := info_set_name]
  } else {
    row_i[, spec := info_set_name]
  }
  
  table_rows[[info_set_name]] <- row_i
  
  cat(
    "Stable share = ",
    sprintf("%.3f", row_i$stable_share[1]),
    " | median stable root = ",
    sprintf("%.3f", row_i$lambda_max_median_stable[1]),
    "\n",
    sep = ""
  )
}

stability <- rbindlist(table_rows, use.names = TRUE, fill = TRUE)

# Put clean columns first, while preserving helper-specific diagnostics.
front_cols <- intersect(
  c(
    "information_set",
    "information_set_label",
    "spec",
    "k",
    "k_minus_1",
    "p",
    "complete_obs",
    "variables_extra",
    "variables_ordered",
    "n_draws",
    "n_stable",
    "stable_share",
    "lambda_max_median_all",
    "lambda_max_p95_all",
    "lambda_max_median_stable",
    "lambda_max_p95_stable"
  ),
  names(stability)
)

other_cols <- setdiff(names(stability), front_cols)
setcolorder(stability, c(front_cols, other_cols))

# =============================================================================
# 5. Write stability diagnostics
# =============================================================================

fwrite(
  stability,
  file.path(out_dir, "var_stability_diagnostics.csv")
)

# Backward-compatible table writer.
write_stability_table_tex(
  stability,
  file.path(out_dir, "var_stability_diagnostics.tex"),
  caption = "Posterior stability diagnostics for the BVAR information sets",
  label   = "tab:var_stability_info_sets"
)

# Machine-readable list of information sets.
info_sets_list <- rbindlist(
  lapply(names(var_info_sets), function(s) {
    vars_extra <- var_info_sets[[s]]
    data.table(
      information_set = s,
      information_set_label = info_set_labels[[s]],
      k_minus_1 = length(vars_extra),
      variables_extra = paste(vars_extra, collapse = ", "),
      variables_ordered = paste(c(impulse_name, vars_extra), collapse = ", ")
    )
  }),
  use.names = TRUE,
  fill = TRUE
)

fwrite(
  info_sets_list,
  file.path(out_dir, "var_information_sets_list.csv")
)

# Backward-compatible old filename.
fwrite(
  info_sets_list,
  file.path(out_dir, "var_information_sets_list.csv")
)

saveRDS(
  results,
  file.path(out_dir, "var_information_set_kernels_index.rds")
)

# Backward-compatible old filename.
saveRDS(
  results,
  file.path(out_dir, "var_information set_kernels_index.rds")
)

# =============================================================================
# 6. Console summary
# =============================================================================

cat("\n", strrep("=", 78), "\n", sep = "")
cat("VAR information-set diagnostics summary\n")
cat(strrep("=", 78), "\n", sep = "")

print(
  stability[, .(
    information_set,
    information_set_label,
    k,
    complete_obs,
    n_draws,
    n_stable,
    stable_share = round(stable_share, 3),
    lambda_max_median_stable = round(lambda_max_median_stable, 4),
    lambda_max_p95_stable = round(lambda_max_p95_stable, 4)
  )]
)

cat("\nVAR information-set diagnostics written to: ", out_dir, "\n", sep = "")
cat("Kernels saved to: ", kernel_dir, "\n", sep = "")