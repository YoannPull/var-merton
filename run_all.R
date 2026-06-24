# =============================================================================
#  run_all.R — Single entry point for the VAR-MERTON replication package
#
#  Run from the repository root:
#    source("run_all.R")
#  or:
#    Rscript run_all.R
#
#  Runs the full pipeline: data -> BVAR application + robustness -> direct-channel
#  test -> short-default-sample -> perfect foresight -> state dependence ->
#  paper deliverables.
#
#  Output layout:
#    output/applications/dralacbn/main_specification/
#    output/robustness/dralacbn/information_set/
#    output/shared/
#    output/paper_replication/   (final paper figures, tables, numbers)
# =============================================================================

t_start <- Sys.time()

# -----------------------------------------------------------------------------
# Execution options
# -----------------------------------------------------------------------------

STOP_ON_ERROR <- TRUE

# -----------------------------------------------------------------------------
# Load configuration
# -----------------------------------------------------------------------------

if (!file.exists("config.R")) {
  stop("Missing config.R. Please run this script from the repository root.")
}

source("config.R")

dir.create(if (exists("DIR_LOGS")) DIR_LOGS else "logs", showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Session log
# -----------------------------------------------------------------------------

run_stamp  <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_file   <- file.path(DIR_LOGS, paste0("run_all_", run_stamp, ".log"))
recap_file <- file.path(DIR_LOGS, paste0("run_all_", run_stamp, "_recap.csv"))

writeLines(capture.output(sessionInfo()), log_file)
cat("Session log written to: ", log_file, "\n\n", sep = "")

# -----------------------------------------------------------------------------
# Pipeline definition
# -----------------------------------------------------------------------------

#  Pipeline producing every figure, table and in-text number of the paper. It
#  includes the direct-channel / orthogonality test (03 -> 04, which relaxes the
#  satellite exclusion restriction and justifies the baseline closed form) and
#  the per-window BMA short-default-sample figures (07 -> 09).
pipeline_steps <- data.frame(
  phase = c(
    "data",
    "data",
    "shared",
    "dralacbn_application_and_robustness",
    "extension_direct_channel",
    "extension_direct_channel",
    "extension_short_sample",
    "extension_short_sample",
    "extension_perfect_foresight",
    "extension_state_dependence",
    "documentation",
    "paper_collation",
    "paper_replication"
  ),
  path = c(
    "code/01_process_gpr.R",
    "code/02_build_dataset.R",
    "code/shared/var_information_set_diagnostics.R",
    "code/dralacbn/02_dralacbn_information_set_robustness.R",
    "code/dralacbn/03_dralacbn_direct_channel.R",
    "code/dralacbn/04_dralacbn_direct_channel_paper_outputs.R",
    "code/dralacbn/07_dralacbn_short_sample_bma.R",
    "code/dralacbn/09_dralacbn_short_sample_single_figures.R",
    "code/dralacbn/08_dralacbn_perfect_foresight.R",
    "code/dralacbn/10_dralacbn_state_dependence.R",
    "code/shared/variable_definitions.R",
    "code/shared/collate_paper_outputs.R",
    "code/shared/make_paper_outputs.R"
  ),
  stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

format_duration <- function(seconds) {
  seconds <- max(0, as.numeric(seconds))
  h <- floor(seconds / 3600)
  m <- floor((seconds - 3600 * h) / 60)
  s <- round(seconds - 3600 * h - 60 * m)
  sprintf("%02d:%02d:%02d", h, m, s)
}

append_to_log <- function(...) {
  cat(paste0(..., collapse = ""), "\n", file = log_file, append = TRUE)
}

run_step <- function(path, phase) {
  cat("\n", strrep("=", 78), "\n", sep = "")
  cat(">>> ", phase, " | ", path, "  (", format(Sys.time(), "%H:%M:%S"), ")\n", sep = "")
  cat(strrep("=", 78), "\n", sep = "")

  append_to_log("")
  append_to_log(strrep("=", 78))
  append_to_log(">>> ", phase, " | ", path, "  (", format(Sys.time(), "%H:%M:%S"), ")")
  append_to_log(strrep("=", 78))

  t0 <- Sys.time()
  warning_messages <- character(0)
  error_message <- ""

  if (!file.exists(path)) {
    error_message <- paste0("File not found: ", path)
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    return(data.frame(
      phase = phase, path = path, ok = FALSE,
      elapsed_sec = elapsed, elapsed_hms = format_duration(elapsed),
      n_warnings = 0L, warnings = "", error = error_message,
      stringsAsFactors = FALSE
    ))
  }

  ok <- tryCatch(
    {
      withCallingHandlers(
        { source(path, local = .GlobalEnv) },
        warning = function(w) {
          msg <- conditionMessage(w)
          warning_messages <<- c(warning_messages, msg)
          message("WARNING in ", path, ": ", msg)
          append_to_log("WARNING in ", path, ": ", msg)
          invokeRestart("muffleWarning")
        }
      )
      TRUE
    },
    error = function(e) {
      error_message <<- conditionMessage(e)
      message("ERROR in ", path, ": ", error_message)
      append_to_log("ERROR in ", path, ": ", error_message)
      FALSE
    }
  )

  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  status <- if (ok) "OK" else "FAILED"
  cat(sprintf("<<< %s finished in %.1f s — %s\n", path, elapsed, status))
  append_to_log(sprintf("<<< %s finished in %.1f s — %s", path, elapsed, status))

  data.frame(
    phase = phase,
    path = path,
    ok = isTRUE(ok),
    elapsed_sec = elapsed,
    elapsed_hms = format_duration(elapsed),
    n_warnings = length(warning_messages),
    warnings = paste(unique(warning_messages), collapse = " | "),
    error = error_message,
    stringsAsFactors = FALSE
  )
}

print_recap <- function(results_df) {
  results_df$ok <- as.logical(results_df$ok)
  results_df$ok[is.na(results_df$ok)] <- FALSE

  cat("\n", strrep("=", 78), "\n", sep = "")
  cat("RUN_ALL RECAP\n")
  cat(strrep("=", 78), "\n", sep = "")

  for (i in seq_len(nrow(results_df))) {
    status <- if (isTRUE(results_df$ok[i])) "OK " else "KO "
    cat(sprintf("  [%s] %-36s %s\n", status, results_df$phase[i], results_df$path[i]))
    if (!isTRUE(results_df$ok[i]) && nzchar(results_df$error[i])) {
      cat("       Error: ", results_df$error[i], "\n", sep = "")
    }
    if (!is.na(results_df$n_warnings[i]) && results_df$n_warnings[i] > 0) {
      cat("       Warnings: ", results_df$n_warnings[i], "\n", sep = "")
    }
  }

  n_ok <- sum(results_df$ok)
  n_total <- nrow(results_df)
  duration_min <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))

  cat("\n")
  cat(sprintf("Status: %d/%d steps OK, %d KO.\n", n_ok, n_total, n_total - n_ok))
  cat(sprintf("Total duration: %.1f min\n", duration_min))
  cat("Recap CSV written to: ", recap_file, "\n", sep = "")
  cat("Session log written to: ", log_file, "\n", sep = "")

  cat("\nMain application outputs:\n")
  cat("  - ", file.path("output", "applications", "dralacbn", "main_specification"), "\n", sep = "")
  cat("Robustness outputs:\n")
  cat("  - ", file.path("output", "robustness", "dralacbn", "information_set"), "\n", sep = "")
  cat("Paper deliverables:\n")
  cat("  - ", file.path("output", "paper_replication"), "\n", sep = "")
}

# -----------------------------------------------------------------------------
# Run pipeline
# -----------------------------------------------------------------------------

pipeline_results <- vector("list", nrow(pipeline_steps))

for (i in seq_len(nrow(pipeline_steps))) {
  pipeline_results[[i]] <- run_step(
    pipeline_steps$path[i],
    pipeline_steps$phase[i]
  )
  
  if (!isTRUE(pipeline_results[[i]]$ok[1]) && isTRUE(STOP_ON_ERROR)) {
    cat("\nStopping run_all because STOP_ON_ERROR = TRUE.\n")
    break
  }
}

pipeline_results <- pipeline_results[
  !vapply(pipeline_results, is.null, logical(1))
]

if (length(pipeline_results) == 0L) {
  stop("run_all.R: no pipeline results to summarize.")
}

results_df <- data.table::rbindlist(
  pipeline_results,
  use.names = TRUE,
  fill = TRUE
)

results_df <- as.data.frame(results_df)

write.csv(results_df, recap_file, row.names = FALSE)
print_recap(results_df)

if (any(!results_df$ok)) {
  cat("\nAt least one step failed. Check the recap above and the log file.\n")
  if (!interactive()) quit(status = 1)
} else {
  cat("\nAll requested steps completed successfully.\n")
}