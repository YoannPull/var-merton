# =============================================================================
#  tests/smoke_test.R — Lightweight post-run sanity checks.
#
#  Run from the repository root AFTER run_all.R:
#     Rscript tests/smoke_test.R
#
#  It does not re-estimate anything; it verifies that the pipeline produced the
#  expected article deliverables with sane shapes (columns present, no NA in the
#  median response, correct horizon length). Exits with status 1 if any check
#  fails, so it can gate a replication run in CI.
# =============================================================================

if (!file.exists("config.R")) {
  stop("smoke_test.R: run from the repository root.")
}
source("config.R")

n_checks <- 0L
n_fail <- 0L

chk <- function(cond, msg) {
  n_checks <<- n_checks + 1L
  if (isTRUE(cond)) {
    cat("  ok   ", msg, "\n", sep = "")
  } else {
    n_fail <<- n_fail + 1L
    cat("  FAIL ", msg, "\n", sep = "")
  }
}

H <- if (exists("HORIZON")) HORIZON else 12L
expected_band_cols <- c("horizon", "lower90", "lower68", "median",
                        "upper68", "upper90")

# The DRALACBN application is the one checked.
apps <- list(
  dralacbn = if (exists("DRALACBN_INFO_SETS_TO_COMPARE")) DRALACBN_INFO_SETS_TO_COMPARE else "real_side"
)

cat("== Paper deliverables ==\n")
chk(file.exists(file.path(DIR_PAPER, "INDEX.csv")), "paper/INDEX.csv exists")

check_bands <- function(path, label) {
  if (!file.exists(path)) { chk(FALSE, paste0(label, " exists")); return(invisible()) }
  d <- tryCatch(read.csv(path), error = function(e) NULL)
  chk(!is.null(d), paste0(label, " readable"))
  if (is.null(d)) return(invisible())
  chk(all(expected_band_cols %in% names(d)), paste0(label, " has band columns"))
  chk(nrow(d) == H + 1L, paste0(label, " has H+1 = ", H + 1L, " rows"))
  chk(all(is.finite(d$median)), paste0(label, " median has no NA/Inf"))
}

for (app in names(apps)) {
  cat("== ", app, " ==\n", sep = "")
  fig_dir <- file.path(DIR_PAPER_FIG, app)
  dat_dir <- file.path(DIR_PAPER_DATA, app)
  tab_dir <- file.path(DIR_PAPER_TAB, app)

  chk(dir.exists(fig_dir) && length(list.files(fig_dir)) > 0,
      paste0(app, ": figures present"))
  chk(file.exists(file.path(tab_dir, paste0("bma_", app, "_real_side.tex"))),
      paste0(app, ": main-spec BMA table present"))

  for (spec in apps[[app]]) {
    tag <- paste0(app, "_", spec)
    check_bands(file.path(dat_dir, paste0("pd_mean_bands_", tag, "_1sd.csv")),
                paste0(tag, " mean PD bands (1sd)"))
    pv <- file.path(dat_dir, paste0("pd_var_peaks_", tag, ".csv"))
    if (file.exists(pv)) {
      d <- tryCatch(read.csv(pv), error = function(e) NULL)
      chk(!is.null(d) && nrow(d) > 0, paste0(tag, " PD-VaR peaks non-empty"))
    } else {
      chk(FALSE, paste0(tag, " PD-VaR peaks exist"))
    }
  }
}

# -----------------------------------------------------------------------------
#  Paper replication package: every figure the LaTeX source includes must be
#  present, under the exact name the .tex references. This list mirrors the
#  \includegraphics calls of the manuscript; keep the two in sync.
# -----------------------------------------------------------------------------

cat("== paper_replication figures (names referenced by the .tex) ==\n")

REP_FIG <- file.path("output", "paper_replication", "figures")

tex_figures <- c(
  "main/graph_dralacbn.pdf",
  "main/graph_e1.pdf",
  "main/graph_gpr.pdf",
  "main/graph_pd.pdf",
  "main/graph_var_irf.pdf",
  "main/pd_girf_window_2005_1sd.pdf",
  "main/pd_girf_window_2015_1sd.pdf",
  "main/pd_mean_vs_es_dralacbn_real_side_1sd.pdf",
  "main/pd_mean_vs_var_dralacbn_real_side_1sd.pdf",
  "main/simulation_benchmark.pdf",
  "robustness/monetary/pd_girf_1sd.pdf",
  "robustness/monetary/var_irf_1sd.pdf",
  "robustness/real_side_light/pd_girf_1sd.pdf",
  "robustness/real_side_light/var_irf_1sd.pdf",
  "robustness/uncertainty/pd_girf_1sd.pdf",
  "robustness/uncertainty/var_irf_1sd.pdf"
)

for (f in tex_figures) {
  p <- file.path(REP_FIG, f)
  chk(file.exists(p) && file.info(p)$size > 0, paste0("figure ", f))
}

chk(file.exists(file.path("output", "paper_replication", "paper_numbers.csv")),
    "paper_numbers.csv exists")

cat("\n", strrep("=", 50), "\n", sep = "")
cat(sprintf("Smoke test: %d/%d checks passed, %d failed.\n",
            n_checks - n_fail, n_checks, n_fail))

if (n_fail > 0L) {
  if (!interactive()) quit(status = 1L)
} else {
  cat("All smoke checks passed.\n")
}
