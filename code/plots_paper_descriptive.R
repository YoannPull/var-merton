# =============================================================================
#  code/plots_paper_descriptive.R — All paper figures, in vector PDF + PNG.
#
#  Run from the repository root, AFTER run_all.R:
#     source("code/plots_paper_descriptive.R")
#
#  Output (split by format then series):
#     output/paper/pdf/<series>/<name>.pdf
#     output/paper/png/<series>/<name>.png
#
#  Series:
#     gpr/      graph_gpr (with events), graph_e1 (with events)
#     eba/      graph_eba (EBA corporate default rate)
#     dralacbn/ graph_dralacbn (delinquency, NBER bands)
#     results/  graph_var_irf, graph_z, graph_pd (one-s.d. shock)
#               graph_z_scenario, graph_pd_scenario (2001Q3 shock)
#
#  Result figures are re-rendered from the main-specification band CSVs, so they
#  are vector PDFs in the same clean (grid-free) theme as the descriptive plots.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(scales)
})

if (file.exists("config.R")) suppressWarnings(try(source("config.R"), silent = TRUE))

paper_dir <- if (exists("DIR_PAPER")) DIR_PAPER else file.path("output", "paper")
pdf_base <- file.path(paper_dir, "pdf")
png_base <- file.path(paper_dir, "png")

main_dir <- if (exists("DIR_DRALACBN_MAIN")) DIR_DRALACBN_MAIN else
  file.path("output", "applications", "dralacbn", "main_specification")

BURGUNDY   <- if (exists("SQUARE_BURGUNDY2")) SQUARE_BURGUNDY2 else "#6F1732"
ROSE       <- if (exists("SQUARE_ROSE2")) SQUARE_ROSE2 else "#AB4A7D"
ROSE_LIGHT <- if (exists("SQUARE_ROSE")) SQUARE_ROSE else "#F3D6E3"

DATE_MAX <- as.Date("2024-12-31")
GPR_DATE_MIN <- as.Date("1986-01-01")

theme_paper <- function() {
  theme_minimal(base_size = 13) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_text(color = "grey20"),
      strip.text = element_text(face = "bold", color = "grey20"),
      plot.margin = margin(6, 12, 6, 6)
    )
}

save_fig <- function(plot, series, stem, width = 9, height = 4.2, dpi = 300) {
  pdir <- file.path(pdf_base, series)
  ndir <- file.path(png_base, series)
  dir.create(pdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(ndir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(pdir, paste0(stem, ".pdf")), plot, device = grDevices::pdf,
         width = width, height = height, useDingbats = FALSE)
  ggsave(file.path(ndir, paste0(stem, ".png")), plot,
         width = width, height = height, dpi = dpi)
  cat("  wrote ", series, "/", stem, ".{pdf,png}\n", sep = "")
}

make_quarter <- function(date) {
  paste0(format(date, "%Y"), "Q", ((as.integer(format(date, "%m")) - 1L) %/% 3L) + 1L)
}

events_dt <- data.table(
  quarter = c("1990Q3", "1991Q1", "2001Q3", "2003Q1", "2005Q3", "2006Q3",
              "2014Q3", "2022Q1", "2023Q4"),
  label = c("Kuwait invasion", "Gulf War", "11 September", "Iraq invasion",
            "London bombings", "Transatlantic aircraft plot",
            "Ukraine/Russia escalation", "Russia invades Ukraine",
            "Israel--Hamas war")
)

recessions <- data.table(
  start = as.Date(c("1990-07-01", "2001-03-01", "2007-12-01", "2020-02-01")),
  end   = as.Date(c("1991-03-31", "2001-11-30", "2009-06-30", "2020-04-30"))
)

var_labels <- c(
  log_GPRD       = "Geopolitical Risk (log)",
  log_inv_pc     = "Investment p.c. (log)",
  log_gdp_pc     = "Real GDP p.c. (log)",
  log_private_pc = "Private employment p.c. (log)",
  log_oil_real   = "Real oil price (log)",
  infl_yoy_pct   = "Inflation YoY (%)"
)

add_events <- function(p, series_dt, value_col, label_seed = 1) {
  d <- as.data.table(copy(series_dt))
  d[, quarter := make_quarter(date)]
  ev <- merge(events_dt, unique(d[, .(quarter, date)]), by = "quarter")
  ev <- ev[!is.na(date)]
  yr <- range(d[[value_col]], na.rm = TRUE)
  off <- diff(yr) * 0.03
  pos <- merge(ev, d[, .(quarter, val = get(value_col))], by = "quarter")
  pos[, y := pmin(val + off, yr[2])]
  p +
    geom_vline(data = ev, aes(xintercept = date), linetype = "dashed",
               alpha = 0.5, color = ROSE) +
    geom_label_repel(data = pos, aes(x = date, y = y, label = label),
                     inherit.aes = FALSE, direction = "y", size = 3.0,
                     fill = "white", label.size = 0, box.padding = 0.2,
                     point.padding = 0.15, min.segment.length = 0,
                     seed = label_seed, segment.color = ROSE,
                     segment.size = 0.2, segment.alpha = 0.5)
}

# Single-panel response band plot (median + 68/90% ribbons), grid-free.
plot_band <- function(path, ylab_txt, scale = 1) {
  if (!file.exists(path)) return(NULL)
  b <- as.data.table(fread(path))
  cols <- c("horizon", "lower90", "lower68", "median", "upper68", "upper90")
  if (!all(cols %in% names(b))) return(NULL)
  for (cc in c("lower90", "lower68", "median", "upper68", "upper90")) {
    b[[cc]] <- b[[cc]] * scale
  }
  ggplot(b, aes(horizon, median)) +
    geom_ribbon(aes(ymin = lower90, ymax = upper90), fill = ROSE_LIGHT, alpha = 0.65) +
    geom_ribbon(aes(ymin = lower68, ymax = upper68), fill = ROSE, alpha = 0.42) +
    geom_line(color = BURGUNDY, linewidth = 1.05) +
    geom_hline(yintercept = 0, color = "grey55", linetype = "dashed", linewidth = 0.4) +
    scale_x_continuous(breaks = function(x) seq(max(0, ceiling(x[1])), floor(x[2]), by = 1),
                       minor_breaks = NULL) +
    labs(x = "Horizon (quarters)", y = ylab_txt) +
    theme_paper()
}

# Faceted macro IRF (one panel per variable).
plot_band_facet <- function(path) {
  if (!file.exists(path)) return(NULL)
  b <- as.data.table(fread(path))
  if (!"variable" %in% names(b)) return(NULL)
  b[, vlab := var_labels[as.character(variable)]]
  b[is.na(vlab), vlab := as.character(variable)]
  ggplot(b, aes(horizon, median)) +
    geom_ribbon(aes(ymin = lower90, ymax = upper90), fill = ROSE_LIGHT, alpha = 0.65) +
    geom_ribbon(aes(ymin = lower68, ymax = upper68), fill = ROSE, alpha = 0.42) +
    geom_line(color = BURGUNDY, linewidth = 0.95) +
    geom_hline(yintercept = 0, color = "grey55", linetype = "dashed", linewidth = 0.35) +
    facet_wrap(~ vlab, scales = "free_y") +
    scale_x_continuous(breaks = function(x) seq(max(0, ceiling(x[1])), floor(x[2]), by = 2),
                       minor_breaks = NULL) +
    labs(x = "Horizon (quarters)", y = "Response") +
    theme_paper()
}

# -----------------------------------------------------------------------------
# Descriptive figures
# -----------------------------------------------------------------------------
gpr_path <- if (exists("PATH_DATA_VAR")) PATH_DATA_VAR else "data/processed/data_var_for_model.csv"
if (file.exists(gpr_path)) {
  g <- fread(gpr_path)
  g[, date := as.Date(Date)]
  g <- g[date >= GPR_DATE_MIN & date <= DATE_MAX & is.finite(GPRD)]
  setorder(g, date)
  p_gpr <- ggplot(g, aes(date, GPRD)) +
    geom_line(color = BURGUNDY, linewidth = 0.7) +
    labs(x = NULL, y = "Geopolitical Risk Index") +
    scale_x_date(date_breaks = "5 years", date_labels = "%Y") + theme_paper()
  save_fig(add_events(p_gpr, g[, .(date, GPRD)], "GPRD", 1), "gpr", "graph_gpr")
} else cat("  skip graph_gpr: ", gpr_path, " not found.\n", sep = "")

e1_path <- file.path(main_dir, "irf", "VAR", "standardized_gpr_innovation.csv")
if (file.exists(e1_path)) {
  e <- fread(e1_path); e[, date := as.Date(date)]
  p_e1 <- ggplot(e, aes(date, e1)) +
    geom_hline(yintercept = 0, color = "grey55", linewidth = 0.4) +
    geom_line(color = BURGUNDY, linewidth = 0.6) +
    labs(x = NULL, y = expression(epsilon[GPR][","][t] ~ "(s.d. units)")) +
    scale_x_date(date_breaks = "5 years", date_labels = "%Y") + theme_paper()
  save_fig(add_events(p_e1, e[, .(date, e1)], "e1", 123), "gpr", "graph_e1")
} else cat("  skip graph_e1: run run_all.R first.\n")

eba_path <- if (exists("PATH_RAW_RISK_EBA")) PATH_RAW_RISK_EBA else "data/raw/data_risk_EBA.csv"
if (file.exists(eba_path)) {
  r <- fread(eba_path, sep = ";")
  r <- r[Country == "United States"]
  r <- data.table(date = as.Date(r$Date, format = "%d/%m/%Y"),
                  rate = as.numeric(r$Corpo_DR_WA) * 100)
  r <- r[is.finite(rate)]; setorder(r, date); ttc <- mean(r$rate, na.rm = TRUE)
  p_eba <- ggplot(r, aes(date, rate)) +
    geom_hline(yintercept = ttc, color = "grey55", linetype = "dashed", linewidth = 0.4) +
    geom_line(color = BURGUNDY, linewidth = 0.8) +
    annotate("text", x = max(r$date), y = ttc, label = sprintf("TTC = %.2f%%", ttc),
             hjust = 1.0, vjust = -0.6, size = 3.4, color = "grey35") +
    labs(x = NULL, y = "Corporate default rate (%)") +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") + theme_paper()
  save_fig(p_eba, "eba", "graph_eba")
} else cat("  skip graph_eba: ", eba_path, " not found.\n", sep = "")

del_path <- file.path(if (exists("DIR_RAW")) DIR_RAW else "data/raw", "default", "DRALACBN.csv")
if (file.exists(del_path)) {
  dd <- fread(del_path)
  dcol <- intersect(c("observation_date", "DATE", "Date", "date"), names(dd))[1]
  vcol <- intersect(c("DRALACBN", "value", "VALUE"), names(dd))[1]
  dd <- data.table(date = as.Date(dd[[dcol]]), rate = as.numeric(dd[[vcol]]))
  dd <- dd[date <= DATE_MAX & is.finite(rate)]; setorder(dd, date)
  ttc <- mean(dd$rate, na.rm = TRUE)
  rec <- recessions[end >= min(dd$date) & start <= max(dd$date)]
  p_dra <- ggplot(dd, aes(date, rate)) +
    geom_rect(data = rec, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = "grey75", alpha = 0.35) +
    geom_hline(yintercept = ttc, color = "grey55", linetype = "dashed", linewidth = 0.4) +
    geom_line(color = BURGUNDY, linewidth = 0.8) +
    annotate("text", x = max(dd$date), y = ttc, label = sprintf("TTC = %.2f%%", ttc),
             hjust = 1.0, vjust = -0.6, size = 3.4, color = "grey35") +
    labs(x = NULL, y = "Delinquency rate (%)") +
    scale_x_date(date_breaks = "5 years", date_labels = "%Y") + theme_paper()
  save_fig(p_dra, "dralacbn", "graph_dralacbn")
} else cat("  skip graph_dralacbn: ", del_path, " not found.\n", sep = "")

# -----------------------------------------------------------------------------
# Result figures (main specification), re-rendered from band CSVs
# -----------------------------------------------------------------------------
render_result <- function(rel_path, stem, ylab_txt, scale = 1, facet = FALSE,
                          width = 9, height = 4.2) {
  path <- file.path(main_dir, rel_path)
  p <- if (facet) plot_band_facet(path) else plot_band(path, ylab_txt, scale)
  if (is.null(p)) {
    cat("  skip results/", stem, " (", path, " missing).\n", sep = "")
  } else {
    save_fig(p, "results", stem, width = width, height = height)
  }
}

render_result(file.path("irf", "VAR", "macro_girf_bands_1sd.csv"), "graph_var_irf",
              NULL, facet = TRUE, width = 9, height = 5.4)
render_result(file.path("irf", "Z", "psiZ_bands_1sd.csv"), "graph_z",
              "Systematic factor Z (s.d.)")
render_result(file.path("irf", "PD", "pd_girf_bands_1sd.csv"), "graph_pd",
              expression(Delta ~ "PD (pp)"), scale = 100)
render_result(file.path("irf", "Z", "psiZ_bands_2001Q3.csv"), "graph_z_scenario",
              "Systematic factor Z (s.d.)")
render_result(file.path("irf", "PD", "pd_girf_bands_2001Q3.csv"), "graph_pd_scenario",
              expression(Delta ~ "PD (pp)"), scale = 100)

cat("Paper figures (PDF + PNG) written under ", paper_dir, "/{pdf,png}/\n", sep = "")
