# =============================================================================
#  code/dralacbn/10_dralacbn_state_dependence.R
#
#  State dependence of the default-probability response.
#
#  The PD GIRF of Proposition 3 is history-dependent through the conditional
#  mean of the systematic factor: the same innovation produces different PD
#  responses depending on the state of the credit cycle when it strikes. This
#  script makes that dependence explicit, using only geopolitical episodes:
#
#  (1) NAMED STATES, full posterior. The GIRF is evaluated conditional on the
#      histories preceding the three largest positive GPR innovations in the
#      sample (transparent selection rule, |e1| ranking):
#         1990:Q3  Gulf War (delinquency ~5.3%, stressed credit state)
#         2001:Q3  September 11 (delinquency ~2.7%, intermediate)
#         2022:Q1  Russia invades Ukraine (delinquency ~1.5%, benign)
#      plus the end-of-sample state (2024:Q4) used in the main results.
#      For each state: response to a common one-s.d. shock (isolates state
#      dependence) and to the episode's own historical shock size.
#
#  (2) CONVEXITY CURVE, posterior-median parameters. The peak PD response to
#      a one-s.d. shock is computed for EVERY admissible history in the
#      sample and plotted against the baseline expected PD at the shock date.
#      The named episodes are highlighted on the curve. This removes any
#      date-selection discretion.
#
#  Mechanics: the macro block is linear, so psi_Z and the conditional
#  variances are history-invariant; only the baseline mean path mu_{t+h}
#  changes with the conditioning history. Parameters are full-sample
#  estimates; the exercise is a counterfactual within the estimated model.
#
#  STANDALONE: requires the main DRALACBN application (02). Run from the root:
#      source("code/dralacbn/10_dralacbn_state_dependence.R")
#
#  Outputs: output/applications/dralacbn/state_dependence/
# =============================================================================

source("code/00_setup.R")
source("code/shared/_helpers_var_merton.R")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

# -----------------------------------------------------------------------------
# 0. Parameters, episodes, output tree
# -----------------------------------------------------------------------------

SD_INFO_SET  <- if (exists("DELINQUENCY_BASELINE_INFO_SET")) DELINQUENCY_BASELINE_INFO_SET else "real_side"
p_lags       <- if (exists("P_LAGS")) P_LAGS else 2
nrep         <- if (exists("NREP")) NREP else 30000
M_target     <- if (exists("M_TARGET")) M_TARGET else 10000L
H            <- if (exists("HORIZON")) HORIZON else 12
seed         <- if (exists("SEED_ROBUSTNESS")) SEED_ROBUSTNESS else 123
impulse_name <- if (exists("IMPULSE_NAME")) IMPULSE_NAME else "log_GPRD"

EPISODES <- list(
  list(quarter = "1990Q3", label = "Gulf War (1990:Q3)"),
  list(quarter = "2001Q3", label = "September 11 (2001:Q3)"),
  list(quarter = "2022Q1", label = "Ukraine invasion (2022:Q1)")
)
# CURRENT_Q / CURRENT_LABEL (end-of-sample state) are computed below, once the
# effective sample dates are known: the shock quarter is the one FOLLOWING the
# last observation of the kernel sample.

sd_out_dir <- file.path("output", "applications", "dralacbn", "state_dependence")
for (sub in c("", "tables", "figures", "irf")) {
  dir.create(file.path(sd_out_dir, sub), recursive = TRUE, showWarnings = FALSE)
}

theme_paper <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      axis.title = element_text(color = "#1A1A1A"),
      axis.text = element_text(color = "grey20"),
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linetype = "dotted", linewidth = 0.25),
      plot.caption = element_blank()
    )
}

save_fig <- function(plot, stem, width, height, dpi = 300) {
  ggsave(file.path(sd_out_dir, "figures", paste0(stem, ".png")),
         plot, width = width, height = height, dpi = dpi)
  tryCatch(
    ggsave(file.path(sd_out_dir, "figures", paste0(stem, ".pdf")),
           plot, width = width, height = height, device = grDevices::pdf),
    error = function(e) tryCatch(
      ggsave(file.path(sd_out_dir, "figures", paste0(stem, ".pdf")),
             plot, width = width, height = height),
      error = function(e2) message("PDF export failed for ", stem)
    )
  )
}

fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)

cat("\n========== DRALACBN: state dependence of the PD response ==========\n")

# -----------------------------------------------------------------------------
# 1. Inputs: Merton parameters, data, kernel, Bayesian satellite
# -----------------------------------------------------------------------------

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
delinq[, quarter := make_quarter(Date)]

z_out <- f_Z_estimation(delinq$value / 100)
P_TTC <- z_out$p_ttc
RHO   <- z_out$rho

DT_raw <- safe_fread(if (exists("PATH_DATA_VAR")) PATH_DATA_VAR else "data/processed/data_var_for_model.csv")
DT_raw[, Date := as.Date(Date)]

kern <- get_or_estimate_kernel(
  spec_name = SD_INFO_SET,
  vars_extra = VAR_SPECS[[SD_INFO_SET]],
  DT_raw = DT_raw,
  p_lags = p_lags, nrep = nrep, H = H, seed = seed, M_target = M_target,
  kernel_dir = if (exists("DIR_VAR_KERNELS")) DIR_VAR_KERNELS else file.path("output", "shared", "var_kernels"),
  impulse_name = impulse_name
)
variables  <- kern$variables
impulse_ix <- match(impulse_name, variables)
if (is.na(impulse_ix)) stop("Impulse variable not found in VAR kernel.", call. = FALSE)
M <- dim(kern$Psi_draws)[4]

bsat_path <- file.path("output", "applications", "dralacbn",
                       "main_specification", "satellite", "bayesian_satellite.rds")
if (!file.exists(bsat_path)) stop("Missing ", bsat_path, call. = FALSE)
bsat <- readRDS(bsat_path)
terms_df <- bsat$terms
n_sat <- nrow(bsat$beta_term_draws)
Lmax <- max(terms_df$lag)

# Dates of the kernel rows (complete cases, in order).
cc_mask <- stats::complete.cases(as.data.frame(DT_raw[, ..variables]))
dates_kept <- as.Date(DT_raw$Date)[cc_mask]
if (length(dates_kept) != nrow(kern$DT)) {
  stop("Date alignment failed between DT_raw and the kernel.", call. = FALSE)
}
quarters_kept <- make_quarter(dates_kept)
Y_all <- as.matrix(kern$DT[, ..variables])
Tn <- nrow(Y_all)

# End-of-sample state: the shock arrives in the quarter after the last
# observation of the kernel sample.
CURRENT_Q <- make_quarter(dates_kept[Tn] + 92)
CURRENT_LABEL <- paste0("End of sample (",
                        sub("Q", ":Q", CURRENT_Q), ")")
cat("Kernel sample ends ", quarters_kept[Tn],
    "; end-of-sample shock quarter: ", CURRENT_Q, "\n", sep = "")

# History end row for a shock occurring in quarter q: the row of q-1, i.e.
# the last row strictly BEFORE q. The end-of-sample state uses row Tn.
history_end_row <- function(shock_quarter) {
  idx <- match(shock_quarter, quarters_kept)
  if (is.na(idx)) stop("Quarter not in sample: ", shock_quarter, call. = FALSE)
  if (idx - 1L < max(p_lags, Lmax + 1L)) {
    stop("Not enough history before ", shock_quarter, call. = FALSE)
  }
  idx - 1L
}

# -----------------------------------------------------------------------------
# 2. History-invariant objects (computed once on the full posterior)
# -----------------------------------------------------------------------------

cat("Computing conditional variances and psi_Z (history-invariant) ...\n")
mom_end <- compute_factor_moment_draws_bayes(kern, bsat, impulse_ix, p_lags, H)
psiZ_1sd <- inject_girf_into_Z_bayes(kern$Psi_draws, kern$Sigma_kept, bsat,
                                     impulse_ix, shock_scale = 1)

# Episode shock sizes (standardized innovations at the episode quarters).
episode_sizes <- vapply(EPISODES, function(ep) {
  compute_structural_e1_median(
    kernel = kern, dates_raw = dates_kept,
    target_quarter = ep$quarter, impulse_idx = impulse_ix
  )$shock_scale
}, numeric(1))

# -----------------------------------------------------------------------------
# 3. Baseline mean path of Z conditional on an arbitrary history (per draw)
# -----------------------------------------------------------------------------

compute_mu_draws_at <- function(end_row) {
  Y_hist <- Y_all[(end_row - p_lags + 1L):end_row, , drop = FALSE]
  hist_tail <- Y_all[(end_row - Lmax):end_row, , drop = FALSE]
  mu_draws <- matrix(NA_real_, nrow = H + 1L, ncol = M)
  for (m in seq_len(M)) {
    sm <- pair_satellite_draw(m, n_sat)
    beta_m <- as.numeric(bsat$beta_term_draws[sm, terms_df$term, drop = TRUE])
    beta0_m <- bsat$beta0_draws[sm]
    Y_fore <- forecast_baseline_path(kern$B_kept[, , m], Y_hist, p_lags, H)
    mu_h <- rep(beta0_m, H + 1L)
    for (r in seq_len(nrow(terms_df))) {
      base_j <- terms_df$col[r]
      lag_r <- terms_df$lag[r]
      for (hh in 0:H) {
        val <- if (hh >= lag_r) {
          Y_fore[hh - lag_r + 1L, base_j]
        } else {
          hist_tail[(Lmax + 1L) - (lag_r - hh), base_j]
        }
        mu_h[hh + 1L] <- mu_h[hh + 1L] + beta_m[r] * val
      }
    }
    mu_draws[, m] <- mu_h
  }
  mu_draws
}

# -----------------------------------------------------------------------------
# 4. Named states: full-posterior responses
# -----------------------------------------------------------------------------

states <- c(
  lapply(seq_along(EPISODES), function(i) {
    list(quarter = EPISODES[[i]]$quarter, label = EPISODES[[i]]$label,
         end_row = history_end_row(EPISODES[[i]]$quarter),
         own_size = episode_sizes[i])
  }),
  list(list(quarter = CURRENT_Q, label = CURRENT_LABEL,
            end_row = Tn, own_size = NA_real_))
)

bands_from_draws <- function(draw_mat, extra = list()) {
  qs <- c(0.05, 0.16, 0.50, 0.84, 0.95)
  qmat <- t(apply(draw_mat, 1, stats::quantile, probs = qs, na.rm = TRUE))
  out <- data.table(
    horizon = 0:(nrow(draw_mat) - 1L),
    lower90 = qmat[, 1], lower68 = qmat[, 2], median = qmat[, 3],
    upper68 = qmat[, 4], upper90 = qmat[, 5]
  )
  for (nm in names(extra)) out[, (nm) := extra[[nm]]]
  out[]
}

state_rows <- list()
state_bands <- list()
state_draws_own <- list()   # own-size PD draws, kept for paired comparisons

for (st in states) {
  cat("State: ", st$label, " ...\n", sep = "")
  mu_st <- if (st$end_row == Tn) mom_end$mu_draws else compute_mu_draws_at(st$end_row)

  # Baseline expected PD at the shock quarter (h = 0), in pp.
  qpi <- qnorm(P_TTC)
  pd0 <- 100 * pnorm((qpi - sqrt(RHO) * mu_st[1L, ]) /
                       sqrt((1 - RHO) + RHO * pmax(mom_end$s2_draws[1L, ], 0)))

  one_state <- function(scale, shock_lab) {
    pd <- compute_pd_girf(psiZ_1sd$draws * scale, mu_st,
                          mom_end$s2_draws, mom_end$s2_delta_draws,
                          p = P_TTC, rho = RHO)
    dr <- 100 * pd$draws
    b <- bands_from_draws(dr, list(state = st$label, shock = shock_lab))
    pk_h <- b[which.max(abs(median)), horizon]
    list(bands = b, peak = b[horizon == pk_h],
         prob_pos = mean(dr[pk_h + 1L, ] > 0, na.rm = TRUE),
         draws = dr)
  }

  r1 <- one_state(1, "one_sd")
  r2 <- if (is.finite(st$own_size)) one_state(st$own_size, "own_size") else NULL
  if (!is.null(r2)) state_draws_own[[st$label]] <- r2$draws

  d_obs <- delinq[quarter == st$quarter, value]
  state_rows[[st$label]] <- data.table(
    state = st$label, quarter = st$quarter,
    delinquency_obs = if (length(d_obs)) d_obs else NA_real_,
    pd0_median_pp = stats::median(pd0),
    own_shock_sd = st$own_size,
    peak_1sd_pp = r1$peak$median, peak_1sd_h = r1$peak$horizon,
    peak_1sd_lo68 = r1$peak$lower68, peak_1sd_hi68 = r1$peak$upper68,
    prob_pos_1sd = r1$prob_pos,
    peak_own_pp = if (!is.null(r2)) r2$peak$median else NA_real_,
    peak_own_lo68 = if (!is.null(r2)) r2$peak$lower68 else NA_real_,
    peak_own_hi68 = if (!is.null(r2)) r2$peak$upper68 else NA_real_
  )
  state_bands[[st$label]] <- rbind(r1$bands, if (!is.null(r2)) r2$bands)
}

summary_dt <- rbindlist(state_rows)
ref_peak <- summary_dt[state == CURRENT_LABEL, peak_1sd_pp]
summary_dt[, amplification_vs_current := peak_1sd_pp / ref_peak]
fwrite(summary_dt, file.path(sd_out_dir, "tables", "state_dependent_responses.csv"))

# Paired posterior comparison of the two largest episodes under their own
# shock sizes: the draws are paired (identical posterior draws across states),
# so P(Gulf > September 11) is the exact Bayesian comparison and is much
# sharper than the overlap of the marginal bands suggests.
g_dr <- state_draws_own[["Gulf War (1990:Q3)"]]
s_dr <- state_draws_own[["September 11 (2001:Q3)"]]
prob_gulf_gt_911 <- if (!is.null(g_dr) && !is.null(s_dr)) {
  data.table(
    horizon = 0:H,
    prob_gulf_greater = vapply(seq_len(H + 1L), function(i) {
      mean(g_dr[i, ] > s_dr[i, ], na.rm = TRUE)
    }, numeric(1))
  )
} else NULL
if (!is.null(prob_gulf_gt_911)) {
  fwrite(prob_gulf_gt_911,
         file.path(sd_out_dir, "tables", "prob_gulf_exceeds_september11.csv"))
}

bands_all <- rbindlist(state_bands)
fwrite(bands_all, file.path(sd_out_dir, "irf", "pd_girf_bands_by_state_pp.csv"))

# -----------------------------------------------------------------------------
# 5. Convexity curve at posterior-median parameters, all admissible histories
# -----------------------------------------------------------------------------

cat("Computing the convexity curve over all histories (posterior-median parameters) ...\n")

B_bar <- apply(kern$B_kept, c(1, 2), median)
Sigma_bar <- apply(kern$Sigma_kept, c(1, 2), median)
Sigma_bar <- 0.5 * (Sigma_bar + t(Sigma_bar))
Psi_bar <- compute_ma_coefficients(B_bar, p = p_lags, H = H)
dimnames(Psi_bar) <- list(horizon = 0:H, variable = variables, shock_col = variables)

beta0_bar <- unname(bsat$beta_ols["(Intercept)"])
beta_bar <- as.numeric(bsat$beta_ols[terms_df$term])
sigma2_bar <- bsat$sigma2_ols

# psi_Z and conditional variances at the median parameters (history-invariant).
sgg <- Sigma_bar[impulse_ix, impulse_ix]
delta_bar <- Sigma_bar[, impulse_ix] / sqrt(sgg)
girfY_bar <- matrix(0.0, H + 1L, length(variables))
for (hh in 0:H) girfY_bar[hh + 1L, ] <- Psi_bar[hh + 1L, , ] %*% delta_bar
psiZ_bar <- numeric(H + 1L)
for (r in seq_len(nrow(terms_df))) {
  lag_r <- terms_df$lag[r]; col_r <- terms_df$col[r]
  contrib <- if (lag_r == 0) girfY_bar[, col_r] else
    c(rep(0, lag_r), girfY_bar[1:(H + 1L - lag_r), col_r])
  psiZ_bar <- psiZ_bar + beta_bar[r] * contrib
}
s2_bar <- compute_s2_one_draw(Psi_bar, Sigma_bar, terms_df, beta_bar,
                              sigma2_bar, impulse_idx = impulse_ix)

pd_curve_one <- function(end_row) {
  Y_hist <- Y_all[(end_row - p_lags + 1L):end_row, , drop = FALSE]
  hist_tail <- Y_all[(end_row - Lmax):end_row, , drop = FALSE]
  Y_fore <- forecast_baseline_path(B_bar, Y_hist, p_lags, H)
  mu_h <- rep(beta0_bar, H + 1L)
  for (r in seq_len(nrow(terms_df))) {
    base_j <- terms_df$col[r]; lag_r <- terms_df$lag[r]
    for (hh in 0:H) {
      val <- if (hh >= lag_r) Y_fore[hh - lag_r + 1L, base_j] else
        hist_tail[(Lmax + 1L) - (lag_r - hh), base_j]
      mu_h[hh + 1L] <- mu_h[hh + 1L] + beta_bar[r] * val
    }
  }
  qpi <- qnorm(P_TTC)
  s_b <- sqrt((1 - RHO) + RHO * pmax(s2_bar$s2, 0))
  s_d <- sqrt((1 - RHO) + RHO * pmax(s2_bar$s2_delta, 0))
  pd_b <- pnorm((qpi - sqrt(RHO) * mu_h) / s_b)
  pd_s <- pnorm((qpi - sqrt(RHO) * (mu_h + psiZ_bar)) / s_d)
  dpd <- 100 * (pd_s - pd_b)
  i <- which.max(abs(dpd))
  c(pd0 = 100 * pd_b[1L], peak = dpd[i], peak_h = i - 1L)
}

first_row <- max(p_lags, Lmax + 1L)
curve_rows <- first_row:Tn
curve <- t(vapply(curve_rows, pd_curve_one, numeric(3)))
curve_dt <- data.table(
  end_row = curve_rows,
  shock_quarter = c(quarters_kept[curve_rows[-length(curve_rows)] + 1L], CURRENT_Q),
  pd0_pp = curve[, "pd0"], peak_dpd_pp = curve[, "peak"], peak_h = curve[, "peak_h"]
)
fwrite(curve_dt, file.path(sd_out_dir, "tables", "convexity_curve_all_histories.csv"))

# -----------------------------------------------------------------------------
# 6. Figures
# -----------------------------------------------------------------------------

state_levels <- vapply(states, function(s) s$label, character(1))
COLS <- stats::setNames(c("#6F1732", "#AB4A7D", "#D78FB4", "#6B6B6B"), state_levels)

# (a) PD GIRF by state, common one-s.d. shock.
#     Shown states: Gulf War, September 11, end of sample (Ukraine omitted:
#     its benign state nearly coincides with the end-of-sample one).
sel_states <- c("Gulf War (1990:Q3)", "September 11 (2001:Q3)", CURRENT_LABEL)
b1 <- bands_all[shock == "one_sd" & state %in% sel_states]
b1[, state := factor(state, levels = state_levels)]
p_states <- ggplot(b1, aes(x = horizon, y = median, color = state, fill = state)) +
  geom_ribbon(aes(ymin = lower90, ymax = upper90), alpha = 0.10, color = NA) +
  geom_hline(yintercept = 0, color = "grey45", linewidth = 0.35, linetype = "dashed") +
  geom_line(linewidth = 1.0) +
  geom_point(size = 1.4) +
  scale_color_manual(values = COLS) +
  scale_fill_manual(values = COLS) +
  scale_x_continuous(breaks = seq(0, 12, by = 1), minor_breaks = NULL) +
  labs(
    x = "Horizon (quarters)",
    y = expression(Delta ~ "PD (percentage points)"),
    caption = "Posterior-median responses to a common one-s.d. GPR innovation, conditional on the history preceding each episode, with 90% pointwise credible bands. The macro propagation is identical across states; only the position in the Merton-Vasicek map differs."
  ) +
  theme_paper()
save_fig(p_states, "fig_state_dependence_girf", width = 7.4, height = 4.6)

# (b) Convexity curve with named episodes highlighted.
ggrepel_or_text <- function(dt) {
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    ggrepel::geom_text_repel(data = dt, aes(label = label, color = label),
                             size = 3.1, show.legend = FALSE, seed = 1)
  } else {
    geom_text(data = dt, aes(label = label, color = label),
              size = 3.1, vjust = -0.9, show.legend = FALSE)
  }
}

ep_quarters <- c(vapply(EPISODES, function(e) e$quarter, character(1)), CURRENT_Q)
curve_dt[, named := shock_quarter %in% ep_quarters]
named_dt <- merge(
  curve_dt[named == TRUE],
  data.table(shock_quarter = ep_quarters, label = state_levels),
  by = "shock_quarter"
)

p_curve <- ggplot(curve_dt, aes(x = pd0_pp, y = peak_dpd_pp)) +
  geom_point(color = "grey60", size = 1.1, alpha = 0.55) +
  geom_point(data = named_dt, aes(color = label), size = 2.6) +
  ggrepel_or_text(named_dt) +
  scale_color_manual(values = COLS) +
  labs(
    x = "Baseline expected PD at the shock date (percentage points)",
    y = expression("Peak " * Delta * "PD (percentage points)"),
    caption = "Each grey point is one admissible history in the sample; the peak response to a one-s.d. GPR innovation is computed at posterior-median parameters. Colored points: histories preceding the three largest positive GPR innovations and the end-of-sample state."
  ) +
  guides(color = "none") +
  theme_paper()
save_fig(p_curve, "fig_state_dependence_curve", width = 7.4, height = 4.8)

# (c) The size-state inversion: the two largest episodes under their OWN
#     historical shock sizes. The Gulf War shock is smaller (2.46 s.d. vs
#     3.84) but strikes a weaker credit state, and produces a uniformly
#     larger PD response.
b_own <- bands_all[shock == "own_size" &
                     state %in% c("Gulf War (1990:Q3)",
                                  "September 11 (2001:Q3)")]
b_own[, state := factor(state, levels = state_levels)]

p_inv <- ggplot(b_own, aes(x = horizon, y = median,
                           color = state, fill = state)) +
  geom_ribbon(aes(ymin = lower90, ymax = upper90), alpha = 0.10, color = NA) +
  geom_hline(yintercept = 0, color = "grey45", linewidth = 0.35, linetype = "dashed") +
  geom_line(linewidth = 1.0) +
  geom_point(size = 1.4) +
  scale_color_manual(values = COLS) +
  scale_fill_manual(values = COLS) +
  scale_x_continuous(breaks = seq(0, 12, by = 1), minor_breaks = NULL) +
  labs(
    x = "Horizon (quarters)",
    y = expression(Delta ~ "PD (percentage points)"),
    caption = "Each episode is replayed at its own date with its own historical shock size (Gulf War: 2.46 s.d.; September 11: 3.84 s.d.). Posterior medians with 90% pointwise credible bands. The smaller shock produces the larger response because it strikes a weaker credit state."
  ) +
  theme_paper()
save_fig(p_inv, "fig_state_dependence_inversion", width = 7.4, height = 4.6)

# -----------------------------------------------------------------------------
# 7. Save, recap
# -----------------------------------------------------------------------------

saveRDS(
  list(info_set = SD_INFO_SET, p_ttc = P_TTC, rho = RHO, n_draws = M,
       episodes = EPISODES, episode_sizes = episode_sizes,
       summary = summary_dt, bands_pp = bands_all, curve = curve_dt,
       prob_gulf_gt_911 = prob_gulf_gt_911),
  file.path(sd_out_dir, "state_dependence_results.rds")
)

readme <- c(
  "# DRALACBN -- state dependence of the PD response",
  "",
  "The PD GIRF (Proposition 3) is history-dependent through the conditional",
  "mean of the systematic factor. This exercise evaluates the same shock at",
  "the histories preceding the three largest positive GPR innovations in the",
  "sample (Gulf War 1990Q3, September 11 2001Q3, Ukraine 2022Q1) and at the",
  "end-of-sample state, on the full posterior; and computes the peak response",
  "for EVERY admissible history at posterior-median parameters (convexity",
  "curve). Parameters are full-sample estimates: this is a counterfactual",
  "within the estimated model, not a real-time analysis.",
  "",
  "Generated by code/dralacbn/10_dralacbn_state_dependence.R (standalone).",
  "",
  "Contents:",
  "  tables/state_dependent_responses.csv     per state: peaks, bands, P(>0), amplification",
  "  tables/convexity_curve_all_histories.csv peak response vs initial PD, all histories",
  "  figures/fig_state_dependence_girf        median responses by state, common shock",
  "  figures/fig_state_dependence_curve       convexity curve with named episodes",
  "  irf/pd_girf_bands_by_state_pp.csv        full bands by state and shock",
  "  state_dependence_results.rds             all objects"
)
writeLines(readme, file.path(sd_out_dir, "README.md"))

cat("\n--- State dependence (one-s.d. shock, posterior medians) ----------------\n")
print(summary_dt[, .(state, pd0_pp = round(pd0_median_pp, 2),
                     peak_pp = round(peak_1sd_pp, 4), h = peak_1sd_h,
                     P_pos = round(prob_pos_1sd, 2),
                     amplif_vs_current = round(amplification_vs_current, 2))])
if (!is.null(prob_gulf_gt_911)) {
  cat(sprintf("Paired comparison at the peak (h=3): P(Gulf > September 11) = %.3f\n",
              prob_gulf_gt_911[horizon == 3, prob_gulf_greater]))
}
cat("Outputs written to: ", sd_out_dir, "\n", sep = "")
cat("--------------------------------------------------------------------------\n")

invisible(TRUE)
