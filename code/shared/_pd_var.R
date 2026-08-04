# =============================================================================
#  code/shared/_pd_var.R — Value-at-Risk (stressed PD) response
#
#  Two distinct quantiles live in this framework and must not be conflated:
#    - epistemic quantiles: the posterior credible bands of a response, i.e.
#      parameter/estimation uncertainty (produced by compute_pd_girf);
#    - aleatoric quantiles: the tail of the systematic factor distribution,
#      i.e. portfolio risk, which is what a Value-at-Risk measures.
#
#  Under the asymptotic single-risk-factor (Vasicek) assumption, the portfolio
#  default rate conditional on the systematic state Z is pi(Z). Since pi is
#  decreasing in Z, the alpha-level loss VaR evaluates pi at the lower tail of
#  Z. With Z_{t+h} ~ N(mu, s^2) (s^2 already contains the VAR forecast variance
#  and the satellite residual variance), the stressed (downturn) PD is
#
#     PD^VaR_alpha(h) = Phi( ( Phi^{-1}(p) - sqrt(rho) * mu
#                              + sqrt(rho) * Phi^{-1}(alpha) * s ) / sqrt(1-rho) ).
#
#  This is a different functional of N(mu, s^2) than the mean response
#  E[pi(Z)] = Phi( (Phi^{-1}(p) - sqrt(rho) mu) / sqrt(1-rho+rho s^2) ):
#  the VaR shifts the factor quantile in the numerator, the mean integrates the
#  variance in the denominator.
#
#  compute_pd_var_girf() returns the GIRF of the stressed PD (shock minus
#  baseline), evaluated per posterior draw, so the reported bands are still the
#  epistemic credible bands AROUND the (aleatoric) VaR. The two uncertainty
#  layers are thus kept separate: alpha indexes risk, the bands index estimation
#  uncertainty.
# =============================================================================

# Stressed (downturn) PD at confidence alpha, given factor moments (mu, s^2).
.pd_var_level <- function(mu, s, p, rho, za) {
  stats::pnorm((stats::qnorm(p) - sqrt(rho) * mu + sqrt(rho) * za * s) /
                 sqrt(1 - rho))
}

compute_pd_var_girf <- function(psiZ_draws, mu_draws, s2_draws, s2_delta_draws,
                                p, rho, alpha) {
  stopifnot(alpha > 0, alpha < 1, rho > 0, rho < 1, p > 0, p < 1)
  H <- nrow(psiZ_draws) - 1
  M <- ncol(psiZ_draws)
  za <- stats::qnorm(alpha)

  out <- matrix(NA_real_, nrow = H + 1, ncol = M)

  for (m in seq_len(M)) {
    mu_base <- mu_draws[, m]
    mu_shock <- mu_base + psiZ_draws[, m]
    s_base <- sqrt(pmax(s2_draws[, m], 0))
    s_shock <- sqrt(pmax(s2_delta_draws[, m], 0))

    pd_var_base <- .pd_var_level(mu_base, s_base, p, rho, za)
    pd_var_shock <- .pd_var_level(mu_shock, s_shock, p, rho, za)

    out[, m] <- pd_var_shock - pd_var_base
  }

  qs <- c(0.05, 0.16, 0.50, 0.84, 0.95)
  qmat <- t(apply(out, 1, stats::quantile, probs = qs, na.rm = TRUE))

  bands <- data.table(
    horizon = 0:H,
    lower90 = qmat[, 1],
    lower68 = qmat[, 2],
    median = qmat[, 3],
    upper68 = qmat[, 4],
    upper90 = qmat[, 5],
    p = p,
    rho = rho,
    alpha = alpha
  )

  list(draws = out, bands = bands)
}

# =============================================================================
#  Expected Shortfall (coherent tail PD) response
#
#  The alpha-level expected shortfall of the portfolio default probability is the
#  average PD over its worst 1-alpha tail. Since pi is decreasing in Z, this is
#
#     ES_alpha(h) = E[ pi(Z_{t+h}) | Z_{t+h} <= z_alpha ],
#                   z_alpha = mu + s * qnorm(1 - alpha),   P(Z <= z_alpha) = 1-alpha.
#
#  Closed form (Corollary "Expected-shortfall response", see the paper):
#
#     ES_alpha(h) = Phi2( m, qnorm(1-alpha) ; rho_star ) / (1 - alpha),
#        m        = (qnorm(p) - sqrt(rho) mu) / sqrt((1-rho) + rho s^2)   (= probit of mean PD),
#        rho_star = sqrt(rho) s / sqrt((1-rho) + rho s^2),
#     with Phi2 the standard bivariate normal CDF of correlation rho_star.
#
#  We evaluate the equivalent one-dimensional truncated integral
#     ES_alpha = (1/(1-alpha)) * Int_{-inf}^{w_alpha} pi(mu + s w) phi(w) dw,
#     w_alpha  = qnorm(1 - alpha),
#  by Simpson quadrature. The truncation limit w_alpha is the SAME across
#  horizons/draws, so the routine vectorises over (mu, s) and is fast; it equals
#  the Phi2 closed form to quadrature precision. (A Phi2 cross-check at a single
#  (mu, s) via mvtnorm::pmvnorm is provided in tests/.)
# =============================================================================

# Expected shortfall of the PD at confidence alpha, given factor moments
# (mu, s) supplied as vectors of equal length; returns a vector of ES values.
.pd_es_level <- function(mu, s, p, rho, alpha, n_grid = 513L) {
  if (n_grid %% 2L == 0L) n_grid <- n_grid + 1L  # Simpson needs an even #intervals
  qp  <- stats::qnorm(p)
  wa  <- stats::qnorm(1 - alpha)            # common upper limit
  lo  <- min(wa - 1e-8, -8)                 # lower limit (far left tail)
  w   <- seq(lo, wa, length.out = n_grid)   # Simpson nodes
  dw  <- (wa - lo) / (n_grid - 1)
  sw  <- c(1, rep(c(4, 2), length.out = n_grid - 2), 1)  # Simpson weights
  phw <- stats::dnorm(w)
  # z_{i,j} = mu_i + s_i * w_j ; pi(z) = Phi((qp - sqrt(rho) z)/sqrt(1-rho))
  Z   <- outer(mu, w, function(a, b) a) + outer(s, w, function(a, b) a * b)
  PI  <- stats::pnorm((qp - sqrt(rho) * Z) / sqrt(1 - rho))
  integ <- PI %*% (sw * phw) * (dw / 3)     # length(mu) integral values
  as.numeric(integ) / (1 - alpha)
}

compute_pd_es_girf <- function(psiZ_draws, mu_draws, s2_draws, s2_delta_draws,
                               p, rho, alpha) {
  stopifnot(alpha > 0, alpha < 1, rho > 0, rho < 1, p > 0, p < 1)
  H <- nrow(psiZ_draws) - 1
  M <- ncol(psiZ_draws)

  out <- matrix(NA_real_, nrow = H + 1, ncol = M)

  for (m in seq_len(M)) {
    mu_base  <- mu_draws[, m]
    mu_shock <- mu_base + psiZ_draws[, m]
    s_base   <- sqrt(pmax(s2_draws[, m], 0))
    s_shock  <- sqrt(pmax(s2_delta_draws[, m], 0))

    es_base  <- .pd_es_level(mu_base,  s_base,  p, rho, alpha)
    es_shock <- .pd_es_level(mu_shock, s_shock, p, rho, alpha)

    out[, m] <- es_shock - es_base
  }

  qs <- c(0.05, 0.16, 0.50, 0.84, 0.95)
  qmat <- t(apply(out, 1, stats::quantile, probs = qs, na.rm = TRUE))

  bands <- data.table(
    horizon = 0:H,
    lower90 = qmat[, 1],
    lower68 = qmat[, 2],
    median = qmat[, 3],
    upper68 = qmat[, 4],
    upper90 = qmat[, 5],
    p = p,
    rho = rho,
    alpha = alpha
  )

  list(draws = out, bands = bands)
}

# -----------------------------------------------------------------------------
# Overlay of the mean PD response and the stressed (VaR) PD responses, with the
# 68% credible band of the mean shaded. var_bands is a named list
# (e.g. "VaR 99%" -> bands, "VaR 99.9%" -> bands).
# -----------------------------------------------------------------------------
plot_pd_mean_vs_var <- function(mean_bands, var_bands,
                                ylab_txt = expression(Delta ~ PD)) {
  mb <- as.data.table(copy(mean_bands))

  parts <- c(
    list(Mean = mb[, .(horizon, median)]),
    lapply(var_bands, function(b) as.data.table(b)[, .(horizon, median)])
  )
  long <- rbindlist(
    lapply(names(parts), function(nm) {
      d <- copy(parts[[nm]]); d[, series := nm]; d
    }),
    use.names = TRUE
  )
  long[, series := factor(series, levels = c("Mean", names(var_bands)))]

  # The mean response is drawn in the same colour as the stand-alone PD
  # response figure (#6F1732); the tail (VaR or ES) overlays take distinct
  # colours and line types. Colours are assigned by position so the scheme
  # works for any number of tail series.
  lvls    <- levels(long$series)
  n_tail  <- length(lvls) - 1L
  col_tail <- c("#1F6FB2", "#E08214", "#4DAF4A", "#984EA3")
  lty_tail <- c("longdash", "dotted", "dotdash", "twodash")
  col_vals <- setNames(c("#6F1732", col_tail[seq_len(n_tail)]), lvls)
  lty_vals <- setNames(c("solid",   lty_tail[seq_len(n_tail)]), lvls)

  ggplot() +
    geom_ribbon(
      data = mb, aes(x = horizon, ymin = lower68, ymax = upper68),
      fill = "#6F1732", alpha = 0.18
    ) +
    geom_hline(yintercept = 0, color = "grey45", linewidth = 0.35,
               linetype = "dashed") +
    geom_line(
      data = long,
      aes(x = horizon, y = median, color = series, linetype = series),
      linewidth = 1.0
    ) +
    scale_color_manual(values = col_vals) +
    scale_linetype_manual(values = lty_vals) +
    scale_x_continuous(
      breaks = function(x) seq(max(0, ceiling(x[1])), floor(x[2]), by = 1),
      minor_breaks = NULL
    ) +
    labs(x = "Horizon (quarters)", y = ylab_txt, color = NULL, linetype = NULL) +
    theme_robustness() +
    theme(legend.position = "top")
}

# -----------------------------------------------------------------------------
# Engine driver: compute and write the VaR (stressed-PD) responses for both
# shocks and all alpha levels, with per-shock overlays. Pure side-effect.
# -----------------------------------------------------------------------------
run_pd_var_robustness <- function(out_dir, psiZ_1sd, psiZ_2001,
                                  mean_bands_1sd, mean_bands_2001, moments,
                                  app_name, info_set_name, info_set_label,
                                  p, rho, alphas) {
  tryCatch({
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    var_shocks <- list(
      list(lab = "one_sd", psi = psiZ_1sd, mean_bands = mean_bands_1sd),
      list(lab = "2001Q3", psi = psiZ_2001, mean_bands = mean_bands_2001)
    )

    var_bands_all <- list()
    var_peaks_all <- list()
    es_bands_all  <- list()
    es_peaks_all  <- list()

    for (sh in var_shocks) {
      overlay_var <- list()
      overlay_es  <- list()
      for (a in alphas) {
        a_tag <- gsub("\\.", "p", format(a, trim = TRUE))
        vr <- compute_pd_var_girf(
          psiZ_draws = sh$psi$draws,
          mu_draws = moments$mu_draws,
          s2_draws = moments$s2_draws,
          s2_delta_draws = moments$s2_delta_draws,
          p = p, rho = rho, alpha = a
        )
        # Coherent tail measure: expected-shortfall response (same moments).
        er <- compute_pd_es_girf(
          psiZ_draws = sh$psi$draws,
          mu_draws = moments$mu_draws,
          s2_draws = moments$s2_draws,
          s2_delta_draws = moments$s2_delta_draws,
          p = p, rho = rho, alpha = a
        )

        b <- as.data.table(copy(vr$bands))
        b[, `:=`(application = app_name, information_set = info_set_name,
                 info_set_label = info_set_label, shock = sh$lab)]
        var_bands_all[[paste(a_tag, sh$lab)]] <- b

        be <- as.data.table(copy(er$bands))
        be[, `:=`(application = app_name, information_set = info_set_name,
                  info_set_label = info_set_label, shock = sh$lab)]
        es_bands_all[[paste(a_tag, sh$lab)]] <- be

        pk <- vr$bands[which.max(abs(vr$bands$median))]
        var_peaks_all[[paste(a_tag, sh$lab)]] <- data.table(
          application = app_name, information_set = info_set_name,
          shock = sh$lab, alpha = a,
          peak_median = pk$median, peak_horizon = pk$horizon,
          lower68 = pk$lower68, upper68 = pk$upper68
        )

        pke <- er$bands[which.max(abs(er$bands$median))]
        es_peaks_all[[paste(a_tag, sh$lab)]] <- data.table(
          application = app_name, information_set = info_set_name,
          shock = sh$lab, alpha = a,
          peak_median = pke$median, peak_horizon = pke$horizon,
          lower68 = pke$lower68, upper68 = pke$upper68
        )

        save_png_plot(
          plot_single_bands(vr$bands, ylab_txt = expression(Delta ~ PD ~ VaR)),
          out_dir, sprintf("pd_var_girf_%s_alpha%s", sh$lab, a_tag),
          width = 6.5, height = 4.2
        )
        save_png_plot(
          plot_single_bands(er$bands, ylab_txt = expression(Delta ~ PD ~ ES)),
          out_dir, sprintf("pd_es_girf_%s_alpha%s", sh$lab, a_tag),
          width = 6.5, height = 4.2
        )
        overlay_var[[sprintf("VaR %g%%", 100 * a)]] <- vr$bands
        overlay_es[[sprintf("ES %g%%", 100 * a)]]   <- er$bands
      }

      save_png_plot(
        plot_pd_mean_vs_var(sh$mean_bands, overlay_var,
                            ylab_txt = expression(Delta ~ PD)),
        out_dir, sprintf("pd_mean_vs_var_%s", sh$lab),
        width = 6.5, height = 4.2
      )
      save_png_plot(
        plot_pd_mean_vs_var(sh$mean_bands, overlay_es,
                            ylab_txt = expression(Delta ~ PD)),
        out_dir, sprintf("pd_mean_vs_es_%s", sh$lab),
        width = 6.5, height = 4.2
      )
    }

    fwrite(rbindlist(var_bands_all, use.names = TRUE, fill = TRUE),
           file.path(out_dir, "pd_var_response_bands.csv"))
    fwrite(rbindlist(var_peaks_all, use.names = TRUE, fill = TRUE),
           file.path(out_dir, "pd_var_peaks.csv"))
    fwrite(rbindlist(es_bands_all, use.names = TRUE, fill = TRUE),
           file.path(out_dir, "pd_es_response_bands.csv"))
    fwrite(rbindlist(es_peaks_all, use.names = TRUE, fill = TRUE),
           file.path(out_dir, "pd_es_peaks.csv"))
  }, error = function(e) {
    cat("    PD-VaR skipped for ", info_set_name, ": ",
        conditionMessage(e), "\n", sep = "")
  })
  invisible(TRUE)
}
