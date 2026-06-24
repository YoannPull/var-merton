# =============================================================================
#  code/_bvar_niw_utils.R — Shared VAR / BVAR engine (single source of truth)
#
#  Centralizes the reusable VAR/BVAR building blocks that used to be copied,
#  with small inconsistencies, across:
#     - code/04_girf.R
#     - code/05_stress_scenario.R
#     - code/shared/_helpers_var_merton.R
#
#  This file is sourced ONCE, from code/00_setup.R, so every analysis script
#  that begins with source("code/00_setup.R") gets the same definitions.
#  Scripts must NOT redefine these functions locally.
#
#  It defines, in this order:
#     simulate_bvar_niw()        Normal-Inverse-Wishart Bayesian VAR sampler
#     build_companion()          companion matrix of a reduced-form VAR
#     max_root()                 largest companion eigenvalue modulus
#     is_stable()                stability test (|lambda|_max < 1)
#     extract_c_A_list()         split B into intercept c and lag matrices A_l
#     compute_ma_coefficients()  MA / impulse coefficients Psi(h)
#     forecast_baseline_path()   deterministic baseline forecast path
#
#  Methodological note (identical to the main pipeline and the robustness
#  engine): GIRFs are computed explicitly as Pesaran-Shin / Koop-Pesaran-Potter
#  generalized impulse responses,
#        delta_j = Sigma[, j] / sqrt(Sigma[j, j]),
#  which, for a GPR shock with GPR ordered first, equals the first Cholesky
#  impulse response of a Gaussian VAR without relying on a Cholesky ordering.
#
#  Dependency: the mniw package (mniw::rMNIW) for the posterior sampler.
# =============================================================================

# -----------------------------------------------------------------------------
# Bayesian VAR (Normal-Inverse-Wishart) posterior sampler
#
# Robust version: symmetrizes and ridge-regularizes the cross-product matrices,
# retries failed draws with a growing ridge, and returns rich metadata
# (m, T_eff, nu_post) used by the GIRF and stress scripts.
#
# Returns a list with:
#   B        array (m x k x nrep) of coefficient draws, intercept in row 1
#   S        array (k x k x nrep) of residual-covariance draws (Sigma)
#   p, k     lag order and number of variables
#   m        number of regressors per equation (1 + k * p)
#   T_eff    effective sample size after lag truncation
#   nu_post  posterior Inverse-Wishart degrees of freedom
# -----------------------------------------------------------------------------
simulate_bvar_niw <- function(Y, p = 2, nrep = 5000, seed = 123,
                              ridge = 1e-8, max_retry = 5) {
  set.seed(seed)

  Y <- as.matrix(Y)

  if (any(!is.finite(Y))) {
    stop("simulate_bvar_niw(): Y contains non-finite values.")
  }

  T <- nrow(Y)
  k <- ncol(Y)

  if (T <= p + 5) {
    stop("simulate_bvar_niw(): too few observations for VAR estimation.")
  }

  Y_t <- Y[(p + 1):T, , drop = FALSE]

  X_t <- do.call(
    cbind,
    lapply(1:p, function(l) {
      Y[(p + 1 - l):(T - l), , drop = FALSE]
    })
  )

  X_t <- cbind(1, X_t)

  T_eff <- nrow(X_t)
  m <- ncol(X_t)
  nu_post <- T_eff - m

  if (nu_post <= (k - 1)) {
    stop(sprintf(
      "Inverse-Wishart not defined: need T_eff - m > k - 1 (here %d <= %d).",
      nu_post, k - 1
    ))
  }

  make_spd <- function(A, ridge = 1e-8) {
    A <- as.matrix(A)
    A <- 0.5 * (A + t(A))

    eig <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
    min_eig <- min(eig, na.rm = TRUE)

    scale_A <- mean(diag(A), na.rm = TRUE)
    if (!is.finite(scale_A) || scale_A <= 0) {
      scale_A <- 1
    }

    add <- max(ridge * scale_A, -min_eig + ridge * scale_A, 0)
    A <- A + diag(add, nrow(A))
    A <- 0.5 * (A + t(A))

    A
  }

  XtX <- crossprod(X_t)
  XtY <- crossprod(X_t, Y_t)

  XtX_spd <- make_spd(XtX, ridge = ridge)

  B_ols <- tryCatch(
    solve(XtX_spd, XtY),
    error = function(e) {
      qr.coef(qr(XtX_spd), XtY)
    }
  )

  if (any(!is.finite(B_ols))) {
    stop("simulate_bvar_niw(): B_ols contains non-finite values.")
  }

  U_ols <- Y_t - X_t %*% B_ols
  S <- crossprod(U_ols)
  S <- make_spd(S, ridge = ridge)

  XtX_inv <- tryCatch(
    chol2inv(chol(XtX_spd)),
    error = function(e) {
      solve(XtX_spd)
    }
  )

  XtX_inv <- make_spd(XtX_inv, ridge = ridge)

  draws_B <- array(NA_real_, dim = c(m, k, nrep))
  draws_Sigma <- array(NA_real_, dim = c(k, k, nrep))

  for (i in seq_len(nrep)) {

    draw_ok <- FALSE
    local_ridge <- ridge

    for (rr in seq_len(max_retry)) {

      S_try <- make_spd(S, ridge = local_ridge)
      XtX_inv_try <- make_spd(XtX_inv, ridge = local_ridge)

      d <- tryCatch(
        mniw::rMNIW(
          1,
          Lambda = B_ols,
          Sigma  = XtX_inv_try,
          Psi    = S_try,
          nu     = nu_post
        ),
        error = function(e) NULL
      )

      if (!is.null(d) &&
          all(is.finite(d$X)) &&
          all(is.finite(d$V))) {

        draws_B[, , i] <- d$X
        draws_Sigma[, , i] <- 0.5 * (d$V + t(d$V))
        draw_ok <- TRUE
        break
      }

      local_ridge <- local_ridge * 10
    }

    if (!draw_ok) {
      stop(
        "simulate_bvar_niw(): mniw::rMNIW failed after ",
        max_retry,
        " retries at draw ",
        i,
        ". Try increasing ridge or checking collinearity in the VAR variables."
      )
    }
  }

  list(
    B = draws_B,
    S = draws_Sigma,
    p = p,
    k = k,
    m = m,
    T_eff = T_eff,
    nu_post = nu_post
  )
}

# -----------------------------------------------------------------------------
# Companion matrix of a reduced-form VAR(p)
#   B is (m x k) = (1 + k*p) x k, intercept in the first row.
# -----------------------------------------------------------------------------
build_companion <- function(B, k, p) {
  A_comp <- matrix(0, nrow = k * p, ncol = k * p)
  A_stack <- t(B[-1, , drop = FALSE])
  A_comp[1:k, 1:(k * p)] <- A_stack

  if (p > 1) {
    A_comp[(k + 1):(k * p), 1:(k * (p - 1))] <- diag(k * (p - 1))
  }

  A_comp
}

# -----------------------------------------------------------------------------
# Largest modulus of the companion eigenvalues (spectral radius).
# Argument order (B, p, k) matches the call sites in 04/05 and the robustness
# engine.
# -----------------------------------------------------------------------------
max_root <- function(B, p, k) {
  max(Mod(eigen(build_companion(B, k, p), only.values = TRUE)$values))
}

# -----------------------------------------------------------------------------
# Dynamic stability test: TRUE if the largest companion root is below 1
# (with a small numerical margin).
# -----------------------------------------------------------------------------
is_stable <- function(B, k, p) {
  rho <- max(Mod(eigen(build_companion(B, k, p), only.values = TRUE)$values))
  rho < 1 - 1e-10
}

# -----------------------------------------------------------------------------
# Split coefficient matrix B into the intercept vector c and the list of lag
# matrices A_1, ..., A_p (each k x k).
# -----------------------------------------------------------------------------
extract_c_A_list <- function(B, k, p) {
  c_vec <- as.numeric(B[1, ])
  A_stack <- t(B[-1, , drop = FALSE])

  A_list <- vector("list", p)

  for (i in seq_len(p)) {
    A_list[[i]] <- A_stack[, ((i - 1) * k + 1):(i * k), drop = FALSE]
  }

  list(c = c_vec, A = A_list)
}

# -----------------------------------------------------------------------------
# Moving-average / impulse coefficients Psi(h), h = 0..H.
# Psi[h+1, , j] is the response of all variables, h periods after a unit shock
# to variable j, obtained by iterating the companion form.
# -----------------------------------------------------------------------------
compute_ma_coefficients <- function(B, p, H) {
  k <- ncol(B)
  A_comp <- build_companion(B, k, p)

  Psi <- array(
    0.0,
    dim = c(H + 1, k, k),
    dimnames = list(horizon = 0:H, variable = NULL, shock_col = NULL)
  )

  for (j in seq_len(k)) {
    delta <- rep(0, k)
    delta[j] <- 1

    state <- c(delta, rep(0, k * (p - 1)))
    Psi[1, , j] <- delta

    for (hh in 1:H) {
      state <- A_comp %*% state
      Psi[hh + 1, , j] <- state[1:k]
    }
  }

  Psi
}

# -----------------------------------------------------------------------------
# Deterministic baseline forecast path of the VAR (no shocks), horizons 0..H,
# starting from the last p observations of Y_hist.
# -----------------------------------------------------------------------------
forecast_baseline_path <- function(B, Y_hist, p, H) {
  k <- ncol(Y_hist)
  decomp <- extract_c_A_list(B, k, p)

  c_vec <- decomp$c
  A_list <- decomp$A

  lag_list <- lapply(1:p, function(ell) {
    Y_hist[nrow(Y_hist) - ell + 1, , drop = FALSE]
  })

  Y_fore <- matrix(NA_real_, nrow = H + 1, ncol = k)

  for (hh in 0:H) {
    y_next <- c_vec

    for (ell in 1:p) {
      y_next <- y_next + A_list[[ell]] %*% as.numeric(lag_list[[ell]])
    }

    Y_fore[hh + 1, ] <- y_next

    if (p > 1) {
      for (ell in seq(p, 2)) {
        lag_list[[ell]] <- lag_list[[ell - 1]]
      }
    }

    lag_list[[1]] <- matrix(y_next, nrow = 1)
  }

  Y_fore
}

invisible(TRUE)
