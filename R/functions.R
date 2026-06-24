# =============================================================================
#  R/functions.R — Fonctions utilitaires partagées (dédupliquées)
#
#  Ces fonctions étaient copiées-collées à l'identique dans plusieurs scripts.
#  Elles sont désormais regroupées ici et chargées une seule fois via
#  code/00_setup.R.
#
#  NOTE de transition : pour ne rien casser, les scripts d'analyse conservent
#  pour l'instant leurs définitions locales (qui priment). L'étape de
#  refactoring recommandée consiste à supprimer ces copies locales et à
#  s'appuyer uniquement sur ce fichier.
# =============================================================================

# --- Dates : yearqtr -> année décimale (ex. 2001Q3 -> 2001.50) --------------
convert_qtr_to_decimal <- function(qtr) {
  year <- floor(as.numeric(format(qtr, "%Y")))
  qnum <- as.numeric(cycle(qtr))            # 1 à 4
  year + (qnum - 1) * 0.25
}

# --- Stationnarité : ADF + KPSS sur une série, conclusion croisée -----------
test_stationarity <- function(ts_vector, name = NULL) {
  ts_clean <- stats::na.omit(ts_vector)
  adf  <- tryCatch(tseries::adf.test(ts_clean, k = 2, alternative = "stationary"),
                   error = function(e) NULL)
  kpss <- tryCatch(tseries::kpss.test(ts_clean, null = "Level"),
                   error = function(e) NULL)

  conclusion <- if (!is.null(adf) && !is.null(kpss)) {
    if (adf$p.value < 0.05 && kpss$p.value >= 0.05) "Stationnaire"
    else if (adf$p.value >= 0.05 && kpss$p.value < 0.05) "Non stationnaire"
    else "Ambigu"
  } else "Erreur"

  list(
    Variable   = name,
    ADF_Stat   = if (!is.null(adf))  round(adf$statistic, 3)  else NA,
    ADF_p      = if (!is.null(adf))  round(adf$p.value, 3)    else NA,
    KPSS_Stat  = if (!is.null(kpss)) round(kpss$statistic, 3) else NA,
    KPSS_p     = if (!is.null(kpss)) round(kpss$p.value, 3)   else NA,
    Conclusion = conclusion
  )
}

# --- Merton–Vasicek : reconstruction du facteur systémique Z ----------------
# Inverse la formule de Vasicek sur les taux de défaut observés et calibre
# rho de sorte que Var(Z) = 1 (cf. eq. 4.4 du papier).
f_Z_estimation <- function(vect_TD) {
  var_Z_rho <- function(rho) {
    if (rho <= 0 || rho >= 1) return(Inf)
    Z_estim <- (qnorm(mean(vect_TD, na.rm = TRUE)) -
                  qnorm(vect_TD) * sqrt(1 - rho)) / sqrt(rho)
    var(Z_estim, na.rm = TRUE) - 1
  }
  result <- tryCatch(uniroot(var_Z_rho, lower = 1e-6, upper = 1 - 1e-6),
                     error = function(e) NULL)
  if (is.null(result)) return(list(Z = rep(NA, length(vect_TD)), rho = NA))
  rho_opt <- result$root
  Z_estim <- (qnorm(mean(vect_TD, na.rm = TRUE)) -
                qnorm(vect_TD) * sqrt(1 - rho_opt)) / sqrt(rho_opt)
  list(Z = Z_estim, rho = rho_opt)
}

# --- Merton–Vasicek : PD point-in-time conditionnelle à Z (eq. 4.2) ---------
pd_vasicek <- function(Z, p = PD_UNCOND, rho = RHO_ASSET) {
  pnorm((qnorm(p) - sqrt(rho) * Z) / sqrt(1 - rho))
}

# --- Nettoyage de l'environnement parallèle (foreach) -----------------------
unregister_dopar <- function() {
  if ("foreach" %in% rownames(installed.packages())) {
    env <- foreach:::.foreachGlobals
    rm(list = ls(name = env), pos = env)
  }
}

invisible(TRUE)
