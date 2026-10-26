#' Study methods
#'
#' The workflows of the final comparison: pooled lasso as the reference,
#' DualAvg, ODAL1/ODAL2, ADAP1/ADAP2, sample-size-weighted coefficient
#' averaging and probability ensembling, the largest training database's
#' model, and each training database's own model.
#'
#' @return Method names.
#' @export
studyMethods <- function() {
  c("PooledLasso", "DualAvg", "ODAL1", "ODAL2", "ADAP1", "ADAP2", "LocalAvgLasso",
    "LocalEnsembleLasso", "BiggestSiteLasso", "LocalSiteLasso")
}

#' Study method settings
#'
#' Penalty tuning uses the framework defaults agreed for the final run: pooled
#' lasso and DualAvg select their penalties by leave-one-training-database-out
#' log loss, and local lasso and ADAP (at the lead site) by outcome-stratified
#' three-fold CV within the database. Settings the framework does not default
#' to are those of the previous production run: DualAvg with one local step
#' per round (`k = 1`), unit step sizes and at most 10,000 rounds; ODAL
#' initialization with ridge fallback (lambda 1e-8); and a 1e-10 convergence
#' threshold for ODAL and ADAP.
#'
#' @param methods Methods from [studyMethods()].
#' @param dualAvgEpsilon DualAvg stopping threshold on the successive change in
#'   loss. 1e-10 as before until the threshold pilot decides.
#' @param dualAvgLambda A fixed DualAvg penalty instead of tuning, for pilots.
#' @param dualAvgRounds DualAvg round limit; 10,000 as before, higher only to
#'   check how far round-limited fits are from converged.
#' @return A list of `methodSettings` objects.
#' @export
studyMethodSettings <- function(methods = studyMethods(), dualAvgEpsilon = 1e-10,
                                dualAvgLambda = NULL, dualAvgRounds = 10000L) {
  if (!length(methods) || !all(methods %in% studyMethods())) {
    stop("methods must name studyMethods(): ", paste(studyMethods(), collapse = ", "), call. = FALSE)
  }
  lapply(methods, function(method) {
    options <- switch(method,
      DualAvg = c(list(k = 1L, etaClient = 1, etaServer = 1, epsilon = dualAvgEpsilon,
        dualavgRounds = as.integer(dualAvgRounds)),
        if (!is.null(dualAvgLambda)) list(dualavgTuneLambda = FALSE, dualavgLambda = dualAvgLambda)),
      ODAL1 = , ODAL2 = list(odalInit = "ridgeFallback", odalRidgeLambda = 1e-8, epsilon = 1e-10),
      ADAP1 = , ADAP2 = list(epsilon = 1e-10),
      list()
    )
    do.call(FederatedLearning::createMethodSettings, c(list(method), options))
  })
}
