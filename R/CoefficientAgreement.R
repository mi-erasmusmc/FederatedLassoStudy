#' Coefficient agreement with pooled lasso
#'
#' Compares each saved model that has a single coefficient vector with the
#' PooledLasso model of the same run, task, predictor set and held-out
#' database. Each local model is compared separately, by source database. The
#' probability ensemble has no single coefficient vector and is skipped. Only
#' saved models are read, so no site data or workers are needed.
#'
#' Coefficients are aligned by covariate id on the shared matrix scale. Age and
#' the intercept are left out, so the compared coefficients are those of the
#' binary predictors (phenotypes or conditions, and the female indicator).
#' Predictors that one model's feature filter removed are counted, not imputed.
#'
#' The `export` directory then holds:
#' - `coefficientAgreement.csv`: one row per compared model, with predictors in
#'   both models and in only one, Pearson correlation and root mean squared and
#'   largest absolute differences over the predictors in both, exact nonzero
#'   counts, and the Jaccard overlap of the nonzero predictor sets (`NA` when
#'   neither model has any);
#' - `coefficients.csv`: every compared model's coefficients, PooledLasso's
#'   included, by covariate id, with age and the intercept.
#'
#' @param resultsDirectories Result directories of [runStudy()]; their folds
#'   are compared within each directory.
#' @param outputDirectory Directory for the `export` tables.
#' @return The agreement table, invisibly.
#' @export
runCoefficientAgreement <- function(resultsDirectories, outputDirectory) {
  models <- do.call(rbind, lapply(resultsDirectories, savedModels))
  coefficients <- list()
  agreement <- list()
  for (key in unique(models$key)) {
    fold <- models[models$key == key, , drop = FALSE]
    if (!"PooledLasso" %in% fold$method) {
      next
    }
    loaded <- lapply(seq_len(nrow(fold)), function(i) FederatedLearning::loadModel(fold$directory[[i]],
      fold$task[[i]], fold$clientId[[i]], fold$featureSet[[i]], fold$method[[i]]))
    components <- do.call(rbind, Map(modelComponents, loaded, fold$method))
    reference <- components$coefficients[[which(components$method == "PooledLasso")]]
    for (i in seq_len(nrow(components))) {
      w <- components$coefficients[[i]]
      labels <- data.frame(resultsDirectory = basename(fold$directory[[1]]), task = fold$task[[1]],
        featureSet = fold$featureSet[[1]], fold = fold$fold[[1]], heldOutDatabase = fold$clientId[[1]],
        method = components$method[[i]], sourceDatabase = components$sourceDatabase[[i]],
        stringsAsFactors = FALSE)
      coefficients[[length(coefficients) + 1L]] <- cbind(labels, covariateId = names(w),
        coefficient = unname(w), row.names = NULL)
      if (components$method[[i]] != "PooledLasso") {
        agreement[[length(agreement) + 1L]] <- cbind(labels, coefficientAgreement(w, reference))
      }
    }
  }
  if (!length(agreement)) {
    stop("No saved model has a PooledLasso model to compare with", call. = FALSE)
  }
  agreement <- do.call(rbind, agreement)
  writeExport(outputDirectory, list(coefficientAgreement = agreement,
    coefficients = do.call(rbind, coefficients)))
  invisible(agreement)
}

# One row per saved model with a single coefficient vector per component.
savedModels <- function(directory) {
  results <- FederatedLearning::readComparisonResults(directory)
  ok <- results$method != "LocalEnsembleLasso" & !is.na(results$modelFile) & nzchar(results$modelFile)
  if ("error" %in% names(results)) {
    ok <- ok & (is.na(results$error) | !nzchar(results$error))
  }
  models <- unique(results[ok, c("task", "featureSet", "fold", "clientId", "method"), drop = FALSE])
  models$directory <- rep(directory, nrow(models))
  models$key <- paste(directory, models$task, models$featureSet, models$fold, sep = "\r")
  models
}

# Each component's coefficients named by covariate id, the intercept as
# "(Intercept)"; local models carry their source database.
modelComponents <- function(model, method) {
  if (method != "LocalSiteLasso" && length(model$models) != 1L) {
    stop(method, " has ", length(model$models), " coefficient vectors; expected one", call. = FALSE)
  }
  out <- data.frame(method = rep(method, length(model$models)),
    sourceDatabase = vapply(model$models, function(m) {
      if (is.null(m$sourceClientId)) NA_character_ else as.character(m$sourceClientId)
    }, ""),
    stringsAsFactors = FALSE)
  out$coefficients <- lapply(model$models, componentCoefficients)
  out
}

coefficientAgreement <- function(w, reference) {
  drop <- c("(Intercept)", as.character(ageCovariateId))
  w <- w[!names(w) %in% drop]
  reference <- reference[!names(reference) %in% drop]
  shared <- intersect(names(w), names(reference))
  a <- w[shared]
  b <- reference[shared]
  support <- function(x) names(x)[x != 0]
  union <- union(support(w), support(reference))
  data.frame(
    sharedPredictors = length(shared),
    onlyModelPredictors = length(setdiff(names(w), shared)),
    onlyPooledPredictors = length(setdiff(names(reference), shared)),
    # NA when either vector is constant over the shared predictors.
    correlation = if (length(shared) > 1L && stats::sd(a) > 0 && stats::sd(b) > 0) stats::cor(a, b) else NA_real_,
    rmsDifference = if (length(shared)) sqrt(mean((a - b)^2)) else NA_real_,
    maxAbsDifference = if (length(shared)) max(abs(a - b)) else NA_real_,
    maxAbsDifferenceCovariateId = if (length(shared)) shared[[which.max(abs(a - b))]] else NA_character_,
    nonzero = length(support(w)),
    pooledNonzero = length(support(reference)),
    bothNonzero = length(intersect(support(w), support(reference))),
    supportJaccard = if (length(union)) length(intersect(support(w), support(reference))) / length(union) else NA_real_,
    stringsAsFactors = FALSE
  )
}
