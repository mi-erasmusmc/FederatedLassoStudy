#' Calibration curve export
#'
#' Collects the calibration curves saved with every model of one or more result
#' directories, including each local model by source database. Each curve is
#' the framework's smoothed observed risk against predicted risk at 50
#' percentiles (1st to 99th) of the held-out predictions, with pointwise 95%
#' intervals, so it holds no patient-level values. Only saved models are read.
#'
#' The `export` directory then holds `calibrationCurves.csv`: one row per curve
#' point with the results directory, task, predictor set, fold, held-out and
#' source database, method, and predicted, observed, lower and upper risk.
#'
#' @param resultsDirectories Result directories of [runStudy()].
#' @param outputDirectory Directory for the `export` table.
#' @return The curve table, invisibly.
#' @export
runCalibrationCurveExport <- function(resultsDirectories, outputDirectory) {
  curves <- list()
  for (directory in resultsDirectories) {
    results <- FederatedLearning::readComparisonResults(directory)
    ok <- !is.na(results$modelFile) & nzchar(results$modelFile)
    if ("error" %in% names(results)) {
      ok <- ok & (is.na(results$error) | !nzchar(results$error))
    }
    models <- unique(results[ok, c("task", "featureSet", "fold", "clientId", "method"), drop = FALSE])
    for (i in seq_len(nrow(models))) {
      model <- FederatedLearning::loadModel(directory, models$task[[i]], models$clientId[[i]],
        models$featureSet[[i]], models$method[[i]])
      curve <- model$calibrationCurves
      if (is.null(curve) || !nrow(curve)) {
        next
      }
      source <- if ("sourceClientId" %in% names(curve)) as.character(curve$sourceClientId) else NA_character_
      curves[[length(curves) + 1L]] <- data.frame(resultsDirectory = basename(directory),
        task = models$task[[i]], featureSet = models$featureSet[[i]], fold = models$fold[[i]],
        heldOutDatabase = models$clientId[[i]], method = models$method[[i]], sourceDatabase = source,
        predicted = curve$predicted, observed = curve$observed, lower = curve$lower, upper = curve$upper,
        stringsAsFactors = FALSE)
    }
  }
  if (!length(curves)) {
    stop("No saved model has a calibration curve", call. = FALSE)
  }
  curves <- do.call(rbind, curves)
  writeExport(outputDirectory, list(calibrationCurves = curves))
  invisible(curves)
}
