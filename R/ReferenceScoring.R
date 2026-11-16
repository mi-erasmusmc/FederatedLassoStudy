#' Score saved models against the converged pooled lasso
#'
#' Scores the saved models of `resultsDirectories` in each held-out database
#' against the converged pooled lasso of [runConvergedReference()], whose
#' coefficients are read from `referenceDirectory`. It needs only the held-out
#' databases' data, because the reference is not refitted. The measures are
#' those of [runConvergedReference()]: AUC with DeLong's 95% interval, the
#' paired difference from the reference, model-based AUC, linear-predictor
#' correlation, absolute differences and top 10% and top 1% overlap of
#' predicted risk, calibration-in-the-large, slope, the calibration curve and
#' the integrated calibration index.
#'
#' The `export` directory holds `referenceScores` and `referenceCurves` in the
#' format of [runConvergedReference()] (the converged pooled lasso first in each
#' evaluation, `kind` converged; every scored model `kind` saved), and
#' `scoredCoefficients`, the scored models' coefficients by covariate id. Only
#' aggregate tables are written.
#'
#' @param dataRoot,siteIds,taskDirectories As for [runStudy()].
#' @param resultsDirectories Result directories of the models to score.
#' @param referenceDirectory Output directory of [runConvergedReference()],
#'   holding `export/referenceCoefficients.csv`.
#' @param outputDirectory Directory for the `export` tables.
#' @param tasks,featureSets,holdouts The evaluations to score; evaluations
#'   without a saved model are skipped.
#' @return The export tables, invisibly.
#' @export
runReferenceScoring <- function(dataRoot, siteIds, resultsDirectories, referenceDirectory, outputDirectory,
                                tasks = names(studyTasks()), featureSets = c("ageSexPhenotypes", "ageSexConditions"),
                                holdouts = siteIds, taskDirectories = NULL) {
  reference <- utils::read.csv(file.path(referenceDirectory, "export", "referenceCoefficients.csv"),
    stringsAsFactors = FALSE, colClasses = c(covariateId = "character"))
  reference <- reference[reference$model == "convergedPooled", , drop = FALSE]
  available <- scoredModels(resultsDirectories)
  tables <- list(referenceScores = list(), referenceCurves = list(), scoredCoefficients = list())
  for (featureSet in featureSets) for (task in tasks) {
    rows <- available[available$task == task & available$featureSet == featureSet, , drop = FALSE]
    if (!nrow(rows)) next
    sites <- studySites(dataRoot, siteIds, task, featureSet, taskDirectories)[[task]]
    population <- studyTasks()[[task]]$populationSettings
    for (holdout in intersect(holdouts, unique(rows$heldOutDatabase))) {
      message(format(Sys.time(), "%H:%M"), "  scoring ", task, ", ", featureSet, ", ", holdout, " held out")
      selected <- rows[rows$heldOutDatabase == holdout & !duplicated(rows[c("heldOutDatabase", "method")]), , drop = FALSE]
      saved <- stats::setNames(lapply(seq_len(nrow(selected)), function(i) {
        FederatedLearning::loadModel(selected$directory[[i]], task, holdout, featureSet, selected$method[[i]])
      }), selected$method)
      key <- data.frame(task = task, featureSet = featureSet, heldOutDatabase = holdout, stringsAsFactors = FALSE)
      w <- reference[reference$task == task & reference$featureSet == featureSet & reference$heldOutDatabase == holdout, ]
      if (!nrow(w)) stop("No converged pooled lasso for ", task, " with ", featureSet, " and ", holdout, " held out",
        call. = FALSE)
      config <- saved[[1]]$config
      # Held-out entries are all scored, as FederatedLearning::predictModel() scores them.
      config$diagnosticControlsPerCase <- NULL
      plpData <- FederatedLearning::loadClientData(sites[[holdout]], population)
      data <- FederatedLearning::createClientMatrix(plpData, config)
      test <- list(xMatrix = data$xMatrix, yLabels = data$yLabels)
      rm(plpData, data)
      columns <- names(componentCoefficients(saved[[1]]$models[[1]]))
      if (length(columns) != ncol(test$xMatrix) || !setequal(columns, w$covariateId)) {
        stop("The scored models and the converged pooled lasso must use the same predictors (", task, ", ",
          featureSet, ", ", holdout, ")", call. = FALSE)
      }
      referenceW <- w$coefficient[match(columns, w$covariateId)]
      eta <- as.numeric(test$xMatrix %*% referenceW)
      candidates <- list(list(model = "convergedPooled", kind = "converged", source = NA_character_,
        probability = stats::plogis(eta), linearPredictor = eta))
      for (method in names(saved)) {
        if (!identical(names(componentCoefficients(saved[[method]]$models[[1]])), columns)) {
          stop(method, " does not use the same predictors as the other scored models", call. = FALSE)
        }
        predictions <- modelPredictions(saved[[method]], test)
        candidates[[length(candidates) + 1L]] <- list(model = method, kind = "saved", source = NA_character_,
          probability = predictions$probability, linearPredictor = predictions$linearPredictor)
        for (component in saved[[method]]$models) {
          coefficients <- componentCoefficients(component)
          tables$scoredCoefficients[[length(tables$scoredCoefficients) + 1L]] <- cbind(key, model = method,
            covariateId = names(coefficients), coefficient = unname(coefficients), stringsAsFactors = FALSE)
        }
      }
      scored <- scoreAgainstReference(test$yLabels, candidates)
      tables$referenceScores[[length(tables$referenceScores) + 1L]] <- cbind(key, scored$scores)
      if (!is.null(scored$curves) && nrow(scored$curves)) {
        tables$referenceCurves[[length(tables$referenceCurves) + 1L]] <- cbind(key, scored$curves)
      }
      rm(test, candidates, scored, saved)
      gc()
    }
  }
  out <- lapply(tables, function(rows) if (length(rows)) do.call(rbind, rows) else data.frame())
  writeExport(outputDirectory, out)
  invisible(out)
}
