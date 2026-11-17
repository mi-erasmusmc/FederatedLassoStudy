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
#' With `scoreConverged`, the converged models of [runConvergedReference()] are
#' scored again from their exported coefficients: the pooled data at DualAvg's
#' penalty, each converged local lasso, and the coefficient average, probability
#' ensemble and largest-database model built from the converged local lassos
#' with the saved models' weights and databases. This recomputes their
#' calibration curves with the installed framework without refitting.
#'
#' The `export` directory holds `referenceScores` and `referenceCurves` in the
#' format of [runConvergedReference()] (the converged pooled lasso first in each
#' evaluation, `kind` converged; every scored model `kind` saved, and the
#' rescored converged models `kind` converged), and
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
#' @param scoreConverged Also score the converged models whose coefficients
#'   `referenceDirectory` exports.
#' @return The export tables, invisibly.
#' @export
runReferenceScoring <- function(dataRoot, siteIds, resultsDirectories, referenceDirectory, outputDirectory,
                                tasks = names(studyTasks()), featureSets = c("ageSexPhenotypes", "ageSexConditions"),
                                holdouts = siteIds, taskDirectories = NULL, scoreConverged = TRUE) {
  exported <- utils::read.csv(file.path(referenceDirectory, "export", "referenceCoefficients.csv"),
    stringsAsFactors = FALSE, colClasses = c(covariateId = "character", sourceDatabase = "character"))
  reference <- exported[exported$model == "convergedPooled", , drop = FALSE]
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
      if (isTRUE(scoreConverged)) {
        converged <- exported[exported$task == task & exported$featureSet == featureSet &
          exported$heldOutDatabase == holdout & exported$model != "convergedPooled", , drop = FALSE]
        candidates <- c(candidates, convergedCandidates(converged, saved, columns, test$xMatrix))
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

# Candidates for scoreAgainstReference() from the exported converged coefficients of one evaluation: the
# pooled data at DualAvg's penalty, the coefficient average, the probability ensemble and largest-database
# model rebuilt from the converged local lassos with the saved models' weights and databases, and each
# converged local lasso, in the order runConvergedReference() scores them.
convergedCandidates <- function(rows, saved, columns, x) {
  coefficientsOf <- function(model, source = NA_character_) {
    keep <- rows$model == model & (if (is.na(source)) is.na(rows$sourceDatabase) | rows$sourceDatabase == "" else
      !is.na(rows$sourceDatabase) & rows$sourceDatabase == source)
    r <- rows[keep, , drop = FALSE]
    if (!nrow(r)) return(NULL)
    if (nrow(r) != length(columns) || !setequal(r$covariateId, columns)) {
      stop("The exported ", model, if (!is.na(source)) paste0(" (", source, ")"),
        " coefficients do not match the scored predictors", call. = FALSE)
    }
    r$coefficient[match(columns, r$covariateId)]
  }
  lp <- function(w) as.numeric(x %*% w)
  candidate <- function(model, linearPredictor, probability = stats::plogis(linearPredictor), source = NA_character_) {
    list(model = model, kind = "converged", source = source, probability = probability,
      linearPredictor = linearPredictor)
  }
  out <- list()
  for (model in c("convergedDualAvgPenalty", "convergedLocalAvgLasso")) {
    w <- coefficientsOf(model)
    if (!is.null(w)) out[[length(out) + 1L]] <- candidate(model, lp(w))
  }
  sources <- unique(rows$sourceDatabase[rows$model == "convergedLocal"])
  local <- stats::setNames(lapply(sources, function(site) coefficientsOf("convergedLocal", site)), sources)
  siteOf <- function(component) as.character(component$sourceClientId)
  if (!is.null(saved$BiggestSiteLasso)) {
    site <- siteOf(saved$BiggestSiteLasso$models[[1]])
    if (site %in% sources) out[[length(out) + 1L]] <- candidate("convergedBiggestSiteLasso", lp(local[[site]]))
  }
  if (!is.null(saved$LocalEnsembleLasso)) {
    models <- saved$LocalEnsembleLasso$models
    sites <- vapply(models, siteOf, character(1))
    if (all(sites %in% sources)) {
      weights <- vapply(models, `[[`, numeric(1), "weight")
      probability <- as.numeric(do.call(cbind, lapply(local[sites], function(w) stats::plogis(lp(w)))) %*% weights)
      out[[length(out) + 1L]] <- candidate("convergedLocalEnsembleLasso",
        stats::qlogis(pmin(pmax(probability, 1e-15), 1 - 1e-15)), probability)
    }
  }
  for (site in sources) out[[length(out) + 1L]] <- candidate("convergedLocal", lp(local[[site]]), source = site)
  out
}
