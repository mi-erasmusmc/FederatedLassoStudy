#' Matched pooled and DualAvg tuning pilot
#'
#' Fits pooled lasso and DualAvg for one task with one database held out. Both
#' select their penalty by the agreed leave-one-training-database-out log-loss
#' search. The `export` directory then holds, without held-out results or
#' patient counts:
#' - `tuningSelected.csv`: per method, the selected prior variance and
#'   penalty, search stop reason and evaluations, final-fit rounds, training
#'   objective, KKT residual and nonzero predictors, and checks that both
#'   methods used the same inner databases, that the selected candidate has the
#'   best mean score, and that each DualAvg inner penalty follows from its
#'   training size;
#' - `tuningTraces.csv`: every candidate's validation score.
#'
#' @inheritParams runStudy
#' @param task A task from [studyTasks()].
#' @param holdout The outer held-out database.
#' @param featureSet A specification from [studyFeatureSettings()].
#' @return The selected-penalty summary, invisibly.
#' @export
runTuningPilot <- function(dataRoot, siteIds, outputDirectory, task = "DvtTetraplegia",
                           holdout = "JMDC", featureSet = "ageSexPhenotypes",
                           taskDirectories = NULL, hosts = "localhost") {
  runStudy(dataRoot, siteIds, outputDirectory, tasks = task, holdouts = holdout,
    methods = c("PooledLasso", "DualAvg"), featureSets = featureSet,
    taskDirectories = taskDirectories, hosts = hosts)
  invisible(writeTuningPilotExport(outputDirectory, task, holdout, featureSet))
}

#' DualAvg stopping-threshold pilot
#'
#' Refits DualAvg at the penalty it selected in [runTuningPilot()], once per
#' stopping threshold, on the outer-training databases only. Each of
#' `innerHoldouts` is held out in turn for validation; the outer held-out
#' database is never loaded. The penalty, training data, initialization and
#' other settings are the same for every threshold. The `export` directory then
#' holds, without patient counts:
#' - `thresholdFits.csv`: per threshold and validation database, rounds, time,
#'   training objective, KKT residual, nonzero predictors and validation
#'   discrimination and calibration;
#' - `thresholdAgreement.csv`: per validation database, each threshold against
#'   the first: coefficient differences, support overlap and differences in
#'   validation predictions.
#'
#' @inheritParams runTuningPilot
#' @param tuningDirectory The [runTuningPilot()] output directory.
#' @param outputDirectory A new result directory.
#' @param innerHoldouts Training databases to hold out for validation.
#' @param epsilons Stopping thresholds; the first is the reference.
#' @return The agreement summary, invisibly.
#' @export
runDualAvgThresholdPilot <- function(dataRoot, siteIds, tuningDirectory, outputDirectory,
                                     task = "DvtTetraplegia", holdout = "JMDC",
                                     innerHoldouts = setdiff(siteIds, holdout),
                                     epsilons = c(1e-10, 1e-8), featureSet = "ageSexPhenotypes",
                                     taskDirectories = NULL, hosts = "localhost") {
  training <- setdiff(siteIds, holdout)
  if (!length(innerHoldouts) || !all(innerHoldouts %in% training)) {
    stop("innerHoldouts must be training databases, not the outer held-out database", call. = FALSE)
  }
  results <- FederatedLearning::readComparisonResults(tuningDirectory)
  tuned <- results[results$task == task & results$method == "DualAvg" &
    results$featureSet == featureSet & results$clientId == holdout, , drop = FALSE]
  if (nrow(tuned) != 1L || !is.finite(tuned$selectedLambda)) {
    stop("No tuned DualAvg penalty for this task, holdout and feature set in tuningDirectory", call. = FALSE)
  }
  for (epsilon in epsilons) {
    runStudy(dataRoot, training, thresholdDirectory(outputDirectory, epsilon), tasks = task,
      holdouts = innerHoldouts, methods = "DualAvg", featureSets = featureSet,
      dualAvgEpsilon = epsilon, dualAvgLambda = tuned$selectedLambda,
      taskDirectories = taskDirectories, hosts = hosts)
  }
  sites <- studySites(dataRoot, training, task, taskDirectories)[[task]]
  invisible(writeThresholdPilotExport(outputDirectory, sites, task, featureSet, epsilons, innerHoldouts))
}

thresholdDirectory <- function(outputDirectory, epsilon) {
  file.path(outputDirectory, paste0("epsilon", format(epsilon, scientific = TRUE)))
}

writeTuningPilotExport <- function(outputDirectory, task, holdout, featureSet) {
  pooled <- FederatedLearning::loadModel(outputDirectory, task, holdout, featureSet, "PooledLasso")
  dual <- FederatedLearning::loadModel(outputDirectory, task, holdout, featureSet, "DualAvg")
  results <- FederatedLearning::readComparisonResults(outputDirectory)
  results <- results[results$task == task & results$featureSet == featureSet &
    results$clientId == holdout, , drop = FALSE]
  pooledTuning <- pooled$models[[1]]$fittingSettings$tuning
  pooledTrace <- pooledTuning$trace
  dualTrace <- dual$config$lambdaSearchTrace
  dualTrace$validationDatabase <- dual$provenance$trainClientIds[dualTrace$validationClient]
  dualMean <- stats::aggregate(score ~ searchValue, dualTrace, mean)

  traces <- rbind(
    data.frame(method = "PooledLasso", iteration = pooledTrace$iteration, variance = pooledTrace$variance,
      validationDatabase = "mean", lambda = NA_real_, rounds = NA_real_, score = pooledTrace$score),
    data.frame(method = "DualAvg", iteration = dualTrace$iteration, variance = dualTrace$searchValue,
      validationDatabase = dualTrace$validationDatabase, lambda = dualTrace$fitLambda,
      rounds = dualTrace$roundsCompleted, score = dualTrace$score))
  # Each inner fit holds out one training database, so its training rows are
  # the total minus that database.
  innerRows <- sum(pooledTuning$databaseSizes) - pooledTuning$databaseSizes
  row <- function(method) results[results$method == method, , drop = FALSE]
  selected <- data.frame(
    method = c("PooledLasso", "DualAvg"),
    selectedVariance = c(pooled$models[[1]]$fittingSettings$selectedVariance, dual$config$selectedVariance),
    # Pooled lasso reports a prior variance; on DualAvg's mean-loss scale its
    # penalty is sqrt(2 / variance) / training rows.
    selectedLambda = c(sqrt(2 / pooled$models[[1]]$fittingSettings$selectedVariance) /
      sum(pooledTuning$databaseSizes), row("DualAvg")$selectedLambda),
    searchStopReason = c(pooledTuning$stopReason, dual$config$lambdaSearchStopReason),
    candidates = c(length(unique(pooledTrace$variance)), length(unique(dualTrace$searchValue))),
    innerDatabases = c(pooledTuning$databases, length(unique(dualTrace$validationClient))),
    sameInnerDatabases = setequal(innerRows, unique(dualTrace$trainingRows)),
    selectedIsBest = c(
      isTRUE(all.equal(pooledTrace$variance[which.max(pooledTrace$score)],
        pooled$models[[1]]$fittingSettings$selectedVariance)),
      isTRUE(all.equal(dualMean$searchValue[which.max(dualMean$score)], dual$config$selectedVariance))),
    innerPenaltiesFollowSizes = c(NA, isTRUE(all.equal(dualTrace$fitLambda,
      sqrt(2 / dualTrace$searchValue) / dualTrace$trainingRows))),
    finalRounds = c(NA, row("DualAvg")$roundsCompleted),
    trainMeanLogLoss = c(row("PooledLasso")$pooledMeanLogLoss, row("DualAvg")$pooledMeanLogLoss),
    objectiveGapToPooled = c(0, row("DualAvg")$pooledObjectiveGap),
    kktMaxAbs = c(row("PooledLasso")$pooledKktMaxAbs, row("DualAvg")$pooledKktMaxAbs),
    nonzeroPredictors = c(row("PooledLasso")$nonzeroPredictors, row("DualAvg")$nonzeroPredictors),
    elapsedSeconds = c(row("PooledLasso")$elapsedSeconds, row("DualAvg")$elapsedSeconds),
    tuningMessages = c(row("PooledLasso")$tuningMessages, row("DualAvg")$tuningMessages),
    stringsAsFactors = FALSE)
  writeExport(outputDirectory, list(tuningSelected = selected, tuningTraces = traces))
  selected
}

writeThresholdPilotExport <- function(outputDirectory, sites, task, featureSet, epsilons, innerHoldouts) {
  taskSettings <- studyTasks()[[task]]
  fitColumns <- c("roundsCompleted", "elapsedSeconds", "selectedLambda", "pooledMeanLogLoss",
    "pooledKktMaxAbs", "nonzeroPredictors", "auc", "logLoss", "calibrationInTheLarge",
    "calibrationSlope", "calibrationIci")
  fits <- list()
  agreement <- list()
  for (heldOut in innerHoldouts) {
    models <- lapply(epsilons, function(epsilon) FederatedLearning::loadModel(
      thresholdDirectory(outputDirectory, epsilon), task, heldOut, featureSet, "DualAvg"))
    for (i in seq_along(epsilons)) {
      results <- FederatedLearning::readComparisonResults(thresholdDirectory(outputDirectory, epsilons[[i]]))
      result <- results[results$task == task & results$featureSet == featureSet &
        results$method == "DualAvg" & results$clientId == heldOut, , drop = FALSE]
      fits[[length(fits) + 1L]] <- cbind(data.frame(epsilon = epsilons[[i]], validationDatabase = heldOut,
        roundLimit = models[[i]]$config$rounds), result[intersect(fitColumns, names(result))])
    }
    predictions <- lapply(models, FederatedLearning::predictModel, sitePath = sites[[heldOut]],
      taskSettings = taskSettings)
    reference <- coefficientVector(models[[1]])
    for (i in seq_along(epsilons)[-1L]) {
      coefficients <- coefficientVector(models[[i]])
      if (!identical(names(coefficients), names(reference))) {
        stop("Threshold fits for ", heldOut, " use different feature maps", call. = FALSE)
      }
      support <- function(w) names(w)[w != 0 & names(w) != "(Intercept)"]
      union <- union(support(reference), support(coefficients))
      lp <- predictions[[i]]$linearPredictor - predictions[[1]]$linearPredictor
      probability <- abs(predictions[[i]]$probability - predictions[[1]]$probability)
      agreement[[length(agreement) + 1L]] <- data.frame(validationDatabase = heldOut,
        epsilon = epsilons[[i]], referenceEpsilon = epsilons[[1]],
        coefficientMaxAbsDiff = max(abs(coefficients - reference)),
        coefficientRmsDiff = sqrt(mean((coefficients - reference)^2)),
        nonzero = length(support(coefficients)), referenceNonzero = length(support(reference)),
        supportJaccard = if (length(union)) length(intersect(support(reference), support(coefficients))) /
          length(union) else 1,
        linearPredictorMaxAbsDiff = max(abs(lp)),
        # NA when a fit has constant predictions.
        linearPredictorCorrelation = suppressWarnings(
          stats::cor(predictions[[i]]$linearPredictor, predictions[[1]]$linearPredictor)),
        probabilityMaxAbsDiff = max(probability),
        probabilityP99AbsDiff = unname(stats::quantile(probability, 0.99)))
    }
  }
  agreement <- do.call(rbind, agreement)
  writeExport(outputDirectory, list(thresholdFits = do.call(rbind, fits), thresholdAgreement = agreement))
  agreement
}

# Shared-scale coefficients named by covariate id, intercept first.
coefficientVector <- function(model) {
  table <- model$models[[1]]$coefficients
  table <- table[order(table$matrixColumn), , drop = FALSE]
  stats::setNames(table$coefficientSharedScale, ifelse(table$isIntercept, "(Intercept)", table$covariateId))
}

writeExport <- function(outputDirectory, tables) {
  directory <- file.path(outputDirectory, "export")
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  for (name in names(tables)) {
    utils::write.csv(tables[[name]], file.path(directory, paste0(name, ".csv")), row.names = FALSE)
  }
}
