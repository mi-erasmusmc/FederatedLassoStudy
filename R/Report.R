#' Study report tables
#'
#' Builds the manuscript tables from the aggregate result files of the final
#' run, so it runs anywhere the exported CSVs are. Paired differences subtract
#' PooledLasso's value within the same task, predictor set and held-out
#' database. Task means weight the five databases equally and overall means
#' weight every task and database evaluation equally.
#'
#' Fits are flagged as agreed for the final run: a fit fails its own check
#' (ADAP KKT acceptance, ODAL `optim` code 0), or its pooled training log loss
#' exceeds PooledLasso's by more than `trainingLossMargin` nats per entry.
#' Flagged fits stay in every table; summaries also give means without them.
#' DualAvg has a fixed round limit and the baselines no iterative check, so
#' their check is `NA`.
#'
#' Tables written to `outputDirectory`:
#' - `population.csv`: cohort entries and outcomes by task and database;
#' - `predictors.csv`: predictors retained by the shared filter, intercept
#'   excluded, by task, predictor set and held-out database;
#' - `evaluations.csv`: one row per evaluated model with AUC, its difference
#'   from PooledLasso, calibration, check and flag;
#' - `discriminationSummary.csv`: mean AUC and paired difference by method,
#'   predictor set and task, and overall, with evaluable counts;
#' - `calibrationSummary.csv`: median calibration-in-the-large and slope by
#'   method, predictor set and task, and overall;
#' - `localBaselines.csv`: per evaluation, the spread of single-source AUCs
#'   beside the aggregated, largest-site and pooled models;
#' - `communication.csv`: median final-fit and tuning communication by method
#'   and predictor set;
#' - `coefficientSummary.csv` and `localCoefficientSummary.csv`, when
#'   `coefficientAgreementDirectory` is given: median agreement with
#'   PooledLasso and evaluable counts by method or source database;
#' - `roundsCheck.csv`, when `roundsCheckDirectories` are given.
#'
#' @param resultsDirectories Result directories of [runStudy()].
#' @param outputDirectory Directory for the tables.
#' @param coefficientAgreementDirectory Output of [runCoefficientAgreement()].
#' @param roundsCheckDirectories Outputs of [runDualAvgRoundsCheck()].
#' @param trainingLossMargin Flag threshold on the training log loss gap.
#' @return The evaluations table, invisibly.
#' @export
writeStudyReport <- function(resultsDirectories, outputDirectory, coefficientAgreementDirectory = NULL,
                             roundsCheckDirectories = NULL, trainingLossMargin = 0.05) {
  results <- do.call(rbind, lapply(resultsDirectories, function(directory) {
    reportColumns(FederatedLearning::readComparisonResults(directory))
  }))
  key <- c("task", "featureSet", "fold", "method", "configLabel")
  if (anyDuplicated(results[key])) {
    stop("Some evaluations appear in more than one results directory", call. = FALSE)
  }
  evaluations <- evaluationTable(results, trainingLossMargin)
  tables <- list(
    population = populationTable(evaluations),
    predictors = predictorTable(evaluations),
    evaluations = evaluations,
    discriminationSummary = discriminationSummary(evaluations),
    calibrationSummary = calibrationSummary(evaluations),
    localBaselines = localBaselineTable(evaluations),
    communication = communicationTable(results)
  )
  if (!is.null(coefficientAgreementDirectory)) {
    agreement <- utils::read.csv(file.path(coefficientAgreementDirectory, "export", "coefficientAgreement.csv"),
      stringsAsFactors = FALSE, colClasses = c(maxAbsDifferenceCovariateId = "character"))
    tables$coefficientSummary <- coefficientSummary(agreement[is.na(agreement$sourceDatabase) |
      agreement$method == "BiggestSiteLasso", , drop = FALSE], "method")
    tables$localCoefficientSummary <- coefficientSummary(agreement[agreement$method == "LocalSiteLasso", ,
      drop = FALSE], "sourceDatabase")
  }
  if (length(roundsCheckDirectories)) {
    tables$roundsCheck <- do.call(rbind, lapply(roundsCheckDirectories, function(directory) {
      fits <- utils::read.csv(file.path(directory, "export", "roundsFits.csv"), stringsAsFactors = FALSE)
      cbind(featureSet = FederatedLearning::readComparisonResults(directory)$featureSet[[1]], fits)
    }))
  }
  dir.create(outputDirectory, recursive = TRUE, showWarnings = FALSE)
  for (name in names(tables)) {
    utils::write.csv(tables[[name]], file.path(outputDirectory, paste0(name, ".csv")), row.names = FALSE)
  }
  invisible(evaluations)
}

# Columns the report reads; diagnostics a method does not produce become NA.
reportColumns <- function(results) {
  columns <- c("task", "featureSet", "fold", "clientId", "method", "configLabel", "error", "n", "outcomes",
    "p", "nonzeroPredictors", "auc", "logLoss", "calibrationInTheLarge", "calibrationSlope", "calibrationIci",
    "roundsCompleted", "optimConvergence", "adapFinalConverged", "pooledObjectiveGap",
    "dualAvgReferenceObjectiveGap", "communicationBasis", "exchanges", "messages", "numbers", "bytes",
    "tuningExchanges", "tuningMessages", "tuningNumbers", "tuningBytes", "elapsedSeconds")
  for (column in setdiff(columns, names(results))) {
    results[[column]] <- rep(NA, nrow(results))
  }
  results[columns]
}

evaluationTable <- function(results, trainingLossMargin) {
  failed <- !is.na(results$error) & nzchar(results$error)
  out <- data.frame(task = results$task, featureSet = results$featureSet, fold = results$fold,
    heldOutDatabase = results$clientId, method = results$method,
    sourceDatabase = ifelse(startsWith(results$configLabel, "source="), sub("^source=", "", results$configLabel),
      NA_character_),
    available = !failed & is.finite(results$auc), n = results$n, outcomes = results$outcomes,
    retainedPredictors = results$p - 1L, nonzeroPredictors = results$nonzeroPredictors,
    auc = results$auc, stringsAsFactors = FALSE)
  pooled <- results[results$method == "PooledLasso", c("task", "featureSet", "fold", "auc",
    "calibrationInTheLarge", "calibrationSlope")]
  names(pooled)[4:6] <- c("pooledAuc", "pooledCalibrationInTheLarge", "pooledCalibrationSlope")
  out <- cbind(out, pooled[match(paste(out$task, out$featureSet, out$fold),
    paste(pooled$task, pooled$featureSet, pooled$fold)), 4:6])
  out$aucDifference <- out$auc - out$pooledAuc
  out$logLoss <- results$logLoss
  out$calibrationInTheLarge <- results$calibrationInTheLarge
  out$calibrationSlope <- results$calibrationSlope
  out$calibrationIci <- results$calibrationIci
  out$check <- ifelse(results$method %in% c("ADAP1", "ADAP2"), results$adapFinalConverged == 1,
    ifelse(results$method %in% c("ODAL1", "ODAL2"), results$optimConvergence == 0, NA))
  out$trainingLogLossGap <- results$pooledObjectiveGap
  out$roundsCompleted <- results$roundsCompleted
  out$dualAvgObjectiveGap <- results$dualAvgReferenceObjectiveGap
  lossFlag <- !is.na(out$trainingLogLossGap) & out$trainingLogLossGap > trainingLossMargin
  checkFlag <- !is.na(out$check) & !out$check
  out$flagged <- checkFlag | lossFlag
  out$flagReason <- ifelse(checkFlag & lossFlag, "check failed; training log loss",
    ifelse(checkFlag, "check failed", ifelse(lossFlag, "training log loss", "")))
  rownames(out) <- NULL
  out
}

populationTable <- function(evaluations) {
  rows <- evaluations[evaluations$method == "PooledLasso", c("task", "heldOutDatabase", "n", "outcomes")]
  rows <- unique(rows)
  if (anyDuplicated(rows[c("task", "heldOutDatabase")])) {
    stop("Population counts differ between predictor sets", call. = FALSE)
  }
  names(rows) <- c("task", "database", "entries", "outcomes")
  rows$outcomePercent <- 100 * rows$outcomes / rows$entries
  rows[order(rows$task, rows$database), ]
}

predictorTable <- function(evaluations) {
  rows <- unique(evaluations[evaluations$method == "PooledLasso",
    c("task", "featureSet", "heldOutDatabase", "retainedPredictors")])
  rows[order(rows$task, rows$featureSet, rows$heldOutDatabase), ]
}

# Shared models only: each local source model is reported in localBaselines.
sharedEvaluations <- function(evaluations) {
  evaluations[evaluations$method != "LocalSiteLasso" & evaluations$available, , drop = FALSE]
}

discriminationSummary <- function(evaluations) {
  rows <- sharedEvaluations(evaluations)
  summarize <- function(x) {
    paired <- is.finite(x$aucDifference)
    kept <- paired & !x$flagged
    data.frame(evaluations = nrow(x), meanAuc = mean(x$auc), pairs = sum(paired),
      meanAucDifference = mean(x$aucDifference[paired]),
      minAucDifference = min(x$aucDifference[paired]), maxAucDifference = max(x$aucDifference[paired]),
      flagged = sum(x$flagged), unflaggedPairs = sum(kept),
      unflaggedMeanAucDifference = if (any(kept)) mean(x$aucDifference[kept]) else NA_real_)
  }
  bySummary(rows, summarize)
}

calibrationSummary <- function(evaluations) {
  rows <- sharedEvaluations(evaluations)
  summarize <- function(x) {
    data.frame(evaluations = nrow(x),
      medianCalibrationInTheLarge = stats::median(x$calibrationInTheLarge, na.rm = TRUE),
      minCalibrationInTheLarge = min(x$calibrationInTheLarge, na.rm = TRUE),
      maxCalibrationInTheLarge = max(x$calibrationInTheLarge, na.rm = TRUE),
      medianCalibrationSlope = stats::median(x$calibrationSlope, na.rm = TRUE),
      minCalibrationSlope = min(x$calibrationSlope, na.rm = TRUE),
      maxCalibrationSlope = max(x$calibrationSlope, na.rm = TRUE),
      medianCalibrationIci = stats::median(x$calibrationIci, na.rm = TRUE))
  }
  bySummary(rows, summarize)
}

# Applies summarize per method, predictor set and task, then over all tasks.
bySummary <- function(rows, summarize) {
  groups <- unique(rows[c("method", "featureSet")])
  out <- list()
  for (i in seq_len(nrow(groups))) {
    group <- rows[rows$method == groups$method[[i]] & rows$featureSet == groups$featureSet[[i]], , drop = FALSE]
    for (task in c(sort(unique(group$task)), "overall")) {
      x <- if (task == "overall") group else group[group$task == task, , drop = FALSE]
      out[[length(out) + 1L]] <- cbind(groups[i, , drop = FALSE], task = task, summarize(x), row.names = NULL)
    }
  }
  do.call(rbind, out)
}

localBaselineTable <- function(evaluations) {
  local <- evaluations[evaluations$method == "LocalSiteLasso" & evaluations$available, , drop = FALSE]
  keys <- unique(evaluations[c("task", "featureSet", "fold", "heldOutDatabase")])
  out <- lapply(seq_len(nrow(keys)), function(i) {
    same <- function(x) x$task == keys$task[[i]] & x$featureSet == keys$featureSet[[i]] & x$fold == keys$fold[[i]]
    sources <- local[same(local), , drop = FALSE]
    auc <- function(method) {
      x <- evaluations[same(evaluations) & evaluations$method == method & evaluations$available, , drop = FALSE]
      if (nrow(x) == 1L) x$auc else NA_real_
    }
    biggest <- evaluations[same(evaluations) & evaluations$method == "BiggestSiteLasso", , drop = FALSE]
    cbind(keys[i, , drop = FALSE], sourceModels = nrow(sources),
      minSourceAuc = if (nrow(sources)) min(sources$auc) else NA_real_,
      maxSourceAuc = if (nrow(sources)) max(sources$auc) else NA_real_,
      bestSourceDatabase = if (nrow(sources)) sources$sourceDatabase[[which.max(sources$auc)]] else NA_character_,
      coefficientAverageAuc = auc("LocalAvgLasso"), probabilityEnsembleAuc = auc("LocalEnsembleLasso"),
      largestSiteAuc = auc("BiggestSiteLasso"),
      largestSiteDatabase = if (nrow(biggest) == 1L) biggest$sourceDatabase else NA_character_,
      pooledAuc = auc("PooledLasso"), row.names = NULL)
  })
  out <- do.call(rbind, out)
  out[order(out$task, out$featureSet, out$fold), ]
}

communicationTable <- function(results) {
  rows <- results[results$method != "LocalSiteLasso", , drop = FALSE]
  columns <- c("exchanges", "messages", "numbers", "bytes", "tuningExchanges", "tuningMessages",
    "tuningNumbers", "tuningBytes", "elapsedSeconds")
  groups <- unique(rows[c("method", "featureSet")])
  out <- lapply(seq_len(nrow(groups)), function(i) {
    x <- rows[rows$method == groups$method[[i]] & rows$featureSet == groups$featureSet[[i]], , drop = FALSE]
    medians <- lapply(x[columns], stats::median, na.rm = TRUE)
    names(medians) <- paste0("median", toupper(substring(columns, 1, 1)), substring(columns, 2))
    cbind(groups[i, , drop = FALSE], basis = paste(sort(unique(x$communicationBasis)), collapse = "; "),
      evaluations = nrow(x), as.data.frame(medians), row.names = NULL)
  })
  do.call(rbind, out)
}

coefficientSummary <- function(agreement, by) {
  measures <- c("correlation", "rmsDifference", "supportJaccard", "nonzero", "pooledNonzero")
  groups <- unique(agreement[c(by, "featureSet")])
  out <- lapply(seq_len(nrow(groups)), function(i) {
    x <- agreement[agreement[[by]] == groups[[by]][[i]] & agreement$featureSet == groups$featureSet[[i]], ,
      drop = FALSE]
    medians <- lapply(x[measures], stats::median, na.rm = TRUE)
    names(medians) <- paste0("median", toupper(substring(measures, 1, 1)), substring(measures, 2))
    pairs <- lapply(x[measures[1:3]], function(v) sum(!is.na(v)))
    names(pairs) <- paste0(measures[1:3], "Pairs")
    cbind(groups[i, , drop = FALSE], comparisons = nrow(x),
      medianSharedPredictors = stats::median(x$sharedPredictors), as.data.frame(medians), as.data.frame(pairs),
      row.names = NULL)
  })
  do.call(rbind, out)
}
