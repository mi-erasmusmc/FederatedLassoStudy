#' Checks of saved models on the site data, without refitting
#'
#' Answers four review questions from the saved models of one or more result
#' directories, reading each site's data once per task and extraction:
#' - `auc`: each model's held-out AUC with DeLong's 95% interval, and its
#'   paired DeLong difference from pooled lasso on the same patients. The
#'   reference is the pooled lasso of the same task, predictor set and fold in
#'   the first directory that has one, so list the main analysis first. The
#'   table also holds:
#'   - the model-based AUC (van Klaveren et al., 2016): the AUC expected in
#'     the held-out population if the model's own risks were correct, which
#'     reflects case mix alone;
#'   - the agreement of each model's linear predictor with pooled lasso's
#'     (Pearson and Spearman correlation) and the mean and 99th percentile
#'     absolute difference in predicted risk;
#'   - the share of pooled lasso's top 10% and top 1% of predicted risk that
#'     the model also places there, and the share of outcomes in its own top
#'     10%.
#'   For `targetFeatureSets`, a `TargetLasso` row adds a lasso fitted inside
#'   the held-out database by five outcome-stratified folds, each tuned as a
#'   local lasso, and scored on its out-of-fold predictions: the performance a
#'   model developed in that database itself could reach. If its penalty
#'   search fails in any fold, as a local lasso's can, the row keeps the
#'   failure in `error` and its measures are missing.
#' - `covariatePrevalence`: for each task and database, the share of the
#'   evaluated population with each of `prevalenceCovariateIds`, with counts
#'   below `minCellCount` (of patients with or without it) suppressed.
#' - `siteConstantPredictors`: predictors of the phenotype models that are
#'   constant within a database's population (absent from everyone or present
#'   in everyone).
#' - `odalInitialization`: the local unpenalized logistic fit that starts ODAL,
#'   replayed at every training database of each phenotype-predictor ODAL
#'   model as the fit ran it, recording whether it fell back to ridge and why,
#'   or the error that stopped it.
#'
#' Local site models are not scored. Only aggregate tables are written, to the
#' `export` directory; patient-level predictions stay in memory.
#'
#' @param resultsDirectories Result directories of [runStudy()], the main
#'   analysis first.
#' @param dataRoot,siteIds,taskDirectories As for [runStudy()].
#' @param outputDirectory Directory for the `export` tables.
#' @param prevalenceCovariateIds Covariates whose prevalence is reported.
#' @param minCellCount Smallest count reported.
#' @param targetFeatureSets Predictor sets for which the `TargetLasso`
#'   reference is fitted.
#' @return The export tables, invisibly.
#' @export
runProductionChecks <- function(resultsDirectories, dataRoot, siteIds, outputDirectory,
                                prevalenceCovariateIds = atrialFibrillationCovariateId, taskDirectories = NULL,
                                minCellCount = 5L,
                                targetFeatureSets = c("ageSexPhenotypes", "ageSexConditions")) {
  models <- scoredModels(resultsDirectories)
  models$extraction <- unname(extractionSuffixes[models$featureSet])
  odal <- models[models$method %in% c("ODAL1", "ODAL2") & models$extraction == "", , drop = FALSE]
  odal <- odal[!duplicated(odal[, c("task", "featureSet", "fold")]), , drop = FALSE]
  tables <- list(auc = list(), covariatePrevalence = list(), siteConstantPredictors = list(),
    odalInitialization = list())
  add <- function(table, rows) tables[[table]][[length(tables[[table]]) + 1L]] <<- rows
  for (task in unique(models$task)) for (extraction in unique(models$extraction[models$task == task])) {
    featureSet <- names(extractionSuffixes)[extractionSuffixes == extraction][[1]]
    sites <- studySites(dataRoot, siteIds, task, featureSet, taskDirectories)[[task]]
    for (site in siteIds) {
      message(format(Sys.time(), "%H:%M"), "  ", task, " ", if (nzchar(extraction)) extraction else "phenotypes",
        " ", site)
      plpData <- FederatedLearning::loadClientData(sites[[site]], studyTasks()[[task]]$populationSettings)
      matrices <- siteMatrices(plpData)
      heldOut <- models[models$task == task & models$extraction == extraction & models$heldOutDatabase == site, ,
        drop = FALSE]
      if (nrow(heldOut)) add("auc", scoreHeldOutModels(heldOut, matrices, targetFeatureSets))
      if (!nzchar(extraction)) {
        for (i in which(odal$task == task & odal$heldOutDatabase != site)) {
          add("odalInitialization", replayOdalInitialization(odal[i, ], matrices, site))
        }
        counts <- matrices$counts()
        add("covariatePrevalence", covariatePrevalence(counts, prevalenceCovariateIds, task, site, minCellCount))
        constant <- setdiff(names(counts$present)[counts$present %in% c(0, counts$n)], c("(Intercept)", "1002"))
        if (length(constant)) {
          add("siteConstantPredictors", data.frame(task = task, database = site, covariateId = constant,
            value = ifelse(counts$present[constant] == 0, "absent in all", "present in all")))
        }
      }
      rm(plpData, matrices)
      gc()
    }
  }
  tables <- lapply(tables, function(rows) if (length(rows)) do.call(rbind, rows) else data.frame())
  writeExport(outputDirectory, tables)
  invisible(tables)
}

# Every saved model of the directories except local site models, with its location.
scoredModels <- function(resultsDirectories) {
  do.call(rbind, lapply(resultsDirectories, function(directory) {
    results <- FederatedLearning::readComparisonResults(directory)
    ok <- !is.na(results$modelFile) & nzchar(results$modelFile) & results$method != "LocalSiteLasso"
    if ("error" %in% names(results)) ok <- ok & (is.na(results$error) | !nzchar(results$error))
    models <- unique(results[ok, c("task", "featureSet", "fold", "clientId", "method"), drop = FALSE])
    names(models)[names(models) == "clientId"] <- "heldOutDatabase"
    if (nrow(models)) cbind(directory = directory, models, stringsAsFactors = FALSE) else NULL
  }))
}

# One site's model matrices, built once per feature map, with the number of patients having each
# predictor in any matrix built so far.
siteMatrices <- function(plpData) {
  built <- list()
  present <- c()
  n <- NA_integer_
  list(
    get = function(model) {
      mapping <- model$config$mapping
      for (b in built) if (identical(b$mapping, mapping)) return(b$data)
      config <- model$config
      # Every patient, as FederatedLearning::predictModel() scores them.
      config$diagnosticControlsPerCase <- NULL
      data <- FederatedLearning::createClientMatrix(plpData, config)
      columns <- as.character(mapping$covariateId[order(mapping$columnId)])
      if (isTRUE(model$config$intercept)) columns <- c("(Intercept)", columns)
      counts <- Matrix::colSums(data$xMatrix != 0)
      present[columns] <<- counts
      n <<- nrow(data$xMatrix)
      built[[length(built) + 1L]] <<- list(mapping = mapping, data = data)
      data
    },
    counts = function() list(present = present, n = n))
}

# AUC with DeLong intervals and agreement with pooled lasso for the models held out at one site.
scoreHeldOutModels <- function(models, matrices, targetFeatureSets = character()) {
  scored <- lapply(seq_len(nrow(models)), function(i) {
    m <- models[i, ]
    model <- FederatedLearning::loadModel(m$directory, m$task, m$heldOutDatabase, m$featureSet, m$method)
    data <- matrices$get(model)
    out <- c(list(row = m), modelPredictions(model, data))
    if (m$method == "PooledLasso") out$target <- function() targetLassoPredictions(model, data, m$fold)
    out
  })
  groups <- split(seq_along(scored), paste(models$featureSet, models$fold))
  do.call(rbind, lapply(groups, function(index) {
    group <- scored[index]
    pooled <- which(vapply(group, function(s) s$row$method == "PooledLasso", logical(1)))[1]
    targetError <- NULL
    if (!is.na(pooled) && group[[pooled]]$row$featureSet %in% targetFeatureSets) {
      row <- group[[pooled]]$row
      row$method <- "TargetLasso"
      target <- tryCatch(group[[pooled]]$target(), error = function(e) conditionMessage(e))
      if (is.character(target)) {
        targetError <- list(row = row, error = target)
      } else {
        group[[length(group) + 1L]] <- c(list(row = row), target)
      }
    }
    reference <- which(vapply(group, function(s) s$row$method == "PooledLasso", logical(1)))[1]
    y <- group[[1]]$outcome
    # Probabilities, as the framework's AUC ranks them.
    delong <- delongAuc(y, do.call(cbind, lapply(group, `[[`, "probability")))
    z <- stats::qnorm(0.975)
    rows <- do.call(rbind, lapply(seq_along(group), function(k) {
      row <- group[[k]]$row
      se <- sqrt(delong$covariance[k, k])
      out <- data.frame(resultsDirectory = basename(row$directory), task = row$task, featureSet = row$featureSet,
        fold = row$fold, heldOutDatabase = row$heldOutDatabase, method = row$method, n = length(y),
        outcomes = sum(y), auc = delong$auc[k], aucSe = se, aucLower = delong$auc[k] - z * se,
        aucUpper = delong$auc[k] + z * se, referenceDirectory = NA_character_, aucDifference = NA_real_,
        differenceSe = NA_real_, differenceLower = NA_real_, differenceUpper = NA_real_, differenceP = NA_real_,
        modelBasedAuc = modelBasedAuc(group[[k]]$probability),
        top10OutcomeShare = sum(y[topRisk(group[[k]]$probability, 0.1)]) / sum(y),
        linearPredictorPearson = NA_real_, linearPredictorSpearman = NA_real_,
        probabilityMeanAbsDifference = NA_real_, probabilityP99AbsDifference = NA_real_,
        top10Overlap = NA_real_, top1Overlap = NA_real_, error = NA_character_, stringsAsFactors = FALSE)
      if (!is.na(reference)) {
        r <- group[[reference]]
        differenceSe <- sqrt(max(delong$covariance[k, k] + delong$covariance[reference, reference] -
          2 * delong$covariance[k, reference], 0))
        difference <- delong$auc[k] - delong$auc[reference]
        absolute <- abs(group[[k]]$probability - r$probability)
        out$referenceDirectory <- basename(r$row$directory)
        out$aucDifference <- difference
        out$differenceSe <- differenceSe
        out$differenceLower <- difference - z * differenceSe
        out$differenceUpper <- difference + z * differenceSe
        out$differenceP <- if (differenceSe > 0) 2 * stats::pnorm(-abs(difference) / differenceSe) else NA_real_
        out$linearPredictorPearson <- stats::cor(group[[k]]$linearPredictor, r$linearPredictor)
        out$linearPredictorSpearman <- stats::cor(group[[k]]$linearPredictor, r$linearPredictor, method = "spearman")
        out$probabilityMeanAbsDifference <- mean(absolute)
        out$probabilityP99AbsDifference <- unname(stats::quantile(absolute, 0.99))
        overlap <- function(share) {
          length(intersect(topRisk(group[[k]]$probability, share), topRisk(r$probability, share))) /
            length(topRisk(r$probability, share))
        }
        out$top10Overlap <- overlap(0.1)
        out$top1Overlap <- overlap(0.01)
      }
      out
    }))
    if (!is.null(targetError)) {
      failed <- rows[1, ]
      failed[] <- NA
      failed[, c("resultsDirectory", "task", "featureSet", "fold", "heldOutDatabase", "method", "n", "outcomes",
        "error")] <- list(basename(targetError$row$directory), targetError$row$task, targetError$row$featureSet,
        targetError$row$fold, targetError$row$heldOutDatabase, "TargetLasso", length(y), sum(y),
        substr(targetError$error, 1, 300))
      rows <- rbind(rows, failed)
    }
    rows
  }))
}

# Out-of-fold predictions of a lasso developed in the held-out database itself, with pooled lasso's
# retained predictors: five outcome-stratified folds, each fitted and tuned as the framework tunes a
# local lasso (three inner folds, the same penalty search and solver).
targetLassoPredictions <- function(pooled, data, fold, folds = 5L) {
  framework <- asNamespace("FederatedLearning")
  intercept <- isTRUE(pooled$config$intercept)
  # The saved feature map already holds only the retained predictors, unscaled.
  if (any(pooled$preprocessing$normFactors != 1)) stop("TargetLasso does not support normalized predictors", call. = FALSE)
  site <- list(xMatrix = data$xMatrix, yLabels = data$yLabels)
  seed <- 42L + fold
  outer <- framework$.stratifiedFoldId(site$yLabels, folds, seed)
  linearPredictor <- numeric(length(site$yLabels))
  for (k in seq_len(folds)) {
    fit <- framework$tuneLocalCyclopsFit(framework$clientRows(site, which(outer != k)), pooled$args, seed + k,
      intercept)
    held <- which(outer == k)
    linearPredictor[held] <- as.numeric(site$xMatrix[held, , drop = FALSE] %*% fit$w)
  }
  list(outcome = site$yLabels, probability = stats::plogis(linearPredictor), linearPredictor = linearPredictor)
}

# The entries with the highest `share` of predicted risk, ties broken by position.
topRisk <- function(probability, share) {
  order(-probability, seq_along(probability))[seq_len(max(1L, ceiling(share * length(probability))))]
}

#' Model-based AUC
#'
#' The AUC expected if the predicted risks were the true risks (van Klaveren
#' et al., 2016): over all ordered pairs of distinct entries, the share of the
#' expected outcome-discordant pairs that the risks order correctly, with ties
#' counted as one half. It depends on the spread of the risks, that is, on case
#' mix, and not on the observed outcomes.
#'
#' @param probability Predicted risks.
#' @return The model-based AUC.
#' @keywords internal
modelBasedAuc <- function(probability) {
  p <- sort(probability)
  tie <- cumsum(c(TRUE, diff(p) != 0))
  risk <- rowsum(p, tie)[, 1]
  noRisk <- rowsum(1 - p, tie)[, 1]
  both <- rowsum(p * (1 - p), tie)[, 1]
  concordant <- sum(risk * (cumsum(noRisk) - noRisk)) + 0.5 * sum(risk * noRisk - both)
  concordant / (sum(p) * sum(1 - p) - sum(p * (1 - p)))
}

# Predicted risk and linear predictor of a saved model, as FederatedLearning::predictModel() computes
# them; a probability ensemble's linear predictor is the logit of its averaged risk.
modelPredictions <- function(model, data) {
  etas <- lapply(model$models, function(component) {
    coefficients <- component$coefficients[order(component$coefficients$matrixColumn), , drop = FALSE]
    if (nrow(coefficients) != ncol(data$xMatrix)) stop("Saved coefficients do not match the rebuilt matrix", call. = FALSE)
    as.numeric(data$xMatrix %*% coefficients$coefficientSharedScale)
  })
  weights <- if (length(etas) > 1L) vapply(model$models, `[[`, numeric(1), "weight") else 1
  # Summed as the framework's evaluation sums them, so that ties and the AUC match the reported ones.
  probability <- as.numeric(do.call(cbind, lapply(etas, stats::plogis)) %*% weights)
  linearPredictor <- if (length(etas) == 1L) etas[[1]] else stats::qlogis(pmin(pmax(probability, 1e-15), 1 - 1e-15))
  list(outcome = data$yLabels, probability = probability, linearPredictor = linearPredictor)
}

#' DeLong AUCs and their covariance
#'
#' The fast algorithm of Sun and Xu (2014) for DeLong's covariance of the AUCs
#' of several scores on the same patients, with ties counted as one half.
#'
#' @param y Binary outcomes.
#' @param scores A matrix with one column of scores per model.
#' @return A list with `auc` (one per column) and their `covariance` matrix.
#' @keywords internal
delongAuc <- function(y, scores) {
  scores <- as.matrix(scores)
  positive <- y == 1
  m <- sum(positive)
  n <- sum(!positive)
  if (m < 2L || n < 2L) stop("DeLong needs at least two patients with and without the outcome", call. = FALSE)
  components <- lapply(seq_len(ncol(scores)), function(k) {
    x <- scores[positive, k]
    z <- scores[!positive, k]
    combined <- rank(c(x, z))
    list(v10 = (combined[seq_len(m)] - rank(x)) / n, v01 = 1 - (combined[m + seq_len(n)] - rank(z)) / m)
  })
  v10 <- do.call(cbind, lapply(components, `[[`, "v10"))
  v01 <- do.call(cbind, lapply(components, `[[`, "v01"))
  list(auc = colMeans(v10), covariance = stats::cov(v10) / m + stats::cov(v01) / n)
}

# Share of a site's evaluated population with each covariate; small counts suppressed.
covariatePrevalence <- function(counts, covariateIds, task, site, minCellCount) {
  do.call(rbind, lapply(covariateIds, function(id) {
    with <- unname(counts$present[as.character(id)])
    small <- !is.na(with) && ((with > 0 && with < minCellCount) || (counts$n - with > 0 && counts$n - with < minCellCount))
    data.frame(task = task, database = site, covariateId = id, entries = counts$n,
      withCovariate = if (small) NA_real_ else with,
      percent = if (small || is.na(with)) NA_real_ else round(100 * with / counts$n, 2), suppressed = small)
  }))
}

# The local logistic fit with which ODAL starts at one training site, as the fit ran it.
replayOdalInitialization <- function(odalModel, matrices, site) {
  model <- FederatedLearning::loadModel(odalModel$directory, odalModel$task, odalModel$heldOutDatabase,
    odalModel$featureSet, odalModel$method)
  config <- model$config
  data <- matrices$get(model)
  framework <- asNamespace("FederatedLearning")
  xRaw <- framework$.stripInterceptColumn(data$xMatrix, config)
  # glm's warnings are recorded in the diagnostics, so they need not reach the console.
  init <- tryCatch(suppressWarnings(framework$.fitLocalLogisticDiagnosed(xRaw = xRaw,
    xDesign = framework$.addInterceptColumn(xRaw), y = data$yLabels, config = config,
    allowRidgeFallback = identical(framework$.odalInitMode(config), "ridgeFallback"))),
    error = function(e) list(beta = NA_real_, diagnostics = list(initMethod = NA_character_,
      fallbackReason = paste("Error:", conditionMessage(e)), glmWarnings = NA_character_)))
  d <- init$diagnostics
  text <- function(x) {
    x <- x[!is.na(x) & nzchar(x)]
    if (length(x)) paste(unique(x), collapse = "; ") else NA_character_
  }
  data.frame(resultsDirectory = basename(odalModel$directory), task = odalModel$task,
    featureSet = odalModel$featureSet, fold = odalModel$fold, heldOutDatabase = odalModel$heldOutDatabase,
    trainingDatabase = site, n = length(data$yLabels), outcomes = sum(data$yLabels), p = ncol(xRaw),
    initMethod = text(d$initMethod), usedRidgeFallback = isTRUE(d$usedRidgeFallback),
    fallbackReason = text(d$fallbackReason), glmWarnings = text(d$glmWarnings),
    betaMaxAbs = max(abs(init$beta)), stringsAsFactors = FALSE)
}
