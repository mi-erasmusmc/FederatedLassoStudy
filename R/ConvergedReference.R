#' Converged reference models, and the saved models scored against them
#'
#' Refits, for every main-analysis evaluation, the Cyclops lasso fits of the
#' main run to convergence, on the same training data and at the same selected
#' penalties: pooled lasso; the lasso of the pooled training data at DualAvg's
#' penalty, which is the solution DualAvg approximates; and the local lasso of
#' each training database. Like [runAgeAnalysis()], the refits centre age and
#' replace binary predictors present in most entries by their absence, which
#' changes neither the solution nor the predictions but removes the
#' near-collinearity with the intercept along which the main run's solver
#' stopped early, and use a tight tolerance. Converged coefficient averaging,
#' probability ensembling and largest-database models are built from the
#' converged local fits with the saved weights.
#'
#' With `target`, the lasso developed in the held-out database (the
#' `TargetLasso` of [runProductionChecks()]: five outcome-stratified folds, each
#' tuned as a local lasso, scored out of fold) is fitted again, and each fold's
#' final fit is also refitted to convergence at the penalty its tuning
#' selected.
#'
#' Every saved model of the evaluation and every converged model is then scored
#' in the held-out database against the converged pooled lasso, with the
#' measures of [runProductionChecks()]: AUC with DeLong's 95% interval, the
#' paired DeLong difference from the reference, the model-based AUC, the
#' correlation of the linear predictors, the mean and 99th percentile absolute
#' difference in predicted risk and the overlap of the top 10% and top 1% of
#' predicted risk, plus calibration-in-the-large, slope and the framework's
#' calibration curve and integrated calibration index.
#'
#' The `export` directory holds:
#' - `referenceFits`: one row per saved and converged fit (`savedPooled`,
#'   `convergedPooled`, `savedDualAvg`, `convergedDualAvgPenalty`, `savedLocal`
#'   and `convergedLocal` per training database, and `tunedTarget` and
#'   `convergedTarget` per fold of the target-database lasso): its penalty on
#'   the mean log loss scale, lasso objective and largest optimality violation
#'   on its own training data, age coefficient, intercept, nonzero count and
#'   solver status. An `averageCheck` row per evaluation holds in `check` the
#'   largest difference between the saved coefficient average and the weighted
#'   sum of the saved local models, which must be about zero for the converged
#'   average to be built the same way. A failed target-database lasso leaves a
#'   `targetLasso` row with the error.
#' - `referenceScores`: one row per scored model (`kind` saved, converged or
#'   target), the converged pooled lasso first as the reference;
#'   `convergedLocal` rows name their training database in `sourceDatabase`.
#' - `referenceCurves`: every scored model's calibration curve, 50 points from
#'   the 1st to the 99th percentile of predicted risk with pointwise 95%
#'   intervals.
#' - `referenceCoefficients`: the converged models' coefficients by covariate
#'   id, age and the intercept included.
#'
#' Only aggregate tables are written; patient-level data stay in memory. A
#' rerun keeps the task and predictor-set pairs an earlier run finished; delete
#' the export directory to start again.
#'
#' @param dataRoot,siteIds,taskDirectories As for [runStudy()].
#' @param mainDirectories Result directories of the main analysis.
#' @param outputDirectory Directory for the `export` tables.
#' @param tasks,featureSets,holdouts The evaluations to refit.
#' @param tolerance,maxIterations,threads Cyclops controls of the refits.
#' @param target Whether to fit the target-database lasso.
#' @return The export tables, invisibly.
#' @export
runConvergedReference <- function(dataRoot, siteIds, mainDirectories, outputDirectory, tasks = names(studyTasks()),
                                  featureSets = c("ageSexPhenotypes", "ageSexConditions"), holdouts = siteIds,
                                  taskDirectories = NULL, tolerance = 1e-10, maxIterations = 20000L,
                                  threads = 1L, target = TRUE) {
  control <- list(tolerance = tolerance, maxIterations = maxIterations, threads = threads)
  tables <- list(referenceFits = list(), referenceScores = list(), referenceCurves = list(),
    referenceCoefficients = list())
  add <- function(table, rows) if (!is.null(rows) && nrow(rows)) tables[[table]][[length(tables[[table]]) + 1L]] <<- rows
  done <- character()
  files <- file.path(outputDirectory, "export", paste0(names(tables), ".csv"))
  if (all(file.exists(files))) {
    for (k in seq_along(tables)) add(names(tables)[[k]], utils::read.csv(files[[k]], stringsAsFactors = FALSE))
    done <- unique(paste(tables$referenceFits[[1]]$featureSet, tables$referenceFits[[1]]$task))
    message("Keeping ", length(done), " finished task and predictor-set pairs from ", outputDirectory)
  }
  available <- scoredModels(mainDirectories)
  for (featureSet in featureSets) for (task in tasks) {
    if (paste(featureSet, task) %in% done) next
    sites <- studySites(dataRoot, siteIds, task, featureSet, taskDirectories)[[task]]
    population <- studyTasks()[[task]]$populationSettings
    plpData <- list()
    for (holdout in holdouts) {
      message(format(Sys.time(), "%H:%M"), "  ", task, ", ", featureSet, ", ", holdout, " held out")
      rows <- available[available$task == task & available$featureSet == featureSet &
        available$heldOutDatabase == holdout, , drop = FALSE]
      if (!"PooledLasso" %in% rows$method) {
        stop("No main-analysis PooledLasso for ", task, " with ", featureSet, " and ", holdout, " held out",
          call. = FALSE)
      }
      rows <- rows[!duplicated(rows$method), , drop = FALSE]
      saved <- stats::setNames(lapply(seq_len(nrow(rows)), function(i) {
        FederatedLearning::loadModel(rows$directory[[i]], task, holdout, featureSet, rows$method[[i]])
      }), rows$method)
      matrices <- lapply(stats::setNames(siteIds, siteIds), function(site) {
        if (is.null(plpData[[site]])) plpData[[site]] <<- FederatedLearning::loadClientData(sites[[site]], population)
        config <- saved$PooledLasso$config
        # Held-out entries are all scored, as FederatedLearning::predictModel() scores them.
        if (site == holdout) config$diagnosticControlsPerCase <- NULL
        data <- FederatedLearning::createClientMatrix(plpData[[site]], config)
        list(xMatrix = data$xMatrix, yLabels = data$yLabels)
      })
      evaluation <- referenceEvaluation(matrices[setdiff(siteIds, holdout)], matrices[[holdout]], saved, control,
        fold = rows$fold[rows$method == "PooledLasso"][[1]], target = target)
      key <- data.frame(task = task, featureSet = featureSet, heldOutDatabase = holdout, stringsAsFactors = FALSE)
      for (table in names(evaluation)) {
        if (!is.null(evaluation[[table]]) && nrow(evaluation[[table]])) add(table, cbind(key, evaluation[[table]]))
      }
      rm(matrices, evaluation, saved)
      gc()
    }
    rm(plpData)
    gc()
    writeExport(outputDirectory, lapply(tables, function(rows) if (length(rows)) do.call(rbind, rows) else data.frame()))
  }
  invisible(lapply(tables, function(rows) if (length(rows)) do.call(rbind, rows) else data.frame()))
}

# The converged refits and held-out scores of one evaluation. `saved` holds the loaded main-analysis
# models by method; `train` the training databases' matrices by name; `fold` seeds the target lasso.
referenceEvaluation <- function(train, test, saved, control, fold = 1L, target = TRUE) {
  pooled <- saved$PooledLasso
  columns <- names(componentCoefficients(pooled$models[[1]]))
  if (length(columns) != ncol(test$xMatrix) || columns[[1]] != "(Intercept)") {
    stop("Pooled lasso's coefficients must match the rebuilt matrix, intercept first", call. = FALSE)
  }
  for (method in names(saved)) for (component in c(saved[[method]]$models, saved[[method]]$components)) {
    if (!identical(names(componentCoefficients(component)), columns)) {
      stop(method, " does not use pooled lasso's predictors", call. = FALSE)
    }
  }
  age <- match(as.character(ageCovariateId), columns)
  if (is.na(age)) stop("The predictors do not include age", call. = FALSE)
  sparse <- function(m) methods::as(methods::as(methods::as(m, "dMatrix"), "generalMatrix"), "CsparseMatrix")
  x <- sparse(do.call(rbind, lapply(train, `[[`, "xMatrix")))
  y <- unlist(lapply(train, `[[`, "yLabels"), use.names = FALSE)
  local <- lapply(train, function(t) list(x = sparse(t$xMatrix), y = t$yLabels))
  meanScale <- function(variance, n) sqrt(2 / variance) / n

  fits <- list()
  record <- function(fit, source, design, outcome, w, lambda, solver = NULL, flipped = NA_integer_,
                     error = NA_character_, check = NA_real_) {
    lasso <- if (!is.null(w) && is.finite(lambda)) lassoDiagnostics(design, outcome, w, lambda) else NULL
    fits[[length(fits) + 1L]] <<- data.frame(fit = fit, sourceDatabase = source,
      n = if (is.null(outcome)) NA_integer_ else length(outcome), outcomes = if (is.null(outcome)) NA_integer_ else sum(outcome),
      lambda = lambda, objective = if (is.null(lasso)) NA_real_ else lasso$objective,
      kktMax = if (is.null(lasso)) NA_real_ else lasso$kktMax,
      kktMaxCoordinate = if (is.null(lasso)) NA_character_ else columns[[lasso$kktMaxCoordinate]],
      ageCoefficient = if (is.null(w)) NA_real_ else w[[age]], intercept = if (is.null(w)) NA_real_ else w[[1]],
      nonzero = if (is.null(w)) NA_integer_ else sum(w[-1] != 0), flippedPredictors = flipped,
      status = if (is.null(solver)) NA_character_ else solver$status,
      iterations = if (is.null(solver)) NA_integer_ else solver$iterations,
      seconds = if (is.null(solver)) NA_real_ else solver$seconds, check = check, error = error,
      stringsAsFactors = FALSE)
    w
  }
  # A lasso at penalty `lambda` (mean log loss scale, intercept unpenalized) solved to convergence on
  # a design with age centred at its mean and binary predictors present in most entries flipped, its
  # coefficients returned on the original columns; NULL, with the error recorded, if it fails.
  converged <- function(fit, source, design, outcome, lambda) {
    attempt <- tryCatch({
      means <- Matrix::colMeans(design)
      center <- means[[age]]
      flipped <- which(vapply(seq_along(columns), function(j) {
        j > 1L && j != age && means[[j]] > 0.5 && binaryColumn(design, j)
      }, logical(1)))
      solver <- cyclopsFit(reparametrize(design, age, center, flipped), outcome, "laplace", 1L, control, lambda)
      list(w = stats::setNames(restoreCoefficients(solver$w, age, center, flipped), columns), solver = solver,
        flipped = length(flipped))
    }, error = function(e) conditionMessage(e))
    if (is.character(attempt)) {
      return(record(fit, source, design, outcome, NULL, lambda, error = substr(attempt, 1, 300)))
    }
    record(fit, source, design, outcome, attempt$w, lambda, attempt$solver, attempt$flipped)
  }

  savedPooled <- componentCoefficients(pooled$models[[1]])
  pooledLambda <- meanScale(pooled$models[[1]]$selectedLambda, length(y))
  record("savedPooled", NA_character_, x, y, savedPooled, pooledLambda)
  reference <- converged("convergedPooled", NA_character_, x, y, pooledLambda)
  convergedModels <- list()
  if (!is.null(saved$DualAvg)) {
    savedDualAvg <- componentCoefficients(saved$DualAvg$models[[1]])
    # DualAvg's saved penalty is already on the mean log loss scale.
    dualAvgLambda <- saved$DualAvg$models[[1]]$selectedLambda
    record("savedDualAvg", NA_character_, x, y, savedDualAvg, dualAvgLambda)
    convergedModels$convergedDualAvgPenalty <- converged("convergedDualAvgPenalty", NA_character_, x, y, dualAvgLambda)
  }

  # Each distinct local fit (database and selected penalty) is refitted once.
  localFits <- list()
  localFit <- function(component) {
    site <- as.character(component$sourceClientId)
    if (!site %in% names(local)) stop("Local model from ", site, ", which is not a training database", call. = FALSE)
    id <- paste(site, format(component$selectedLambda, digits = 15))
    if (is.null(localFits[[id]])) {
      d <- local[[site]]
      lambda <- meanScale(component$selectedLambda, length(d$y))
      record("savedLocal", site, d$x, d$y, componentCoefficients(component), lambda)
      localFits[[id]] <<- list(site = site, w = converged("convergedLocal", site, d$x, d$y, lambda))
    }
    localFits[[id]]$w
  }
  average <- saved$LocalAvgLasso
  if (!is.null(average) && length(average$components)) {
    weights <- vapply(average$components, `[[`, numeric(1), "weight")
    combined <- Reduce(`+`, Map(function(component, weight) weight * componentCoefficients(component),
      average$components, weights))
    record("averageCheck", NA_character_, NULL, NULL, NULL, NA_real_,
      check = max(abs(combined - componentCoefficients(average$models[[1]]))))
    parts <- lapply(average$components, localFit)
    if (!any(vapply(parts, is.null, logical(1)))) {
      convergedModels$convergedLocalAvgLasso <- Reduce(`+`, Map(`*`, weights, parts))
    }
  }
  ensemble <- NULL
  if (!is.null(saved$LocalEnsembleLasso)) {
    parts <- lapply(saved$LocalEnsembleLasso$models, localFit)
    if (!any(vapply(parts, is.null, logical(1)))) {
      ensemble <- list(w = parts, weights = vapply(saved$LocalEnsembleLasso$models, `[[`, numeric(1), "weight"))
    }
  }
  if (!is.null(saved$BiggestSiteLasso)) {
    convergedModels$convergedBiggestSiteLasso <- localFit(saved$BiggestSiteLasso$models[[1]])
  }

  # The lasso developed in the held-out database, as runProductionChecks() fits it, and the same with
  # each fold's final fit refitted to convergence at the penalty its tuning selected.
  targetPredictions <- NULL
  if (target) {
    targetPredictions <- tryCatch({
      framework <- asNamespace("FederatedLearning")
      site <- list(xMatrix = test$xMatrix, yLabels = test$yLabels)
      design <- sparse(test$xMatrix)
      seed <- 42L + fold
      outer <- framework$.stratifiedFoldId(site$yLabels, 5L, seed)
      tuned <- refitted <- numeric(length(site$yLabels))
      for (k in seq_len(5L)) {
        rows <- which(outer != k)
        held <- which(outer == k)
        fit <- framework$tuneLocalCyclopsFit(framework$clientRows(site, rows), pooled$args, seed + k,
          isTRUE(pooled$config$intercept))
        lambda <- meanScale(fit$selectedLambda, length(rows))
        xr <- design[rows, , drop = FALSE]
        w <- record("tunedTarget", paste("fold", k), xr, site$yLabels[rows], stats::setNames(as.numeric(fit$w), columns),
          lambda)
        v <- converged("convergedTarget", paste("fold", k), xr, site$yLabels[rows], lambda)
        if (is.null(v)) stop("The converged refit of fold ", k, " failed", call. = FALSE)
        tuned[held] <- as.numeric(design[held, , drop = FALSE] %*% w)
        refitted[held] <- as.numeric(design[held, , drop = FALSE] %*% v)
      }
      list(tuned = tuned, converged = refitted)
    }, error = function(e) {
      record("targetLasso", NA_character_, NULL, NULL, NULL, NA_real_, error = substr(conditionMessage(e), 1, 300))
      NULL
    })
  }

  coefficients <- do.call(rbind, c(
    list(if (!is.null(reference)) data.frame(model = "convergedPooled", sourceDatabase = NA_character_,
      covariateId = columns, coefficient = unname(reference), stringsAsFactors = FALSE)),
    lapply(names(convergedModels), function(name) {
      w <- convergedModels[[name]]
      if (is.null(w) || name == "convergedBiggestSiteLasso") return(NULL)
      data.frame(model = name, sourceDatabase = NA_character_, covariateId = columns, coefficient = unname(w),
        stringsAsFactors = FALSE)
    }),
    lapply(localFits, function(f) {
      if (is.null(f$w)) return(NULL)
      data.frame(model = "convergedLocal", sourceDatabase = f$site, covariateId = columns, coefficient = unname(f$w),
        stringsAsFactors = FALSE)
    })))
  rownames(coefficients) <- NULL
  if (is.null(reference)) {
    return(list(referenceFits = do.call(rbind, fits), referenceScores = NULL, referenceCurves = NULL,
      referenceCoefficients = coefficients))
  }

  # Held-out scores: the converged pooled lasso first, then every saved and converged model.
  candidate <- function(model, kind, linearPredictor, probability = stats::plogis(linearPredictor), source = NA_character_) {
    list(model = model, kind = kind, source = source, probability = probability, linearPredictor = linearPredictor)
  }
  lp <- function(w) as.numeric(test$xMatrix %*% w)
  candidates <- list(candidate("convergedPooled", "converged", lp(reference)))
  for (method in names(saved)) {
    predictions <- modelPredictions(saved[[method]], test)
    candidates[[length(candidates) + 1L]] <- candidate(method, "saved", predictions$linearPredictor,
      predictions$probability)
  }
  for (name in names(convergedModels)) {
    if (!is.null(convergedModels[[name]])) {
      candidates[[length(candidates) + 1L]] <- candidate(name, "converged", lp(convergedModels[[name]]))
    }
  }
  if (!is.null(ensemble)) {
    etas <- lapply(ensemble$w, lp)
    probability <- as.numeric(do.call(cbind, lapply(etas, stats::plogis)) %*% ensemble$weights)
    candidates[[length(candidates) + 1L]] <- candidate("convergedLocalEnsembleLasso", "converged",
      stats::qlogis(pmin(pmax(probability, 1e-15), 1 - 1e-15)), probability)
  }
  for (f in localFits) {
    if (!is.null(f$w)) candidates[[length(candidates) + 1L]] <- candidate("convergedLocal", "converged", lp(f$w), source = f$site)
  }
  if (!is.null(targetPredictions)) {
    candidates[[length(candidates) + 1L]] <- candidate("TargetLasso", "target", targetPredictions$tuned)
    candidates[[length(candidates) + 1L]] <- candidate("convergedTargetLasso", "target", targetPredictions$converged)
  }
  scored <- scoreAgainstReference(test$yLabels, candidates)
  list(referenceFits = do.call(rbind, fits), referenceScores = scored$scores, referenceCurves = scored$curves,
    referenceCoefficients = coefficients)
}

# Discrimination, calibration and agreement of each candidate (model, kind, source, probability,
# linear predictor) with the first, which is the reference; and each candidate's calibration curve.
scoreAgainstReference <- function(outcome, candidates) {
  framework <- asNamespace("FederatedLearning")
  delong <- delongAuc(outcome, do.call(cbind, lapply(candidates, `[[`, "probability")))
  z <- stats::qnorm(0.975)
  r <- candidates[[1]]
  curves <- list()
  scores <- do.call(rbind, lapply(seq_along(candidates), function(k) {
    m <- candidates[[k]]
    se <- sqrt(delong$covariance[k, k])
    differenceSe <- sqrt(max(delong$covariance[k, k] + delong$covariance[1, 1] - 2 * delong$covariance[k, 1], 0))
    difference <- delong$auc[k] - delong$auc[1]
    absolute <- abs(m$probability - r$probability)
    overlap <- function(share) {
      length(intersect(topRisk(m$probability, share), topRisk(r$probability, share))) / length(topRisk(r$probability, share))
    }
    calibration <- framework$.calibrationMetrics(outcome, m$linearPredictor)
    curve <- framework$.calibrationCurve(outcome, m$linearPredictor)
    if (!is.null(curve$curve)) {
      curves[[length(curves) + 1L]] <<- data.frame(model = m$model, kind = m$kind, sourceDatabase = m$source, curve$curve,
        stringsAsFactors = FALSE)
    }
    data.frame(model = m$model, kind = m$kind, sourceDatabase = m$source, n = length(outcome), outcomes = sum(outcome),
      auc = delong$auc[k], aucSe = se, aucLower = delong$auc[k] - z * se, aucUpper = delong$auc[k] + z * se,
      aucDifference = difference, differenceSe = differenceSe, differenceLower = difference - z * differenceSe,
      differenceUpper = difference + z * differenceSe, modelBasedAuc = modelBasedAuc(m$probability),
      calibrationInTheLarge = calibration[["calibrationInTheLarge"]], calibrationSlope = calibration[["calibrationSlope"]],
      calibrationIci = curve$ici,
      linearPredictorPearson = stats::cor(m$linearPredictor, r$linearPredictor),
      linearPredictorSpearman = stats::cor(m$linearPredictor, r$linearPredictor, method = "spearman"),
      probabilityMeanAbsDifference = mean(absolute),
      probabilityP99AbsDifference = unname(stats::quantile(absolute, 0.99)),
      top10Overlap = overlap(0.1), top1Overlap = overlap(0.01), stringsAsFactors = FALSE)
  }))
  list(scores = scores, curves = if (length(curves)) do.call(rbind, curves) else NULL)
}
