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
#' Every saved model of the evaluation and every converged model is then scored
#' in the held-out database against the converged pooled lasso, with the
#' measures of [runProductionChecks()]: AUC with DeLong's 95% interval, the
#' paired DeLong difference from the reference, the model-based AUC, the
#' correlation of the linear predictors, the mean and 99th percentile absolute
#' difference in predicted risk and the overlap of the top 10% and top 1% of
#' predicted risk, plus calibration-in-the-large and slope.
#'
#' The `export` directory holds:
#' - `referenceFits`: one row per saved and converged fit (`savedPooled`,
#'   `convergedPooled`, `savedDualAvg`, `convergedDualAvgPenalty`, and
#'   `savedLocal` and `convergedLocal` per training database): its penalty on
#'   the mean log loss scale, lasso objective and largest optimality violation
#'   on its own training data, age coefficient, intercept, nonzero count and
#'   solver status. An `averageCheck` row per evaluation holds in `check` the
#'   largest difference between the saved coefficient average and the weighted
#'   sum of the saved local models, which must be about zero for the converged
#'   average to be built the same way.
#' - `referenceScores`: one row per scored model (`kind` saved or converged),
#'   the converged pooled lasso first as the reference.
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
#' @return The export tables, invisibly.
#' @export
runConvergedReference <- function(dataRoot, siteIds, mainDirectories, outputDirectory, tasks = names(studyTasks()),
                                  featureSets = c("ageSexPhenotypes", "ageSexConditions"), holdouts = siteIds,
                                  taskDirectories = NULL, tolerance = 1e-10, maxIterations = 20000L,
                                  threads = 1L) {
  control <- list(tolerance = tolerance, maxIterations = maxIterations, threads = threads)
  tables <- list(referenceFits = list(), referenceScores = list(), referenceCoefficients = list())
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
      evaluation <- referenceEvaluation(matrices[setdiff(siteIds, holdout)], matrices[[holdout]], saved, control)
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
# models by method; `train` the training databases' matrices by name.
referenceEvaluation <- function(train, test, saved, control) {
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
      localFits[[id]] <<- list(w = converged("convergedLocal", site, d$x, d$y, lambda))
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

  coefficients <- do.call(rbind, c(
    list(if (!is.null(reference)) data.frame(model = "convergedPooled", sourceDatabase = NA_character_,
      covariateId = columns, coefficient = unname(reference), stringsAsFactors = FALSE)),
    lapply(names(convergedModels), function(name) {
      w <- convergedModels[[name]]
      if (is.null(w) || name == "convergedBiggestSiteLasso") return(NULL)
      data.frame(model = name, sourceDatabase = NA_character_, covariateId = columns, coefficient = unname(w),
        stringsAsFactors = FALSE)
    }),
    lapply(names(localFits), function(id) {
      w <- localFits[[id]]$w
      if (is.null(w)) return(NULL)
      data.frame(model = "convergedLocal", sourceDatabase = sub(" .*$", "", id), covariateId = columns,
        coefficient = unname(w), stringsAsFactors = FALSE)
    })))
  if (is.null(reference)) {
    return(list(referenceFits = do.call(rbind, fits), referenceScores = NULL, referenceCoefficients = coefficients))
  }

  # Held-out scores: the converged pooled lasso first, then every saved and converged model.
  single <- function(w) {
    eta <- as.numeric(test$xMatrix %*% w)
    list(probability = stats::plogis(eta), linearPredictor = eta)
  }
  candidates <- list(convergedPooled = c(single(reference), kind = "converged"))
  for (method in names(saved)) {
    predictions <- modelPredictions(saved[[method]], test)
    candidates[[method]] <- list(probability = predictions$probability, linearPredictor = predictions$linearPredictor,
      kind = "saved")
  }
  for (name in names(convergedModels)) {
    if (!is.null(convergedModels[[name]])) candidates[[name]] <- c(single(convergedModels[[name]]), kind = "converged")
  }
  if (!is.null(ensemble)) {
    etas <- lapply(ensemble$w, function(w) as.numeric(test$xMatrix %*% w))
    probability <- as.numeric(do.call(cbind, lapply(etas, stats::plogis)) %*% ensemble$weights)
    candidates$convergedLocalEnsembleLasso <- list(probability = probability,
      linearPredictor = stats::qlogis(pmin(pmax(probability, 1e-15), 1 - 1e-15)), kind = "converged")
  }
  list(referenceFits = do.call(rbind, fits), referenceScores = scoreAgainstReference(test$yLabels, candidates),
    referenceCoefficients = coefficients)
}

# Discrimination, calibration and agreement of each candidate (probability, linear predictor, kind)
# with the first, which is the reference.
scoreAgainstReference <- function(outcome, candidates) {
  delong <- delongAuc(outcome, do.call(cbind, lapply(candidates, `[[`, "probability")))
  z <- stats::qnorm(0.975)
  r <- candidates[[1]]
  do.call(rbind, lapply(seq_along(candidates), function(k) {
    m <- candidates[[k]]
    se <- sqrt(delong$covariance[k, k])
    differenceSe <- sqrt(max(delong$covariance[k, k] + delong$covariance[1, 1] - 2 * delong$covariance[k, 1], 0))
    difference <- delong$auc[k] - delong$auc[1]
    lp <- m$linearPredictor
    absolute <- abs(m$probability - r$probability)
    overlap <- function(share) {
      length(intersect(topRisk(m$probability, share), topRisk(r$probability, share))) / length(topRisk(r$probability, share))
    }
    calibration <- suppressWarnings(c(
      stats::coef(stats::glm(outcome ~ 1 + offset(lp), family = stats::binomial()))[[1]],
      stats::coef(stats::glm(outcome ~ lp, family = stats::binomial()))[[2]]))
    data.frame(model = names(candidates)[[k]], kind = m$kind, n = length(outcome), outcomes = sum(outcome),
      auc = delong$auc[k], aucSe = se, aucLower = delong$auc[k] - z * se, aucUpper = delong$auc[k] + z * se,
      aucDifference = difference, differenceSe = differenceSe, differenceLower = difference - z * differenceSe,
      differenceUpper = difference + z * differenceSe, modelBasedAuc = modelBasedAuc(m$probability),
      calibrationInTheLarge = calibration[[1]], calibrationSlope = calibration[[2]],
      linearPredictorPearson = stats::cor(lp, r$linearPredictor),
      linearPredictorSpearman = stats::cor(lp, r$linearPredictor, method = "spearman"),
      probabilityMeanAbsDifference = mean(absolute),
      probabilityP99AbsDifference = unname(stats::quantile(absolute, 0.99)),
      top10Overlap = overlap(0.1), top1Overlap = overlap(0.01), stringsAsFactors = FALSE)
  }))
}
