#' Age effect: convergence, and within- and between-database components
#'
#' In bleeding after warfarin, workflows with nearly equal training loss
#' disagreed on the age coefficient (age in years / 100) by more than one unit,
#' even in sign. Age enters uncentred, so its coefficient and the intercept form
#' a direction along which the solvers converge slowly. For each task,
#' predictor set and held-out database of the main analysis, this analysis
#' refits on the pooled outer-training data with age centred at its pooled
#' training mean and binary predictors present in most entries replaced by
#' their absence. This reparametrization leaves every solution unchanged,
#' because the intercept is unpenalized and the penalties ignore a
#' coefficient's sign, but makes it fast to reach; coefficients are reported
#' for the original predictors. The fits are:
#' - `tightPooled`: lasso at pooled lasso's selected variance, solved to a
#'   tight tolerance;
#' - `tightDualAvgPenalty`: lasso at DualAvg's selected penalty, likewise;
#' - `logistic`: logistic regression with the intercept, age and the female
#'   indicator unpenalized and a weak normal prior (variance `priorVariance`)
#'   on the other predictors, which keeps rare predictors finite;
#' - `logisticWithin`: the same with unpenalized indicators of the training
#'   databases, so that age is estimated from differences within databases;
#' - `ageOnly` and `ageOnlyWithin`: age alone, without and with the database
#'   indicators;
#' - `logisticSite`: `logistic` in each training database alone.
#' The saved pooled lasso and DualAvg (`savedPooled`, `savedDualAvg`) are
#' evaluated on the same data.
#'
#' The `export` directory holds:
#' - `ageFits.csv`, one row per fit: the age coefficient with its standard
#'   error (from the observed information, including the prior's), the
#'   intercept on the uncentred scale, the training mean log loss, the solver's
#'   status, iterations and time, and for lasso fits the penalty, the objective
#'   and the largest violation of the optimality conditions on the mean log
#'   loss scale, as the framework measures them;
#' - `ageScores.csv`, held-out scores of the refits and of the saved pooled
#'   lasso and DualAvg, also with their age coefficient replaced by another
#'   fit's (`ageFrom`): AUC with the paired DeLong difference from the saved
#'   pooled lasso, calibration slope and the correlation of the linear
#'   predictor with saved pooled lasso's. The intercept affects none of these;
#' - `ageDatabases.csv`, each database's entries, outcomes and mean and
#'   standard deviation of age, overall and by outcome, with means over fewer
#'   than `minCellCount` entries suppressed.
#' The tables are rewritten after each task and predictor set. A rerun keeps
#' the pairs an earlier run wrote and fits only the others; delete the export
#' directory to start over.
#'
#' @param dataRoot,siteIds,taskDirectories As for [runStudy()].
#' @param mainDirectories Result directories of the main analysis.
#' @param outputDirectory Directory for the `export` tables.
#' @param tasks,featureSets,holdouts Evaluations to run, in this order (each
#'   predictor set over all tasks before the next).
#' @param priorVariance Variance of the weak normal prior of the logistic fits.
#' @param tolerance,maxIterations Cyclops convergence settings of every refit.
#' @param threads Cyclops threads.
#' @param minCellCount Smallest number of entries a reported mean may use.
#' @return The export tables, invisibly.
#' @export
runAgeAnalysis <- function(dataRoot, siteIds, mainDirectories, outputDirectory, tasks = names(studyTasks()),
                           featureSets = c("ageSexPhenotypes", "ageSexConditions"), holdouts = siteIds,
                           taskDirectories = NULL, priorVariance = 10, tolerance = 1e-10, maxIterations = 20000L,
                           threads = 1L, minCellCount = 5L) {
  control <- list(tolerance = tolerance, maxIterations = maxIterations, threads = threads)
  tables <- list(ageFits = list(), ageScores = list(), ageDatabases = list())
  add <- function(table, rows) tables[[table]][[length(tables[[table]]) + 1L]] <<- rows
  # A rerun keeps the task and predictor-set pairs an earlier run finished.
  done <- character()
  files <- file.path(outputDirectory, "export", paste0(names(tables), ".csv"))
  if (all(file.exists(files))) {
    for (k in seq_along(tables)) add(names(tables)[[k]], utils::read.csv(files[[k]], stringsAsFactors = FALSE))
    done <- unique(paste(tables$ageFits[[1]]$featureSet, tables$ageFits[[1]]$task))
    message("Keeping ", length(done), " finished task and predictor-set pairs from ", outputDirectory)
  }
  for (featureSet in featureSets) for (task in tasks) {
    if (paste(featureSet, task) %in% done) next
    sites <- studySites(dataRoot, siteIds, task, featureSet, taskDirectories)[[task]]
    population <- studyTasks()[[task]]$populationSettings
    plpData <- list()
    for (holdout in holdouts) {
      message(format(Sys.time(), "%H:%M"), "  ", task, ", ", featureSet, ", ", holdout, " held out")
      load <- function(method) {
        FederatedLearning::loadModel(mainDirectory(mainDirectories, task, featureSet, holdout, method), task,
          holdout, featureSet, method)
      }
      pooled <- load("PooledLasso")
      dualAvg <- load("DualAvg")
      matrices <- lapply(stats::setNames(siteIds, siteIds), function(site) {
        if (is.null(plpData[[site]])) plpData[[site]] <<- FederatedLearning::loadClientData(sites[[site]], population)
        config <- pooled$config
        # Held-out entries are all scored, as FederatedLearning::predictModel() scores them.
        if (site == holdout) config$diagnosticControlsPerCase <- NULL
        data <- FederatedLearning::createClientMatrix(plpData[[site]], config)
        list(xMatrix = data$xMatrix, yLabels = data$yLabels)
      })
      if (featureSet == featureSets[[1]] && holdout == holdouts[[1]]) {
        add("ageDatabases", databaseAges(matrices, pooled, task, minCellCount))
      }
      evaluation <- ageEvaluation(matrices[setdiff(siteIds, holdout)], matrices[[holdout]], pooled, dualAvg,
        priorVariance, control)
      key <- data.frame(task = task, featureSet = featureSet, heldOutDatabase = holdout, stringsAsFactors = FALSE)
      add("ageFits", cbind(key, evaluation$fits))
      add("ageScores", cbind(key, evaluation$scores))
      rm(matrices, evaluation)
      gc()
    }
    rm(plpData)
    gc()
    writeExport(outputDirectory, lapply(tables, function(rows) if (length(rows)) do.call(rbind, rows) else data.frame()))
  }
  invisible(lapply(tables, function(rows) if (length(rows)) do.call(rbind, rows) else data.frame()))
}


# The refits and held-out scores of one evaluation.
ageEvaluation <- function(train, test, pooled, dualAvg, priorVariance, control) {
  savedPooled <- componentCoefficients(pooled$models[[1]])
  savedDualAvg <- componentCoefficients(dualAvg$models[[1]])
  columns <- names(savedPooled)
  if (!identical(columns, names(savedDualAvg)) || length(columns) != ncol(test$xMatrix) ||
    columns[[1]] != "(Intercept)") {
    stop("Pooled lasso and DualAvg must share the rebuilt matrix's columns, intercept first", call. = FALSE)
  }
  age <- match(as.character(ageCovariateId), columns)
  if (is.na(age)) stop("The predictors do not include age", call. = FALSE)
  unpenalized <- c(1L, age, match(as.character(femaleCovariateId), columns))
  unpenalized <- unpenalized[!is.na(unpenalized)]
  x <- methods::as(methods::as(methods::as(do.call(rbind, lapply(train, `[[`, "xMatrix")), "dMatrix"),
    "generalMatrix"), "CsparseMatrix")
  y <- unlist(lapply(train, `[[`, "yLabels"), use.names = FALSE)
  n <- length(y)
  database <- rep(seq_along(train), vapply(train, function(t) length(t$yLabels), numeric(1)))
  indicators <- Matrix::sparseMatrix(i = which(database > 1L), j = database[database > 1L] - 1L, x = 1,
    dims = c(n, length(train) - 1L))
  # The refits use an exact reparametrization: age centred at its pooled training mean, and binary
  # predictors present in most entries replaced by their absence. Neither changes a lasso or normal-prior
  # solution, because the intercept is unpenalized and the penalties ignore a coefficient's sign, but
  # both remove near-collinearity with the intercept, along which the solvers converge slowly.
  means <- Matrix::colMeans(x)
  center <- means[[age]]
  flipped <- which(vapply(seq_along(columns), function(j) {
    j > 1L && j != age && means[[j]] > 0.5 && binaryColumn(x, j)
  }, logical(1)))
  centred <- reparametrize(x, age, center, flipped)
  pooledLambda <- sqrt(2 / pooled$models[[1]]$selectedLambda) / n
  dualAvgLambda <- dualAvg$models[[1]]$selectedLambda

  fits <- list()
  record <- function(name, w = NULL, loss = NA_real_, lambda = NA_real_, solver = NULL, se = NA_real_,
                     site = NA_character_, rows = seq_len(n), error = NA_character_) {
    lasso <- if (!is.null(w) && is.finite(lambda)) lassoDiagnostics(x, y, w, lambda) else NULL
    fits[[length(fits) + 1L]] <<- data.frame(fit = name, trainingDatabase = site, n = length(rows),
      outcomes = sum(y[rows]), ageCenter = center, flippedPredictors = length(flipped),
      ageCoefficient = if (is.null(w)) NA_real_ else w[[age]], ageSe = se, intercept = if (is.null(w)) NA_real_ else w[[1]],
      nonzero = if (is.null(w)) NA_integer_ else sum(w[-1] != 0), meanLogLoss = loss, lambda = lambda,
      objective = if (is.null(lasso)) NA_real_ else lasso$objective,
      kktMax = if (is.null(lasso)) NA_real_ else lasso$kktMax,
      kktMaxCoordinate = if (is.null(lasso)) NA_character_ else columns[[lasso$kktMaxCoordinate]],
      status = if (is.null(solver)) NA_character_ else solver$status,
      iterations = if (is.null(solver)) NA_integer_ else solver$iterations,
      seconds = if (is.null(solver)) NA_real_ else solver$seconds, error = error, stringsAsFactors = FALSE)
    w
  }
  # A refit on a reparametrized design whose first columns are `columns[index]` (intercept first),
  # returned as coefficients of the original `columns`; further design columns (database indicators)
  # are dropped, so a returned intercept is that of the first training database.
  refit <- function(name, design, index, prior, exclude = integer(), lambda = NA_real_, standardError = TRUE,
                    rows = seq_len(n), site = NA_character_) {
    position <- match(age, index)
    attempt <- tryCatch({
      solver <- cyclopsFit(design, y[rows], prior, exclude, control, lambda)
      precision <- if (is.list(prior)) replace(rep(1 / prior[[2]], ncol(design)), exclude, 0) else rep(0, ncol(design))
      se <- if (standardError) ageStandardError(design, y[rows], solver$w, precision, position) else NA_real_
      w <- stats::setNames(numeric(length(columns)), columns)
      w[index] <- solver$w[seq_along(index)]
      w <- restoreCoefficients(w, age, center, flipped)
      list(w = w, solver = solver, se = se, loss = meanLogLoss(design, y[rows], solver$w))
    }, error = function(e) conditionMessage(e))
    if (is.character(attempt)) return(record(name, site = site, rows = rows, error = substr(attempt, 1, 300)))
    record(name, attempt$w, attempt$loss, lambda, attempt$solver, attempt$se, site, rows)
  }

  all <- seq_along(columns)
  free <- c(unpenalized, length(columns) + seq_len(ncol(indicators)))
  weak <- list("normal", priorVariance)
  record("savedPooled", savedPooled, meanLogLoss(x, y, savedPooled), pooledLambda)
  record("savedDualAvg", savedDualAvg, meanLogLoss(x, y, savedDualAvg), dualAvgLambda)
  tight <- refit("tightPooled", centred, all, "laplace", 1L, pooledLambda, standardError = FALSE)
  tightDualAvg <- refit("tightDualAvgPenalty", centred, all, "laplace", 1L, dualAvgLambda, standardError = FALSE)
  logistic <- refit("logistic", centred, all, weak, unpenalized)
  logisticWithin <- refit("logisticWithin", cbind(centred, indicators), all, weak, free)
  refit("ageOnly", centred[, c(1L, age), drop = FALSE], c(1L, age), "none")
  refit("ageOnlyWithin", cbind(centred[, c(1L, age), drop = FALSE], indicators), c(1L, age), "none")
  for (k in seq_along(train)) {
    rows <- which(database == k)
    refit("logisticSite", centred[rows, , drop = FALSE], all, weak, unpenalized, rows = rows, site = names(train)[[k]])
  }

  # Held-out scores, the saved pooled lasso first as the reference; failed fits are skipped.
  ageOf <- function(w) if (is.null(w)) NA_real_ else w[[age]]
  withAge <- function(w, value) if (is.null(w) || is.na(value)) NULL else replace(w, age, value)
  candidates <- list(
    list("savedPooled", NA_character_, savedPooled), list("savedDualAvg", NA_character_, savedDualAvg),
    list("tightPooled", NA_character_, tight), list("tightDualAvgPenalty", NA_character_, tightDualAvg),
    list("logistic", NA_character_, logistic), list("logisticWithin", NA_character_, logisticWithin),
    list("savedPooled", "zero", withAge(savedPooled, 0)),
    list("savedPooled", "tightPooled", withAge(savedPooled, ageOf(tight))),
    list("savedPooled", "logistic", withAge(savedPooled, ageOf(logistic))),
    list("savedPooled", "logisticWithin", withAge(savedPooled, ageOf(logisticWithin))),
    list("savedPooled", "savedDualAvg", withAge(savedPooled, ageOf(savedDualAvg))),
    list("savedDualAvg", "savedPooled", withAge(savedDualAvg, ageOf(savedPooled))),
    list("savedDualAvg", "tightDualAvgPenalty", withAge(savedDualAvg, ageOf(tightDualAvg))))
  candidates <- Filter(function(candidate) !is.null(candidate[[3]]), candidates)
  linearPredictors <- lapply(candidates, function(candidate) as.numeric(test$xMatrix %*% candidate[[3]]))
  outcome <- test$yLabels
  delong <- delongAuc(outcome, do.call(cbind, lapply(linearPredictors, stats::plogis)))
  z <- stats::qnorm(0.975)
  scores <- do.call(rbind, lapply(seq_along(candidates), function(k) {
    se <- sqrt(max(delong$covariance[k, k] + delong$covariance[1, 1] - 2 * delong$covariance[k, 1], 0))
    difference <- delong$auc[k] - delong$auc[1]
    lp <- linearPredictors[[k]]
    slope <- suppressWarnings(stats::coef(stats::glm(outcome ~ lp, family = stats::binomial()))[[2]])
    data.frame(model = candidates[[k]][[1]], ageFrom = candidates[[k]][[2]], ageCoefficient = ageOf(candidates[[k]][[3]]),
      n = length(outcome), outcomes = sum(outcome), auc = delong$auc[k], aucDifference = difference,
      differenceLower = difference - z * se, differenceUpper = difference + z * se, calibrationSlope = slope,
      linearPredictorPearson = stats::cor(lp, linearPredictors[[1]]), stringsAsFactors = FALSE)
  }))
  list(fits = do.call(rbind, fits), scores = scores)
}

# A logistic regression fitted by Cyclops on `design` (intercept column first). `prior` is "none",
# "laplace" (penalty `lambda` on the mean log loss scale, as the framework converts Cyclops variances)
# or list("normal", variance); `exclude` lists the unpenalized design columns.
cyclopsFit <- function(design, y, prior, exclude, control, lambda = NA_real_) {
  data <- Cyclops::createCyclopsData(y = y, sx = design, modelType = "lr")
  ids <- as.numeric(Cyclops::getCovariateIds(data))
  cyclopsPrior <- if (identical(prior, "none")) {
    Cyclops::createPrior("none")
  } else if (identical(prior, "laplace")) {
    Cyclops::createPrior("laplace", variance = 2 / (lambda * length(y))^2, exclude = ids[exclude])
  } else {
    Cyclops::createPrior(prior[[1]], variance = prior[[2]], exclude = ids[exclude])
  }
  started <- Sys.time()
  fit <- Cyclops::fitCyclopsModel(data, prior = cyclopsPrior, control = Cyclops::createControl(
    tolerance = control$tolerance, maxIterations = control$maxIterations, noiseLevel = "silent",
    threads = control$threads, seed = 42L), warnings = FALSE)
  list(w = as.numeric(stats::coef(fit)), status = fit$return_flag, iterations = fit$iterations,
    seconds = as.numeric(difftime(Sys.time(), started, units = "secs")))
}

# Whether column j of a sparse matrix holds only zeros and ones.
binaryColumn <- function(x, j) {
  all(x@x[seq.int(x@p[j] + 1L, length.out = x@p[j + 1L] - x@p[j])] == 1)
}

# The refit design: column `age` minus `center`, and each `flipped` binary column replaced by one minus
# it. Built from the nonzero entries, because assigning a column of a large sparse matrix is slow.
reparametrize <- function(x, age, center, flipped) {
  entries <- Matrix::summary(x)
  keep <- !(entries$j %in% c(age, flipped))
  ageValues <- as.numeric(x[, age]) - center
  i <- list(entries$i[keep], which(ageValues != 0))
  j <- list(entries$j[keep], rep(age, length(i[[2]])))
  values <- list(entries$x[keep], ageValues[i[[2]]])
  for (k in flipped) {
    absent <- which(as.numeric(x[, k]) == 0)
    i[[length(i) + 1L]] <- absent
    j[[length(j) + 1L]] <- rep(k, length(absent))
    values[[length(values) + 1L]] <- rep(1, length(absent))
  }
  Matrix::sparseMatrix(i = unlist(i), j = unlist(j), x = unlist(values), dims = dim(x))
}

# Coefficients of the original columns from those of the reparametrized design.
restoreCoefficients <- function(w, age, center, flipped) {
  w[[1]] <- w[[1]] - w[[age]] * center + sum(w[flipped])
  w[flipped] <- -w[flipped]
  w
}

# Mean log loss of a logistic model.
meanLogLoss <- function(x, y, w) {
  eta <- as.numeric(x %*% w)
  mean(pmax(eta, 0) + log1p(exp(-abs(eta))) - y * eta)
}

# Lasso objective (mean log loss plus lambda times the L1 norm of all but the intercept) and the
# largest violation of its optimality conditions, as the framework measures them.
lassoDiagnostics <- function(x, y, w, lambda) {
  gradient <- as.numeric(Matrix::crossprod(x, stats::plogis(as.numeric(x %*% w)) - y)) / length(y)
  violation <- abs(gradient)
  penalized <- seq_along(w) > 1L
  active <- penalized & abs(w) > 1e-8
  violation[active] <- abs(gradient[active] + lambda * sign(w[active]))
  inactive <- penalized & !active
  violation[inactive] <- pmax(abs(gradient[inactive]) - lambda, 0)
  list(objective = meanLogLoss(x, y, w) + lambda * sum(abs(w[penalized])), kktMax = max(violation),
    kktMaxCoordinate = which.max(violation))
}

# Standard error of coefficient `position` from the observed information of the summed log likelihood
# plus the prior's precision; NA when the information is singular.
ageStandardError <- function(design, y, w, precision, position) {
  p <- stats::plogis(as.numeric(design %*% w))
  information <- as.matrix(Matrix::crossprod(design, Matrix::Diagonal(x = p * (1 - p)) %*% design))
  diag(information) <- diag(information) + precision
  unit <- replace(numeric(ncol(information)), position, 1)
  tryCatch(sqrt(solve(information, unit)[[position]]), error = function(e) NA_real_)
}

# Entries, outcomes and age by database, from the matrices of one evaluation.
databaseAges <- function(matrices, pooled, task, minCellCount) {
  age <- match(as.character(ageCovariateId), names(componentCoefficients(pooled$models[[1]])))
  do.call(rbind, lapply(names(matrices), function(site) {
    a <- matrices[[site]]$xMatrix[, age]
    y <- matrices[[site]]$yLabels
    summary <- function(values) {
      if (length(values) < minCellCount) c(NA_real_, NA_real_) else c(mean(values), stats::sd(values))
    }
    all <- summary(a)
    with <- summary(a[y == 1])
    without <- summary(a[y == 0])
    data.frame(task = task, database = site, entries = length(y), outcomes = sum(y), ageMean = all[1], ageSd = all[2],
      ageMeanOutcome = with[1], ageSdOutcome = with[2], ageMeanNoOutcome = without[1], ageSdNoOutcome = without[2])
  }))
}
