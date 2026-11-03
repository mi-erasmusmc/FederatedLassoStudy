#' Many small sites: local-model combinations when databases are split
#'
#' Splits each outer-training database at random into pseudo-sites holding
#' about `outcomesPerSite` outcomes each, fits a local lasso at every
#' pseudo-site exactly as the main analysis fits a local model (three
#' outcome-stratified folds, the same penalty search and solver), and combines
#' them by sample-size-weighted coefficient averaging and probability
#' ensembling. Pooled lasso is unaffected by the split, so the main analysis's
#' pooled lasso is the reference; DualAvg reproduced it, so it also stands for
#' fully coordinated fitting. Pseudo-sites whose local fit fails, for example
#' because a fold has no outcomes, are left out of the combinations and
#' counted.
#'
#' Two further models are fitted once per evaluation:
#' - the whole-database split (each training database one site), which must
#'   reproduce the main analysis's coefficient average;
#' - pooled lasso with unpenalized indicators of the training databases, at
#'   pooled lasso's selected penalty, which estimates the predictor effects
#'   within databases only. Its intercept for a new database is the
#'   sample-size-weighted mean of the database intercepts.
#'
#' With `siteLabelDirectory`, two further splits use natural sites: `natural`
#' splits each training database by the labels that
#' `extras/ExtractSiteLabels.R` wrote (entries without a label form one
#' `unknown` site; a database without a label file stays whole), and
#' `naturalSizeMatched` keeps the same site sizes but assigns entries to them
#' at random. Their difference isolates what real differences between sites
#' add to sample size.
#'
#' The `export` directory holds `pseudoSites.csv` (one row per task, held-out
#' database, split and model: AUC with DeLong difference from pooled lasso,
#' model-based AUC, log loss, calibration-in-the-large and slope, agreement of
#' the coefficients with pooled lasso) and `pseudoSiteCounts.csv` (sites,
#' failed and intercept-only fits, and the distribution of outcomes per site).
#'
#' @param dataRoot,siteIds,taskDirectories As for [runStudy()]; `siteIds` in
#'   the main analysis's order, which numbers the folds.
#' @param mainDirectories Result directories of the main analysis; the one
#'   holding each task's `featureSet` models is used.
#' @param outputDirectory Directory for the `export` tables.
#' @param outcomesPerSite Expected outcomes per pseudo-site.
#' @param tasks Tasks to run.
#' @param featureSet Predictor set.
#' @param holdouts Held-out databases to run.
#' @param cores Parallel R processes for the local fits.
#' @param seed Seed for the splits.
#' @param siteLabelDirectory Directory of natural-site label files,
#'   `<target cohort id>/<database>.csv`; `NULL` for random splits only.
#' @return The export tables, invisibly. The `check` column holds, for the
#'   whole-database coefficient average, its largest absolute coefficient
#'   difference from the main analysis's, and for the database-indicator model,
#'   that of a refit of pooled lasso through the same call from the saved
#'   pooled lasso; both should be near zero.
#' @export
runPseudoSiteAnalysis <- function(dataRoot, siteIds, mainDirectories, outputDirectory,
                                  outcomesPerSite = c(100, 30, 10), tasks = names(studyTasks()),
                                  featureSet = "ageSexPhenotypes", holdouts = siteIds, cores = 4L, seed = 42L,
                                  taskDirectories = NULL, siteLabelDirectory = NULL) {
  cluster <- if (cores > 1L) parallel::makePSOCKcluster(cores) else NULL
  if (!is.null(cluster)) {
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    parallel::clusterEvalQ(cluster, loadNamespace("FederatedLassoStudy"))
  }
  results <- list()
  counts <- list()
  for (task in tasks) {
    sites <- studySites(dataRoot, siteIds, task, featureSet, taskDirectories)[[task]]
    for (holdout in holdouts) {
      fold <- match(holdout, siteIds)
      training <- setdiff(siteIds, holdout)
      directory <- mainDirectory(mainDirectories, task, featureSet, holdout)
      message(format(Sys.time(), "%H:%M"), "  ", task, ", ", holdout, " held out")
      pooled <- FederatedLearning::loadModel(directory, task, holdout, featureSet, "PooledLasso")
      average <- FederatedLearning::loadModel(directory, task, holdout, featureSet, "LocalAvgLasso")
      load <- function(site) {
        plpData <- FederatedLearning::loadClientData(sites[[site]], studyTasks()[[task]]$populationSettings)
        data <- FederatedLearning::createClientMatrix(plpData, pooled$config)
        list(xMatrix = data$xMatrix, yLabels = data$yLabels, n = length(data$yLabels),
          labels = if (!is.null(siteLabelDirectory)) siteLabels(siteLabelDirectory, plpData$population, site))
      }
      test <- load(holdout)
      train <- stats::setNames(lapply(training, load), training)
      reference <- modelPredictions(pooled, test)
      models <- list(pooled = list(w = componentCoefficients(pooled$models[[1]]), probability = reference$probability))
      splits <- c(list(database = lapply(seq_along(train), function(i) rep(1L, train[[i]]$n))),
        stats::setNames(lapply(outcomesPerSite, function(target) lapply(seq_along(train), function(i) {
          y <- train[[i]]$yLabels
          sites <- max(1L, round(sum(y) / target))
          set.seed(seed + 1000L * fold + i)
          sample(rep_len(seq_len(sites), length(y)))
        })), paste0("outcomes", outcomesPerSite)))
      if (!is.null(siteLabelDirectory)) {
        natural <- lapply(train, function(t) if (is.null(t$labels)) rep(1L, t$n) else t$labels)
        splits$natural <- natural
        splits$naturalSizeMatched <- lapply(seq_along(natural), function(i) {
          set.seed(seed + 1000L * fold + i)
          sample(natural[[i]])
        })
      }
      for (split in names(splits)) {
        local <- fitPseudoSites(train, splits[[split]], pooled$args, fold, cluster)
        counts[[length(counts) + 1L]] <- cbind(task = task, heldOutDatabase = holdout, split = split, local$counts)
        if (!length(local$fits)) next
        weights <- local$n / sum(local$n)
        w <- Reduce(`+`, Map(`*`, lapply(local$fits, `[[`, "w"), weights))
        models[[paste0(split, ":LocalAvgLasso")]] <- list(w = w,
          probability = stats::plogis(as.numeric(test$xMatrix %*% w)))
        models[[paste0(split, ":LocalEnsembleLasso")]] <- list(w = NULL,
          probability = as.numeric(do.call(cbind, lapply(local$fits, function(f) {
            stats::plogis(as.numeric(test$xMatrix %*% f$w))
          })) %*% weights))
        if (split == "database") {
          check <- max(abs(w - componentCoefficients(average$models[[1]])))
          models[[paste0(split, ":LocalAvgLasso")]]$check <- check
        }
      }
      models[["databaseIndicators"]] <- databaseIndicatorLasso(train, pooled, test, fold)
      results[[length(results) + 1L]] <- scorePseudoSiteModels(models, test$yLabels, task, holdout)
      rm(train, test)
      gc()
    }
  }
  tables <- list(pseudoSites = do.call(rbind, results), pseudoSiteCounts = do.call(rbind, counts))
  writeExport(outputDirectory, tables)
  invisible(tables)
}

# Natural-site labels of a database's evaluated population, as integers in population order;
# NULL when the database has no label file. Entries without a label form one site.
siteLabels <- function(directory, population, database) {
  file <- file.path(directory, unique(population$targetId), paste0(database, ".csv"))
  if (length(file) != 1L || !file.exists(file)) return(NULL)
  labels <- utils::read.csv(file, colClasses = "character")
  key <- paste(format(population$subjectId, scientific = FALSE, trim = TRUE), format(population$cohortStartDate))
  matched <- match(key, paste(labels$subjectId, labels$cohortStartDate))
  if (mean(!is.na(matched)) < 0.95) {
    stop("Only ", round(100 * mean(!is.na(matched))), "% of ", database, " entries have a site label in ", file,
      call. = FALSE)
  }
  label <- labels$label[matched]
  label[is.na(label) | !nzchar(label)] <- "unknown"
  as.integer(factor(label))
}

# The main-analysis result directory holding pooled lasso for this task, predictor set and holdout.
mainDirectory <- function(directories, task, featureSet, holdout) {
  for (directory in directories) {
    results <- FederatedLearning::readComparisonResults(directory)
    if (any(results$task == task & results$featureSet == featureSet & results$method == "PooledLasso" &
      results$clientId == holdout)) {
      return(directory)
    }
  }
  stop("No main-analysis pooled lasso for ", task, " with ", featureSet, " and ", holdout, " held out",
    call. = FALSE)
}

# Local lasso at every pseudo-site of every training database, seeded per database as the main
# analysis seeds its local fits (so the whole-database split reproduces them).
fitPseudoSites <- function(train, assignment, args, fold, cluster) {
  jobs <- unlist(lapply(seq_along(train), function(i) {
    lapply(sort(unique(assignment[[i]])), function(s) {
      rows <- which(assignment[[i]] == s)
      list(data = list(xMatrix = train[[i]]$xMatrix[rows, , drop = FALSE], yLabels = train[[i]]$yLabels[rows],
        n = length(rows)), database = names(train)[[i]], seed = 42L + fold + i)
    })
  }), recursive = FALSE)
  fits <- if (is.null(cluster)) {
    lapply(jobs, fitPseudoSite, args = args)
  } else {
    parallel::parLapply(cluster, jobs, fitPseudoSite, args = args)
  }
  ok <- vapply(fits, `[[`, logical(1), "ok")
  outcomes <- vapply(fits, `[[`, numeric(1), "outcomes")
  interceptOnly <- vapply(fits[ok], function(f) all(f$w[-1] == 0), logical(1))
  list(fits = fits[ok], n = vapply(fits[ok], `[[`, numeric(1), "n"),
    counts = data.frame(sites = length(fits), failedSites = sum(!ok),
      failedEntries = sum(vapply(fits[!ok], `[[`, numeric(1), "n")), failedOutcomes = sum(outcomes[!ok]),
      interceptOnlySites = sum(interceptOnly),
      interceptOnlyEntryShare = sum(vapply(fits[ok], `[[`, numeric(1), "n")[interceptOnly]) /
        sum(vapply(fits, `[[`, numeric(1), "n")),
      outcomesMin = min(outcomes), outcomesMedian = stats::median(outcomes), outcomesMax = max(outcomes),
      commonFailure = if (any(!ok)) names(which.max(table(vapply(fits[!ok], `[[`, "", "error")))) else NA_character_))
}

# One pseudo-site's local lasso. A package-level function, so that a parallel worker receives only
# its own job's data.
fitPseudoSite <- function(job, args) {
  attempt <- asNamespace("FederatedLearning")$fitLocalBaselineSafely(list(job$data), 1L, job$database, args,
    job$seed, TRUE)
  list(ok = attempt$ok, w = if (attempt$ok) attempt$fit$w else NULL, n = job$data$n, outcomes = sum(job$data$yLabels),
    error = attempt$error)
}

# Pooled lasso with unpenalized training-database indicators at pooled lasso's selected variance,
# fitted as the framework fits pooled lasso. Also refits pooled lasso itself through the same call,
# which must reproduce the saved model.
databaseIndicatorLasso <- function(train, pooled, test, fold) {
  x <- do.call(rbind, lapply(train, `[[`, "xMatrix"))
  y <- unlist(lapply(train, `[[`, "yLabels"), use.names = FALSE)
  database <- rep(seq_along(train), vapply(train, `[[`, numeric(1), "n"))
  indicators <- Matrix::sparseMatrix(i = which(database > 1L), j = database[database > 1L] - 1L,
    x = 1, dims = c(length(y), length(train) - 1L))
  variance <- pooled$models[[1]]$selectedLambda
  fit <- function(design, unpenalized) {
    data <- Cyclops::createCyclopsData(y = y, sx = design, modelType = "lr")
    ids <- as.numeric(Cyclops::getCovariateIds(data))
    model <- Cyclops::fitCyclopsModel(data,
      prior = Cyclops::createPrior("laplace", variance = variance, exclude = ids[unpenalized]),
      control = Cyclops::createControl(maxIterations = 3000L, tolerance = 2e-6, noiseLevel = "silent",
        threads = 1L, seed = 42L + fold), warnings = FALSE)
    as.numeric(stats::coef(model))
  }
  refit <- fit(x, 1L)
  w <- fit(cbind(x, indicators), c(1L, ncol(x) + seq_len(ncol(indicators))))
  effects <- c(0, w[ncol(x) + seq_len(ncol(indicators))])
  shares <- vapply(train, `[[`, numeric(1), "n") / length(y)
  coefficients <- w[seq_len(ncol(x))]
  coefficients[1] <- coefficients[1] + sum(shares * effects)
  list(w = coefficients, probability = stats::plogis(as.numeric(test$xMatrix %*% coefficients)),
    check = max(abs(refit - componentCoefficients(pooled$models[[1]]))))
}

# Discrimination, calibration and coefficient agreement of each model against pooled lasso.
scorePseudoSiteModels <- function(models, y, task, holdout) {
  delong <- delongAuc(y, do.call(cbind, lapply(models, `[[`, "probability")))
  pooledW <- models$pooled$w
  slopes <- !names(pooledW) %in% c("(Intercept)", "1002")
  z <- stats::qnorm(0.975)
  do.call(rbind, lapply(seq_along(models), function(k) {
    m <- models[[k]]
    lp <- stats::qlogis(pmin(pmax(m$probability, 1e-15), 1 - 1e-15))
    se <- sqrt(max(delong$covariance[k, k] + delong$covariance[1, 1] - 2 * delong$covariance[k, 1], 0))
    difference <- delong$auc[k] - delong$auc[1]
    agreement <- if (is.null(m$w)) list(cor = NA_real_, shrink = NA_real_, jaccard = NA_real_, nonzero = NA_real_) else {
      a <- m$w[slopes]
      b <- pooledW[slopes]
      shared <- a != 0 & b != 0
      list(cor = if (length(a) > 1L && stats::sd(a) > 0 && stats::sd(b) > 0) stats::cor(a, b) else NA_real_,
        shrink = if (any(shared)) stats::median(abs(a[shared]) / abs(b[shared])) else NA_real_,
        jaccard = sum(shared) / max(1, sum(a != 0 | b != 0)), nonzero = sum(a != 0))
    }
    parts <- strsplit(names(models)[[k]], ":", fixed = TRUE)[[1]]
    data.frame(task = task, heldOutDatabase = holdout, split = parts[[1]],
      method = if (length(parts) > 1L) parts[[2]] else parts[[1]], n = length(y), outcomes = sum(y),
      auc = delong$auc[k], aucDifference = difference, differenceLower = difference - z * se,
      differenceUpper = difference + z * se, modelBasedAuc = modelBasedAuc(m$probability),
      logLoss = mean(-(y * log(pmax(m$probability, 1e-15)) + (1 - y) * log(pmax(1 - m$probability, 1e-15)))),
      calibrationInTheLarge = unname(stats::coef(stats::glm(y ~ 1 + offset(lp), family = stats::binomial()))[1]),
      calibrationSlope = unname(stats::coef(stats::glm(y ~ lp, family = stats::binomial()))[2]),
      coefficientCorrelation = agreement$cor, sharedMagnitudeRatio = agreement$shrink,
      supportJaccard = agreement$jaccard, nonzero = agreement$nonzero,
      check = if (is.null(m$check)) NA_real_ else m$check, stringsAsFactors = FALSE)
  }))
}
