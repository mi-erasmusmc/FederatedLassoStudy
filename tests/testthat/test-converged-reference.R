simulatedDatabase <- function(n, seed) {
  set.seed(seed)
  age <- stats::runif(n, 0.55, 0.85)
  prevalent <- stats::rbinom(n, 1, 0.97)
  other <- matrix(stats::rbinom(n * 4, 1, 0.2), n)
  female <- stats::rbinom(n, 1, 0.5)
  eta <- -6 + 3 * age + 0.4 * female + 1.2 * prevalent + as.numeric(other %*% c(0.6, -0.4, 0, 0.3))
  x <- Matrix::Matrix(cbind(1, age, female, prevalent, other), sparse = TRUE)
  list(xMatrix = methods::as(methods::as(x, "generalMatrix"), "CsparseMatrix"),
    yLabels = stats::rbinom(n, 1, stats::plogis(eta)))
}
referenceColumns <- c("(Intercept)", "1002", "8532001", "1160049", "101", "102", "103", "104")
savedComponent <- function(w, variance = NA_real_, weight = NA_real_, source = NA_character_) {
  list(sourceClientId = source, weight = weight, selectedLambda = variance,
    coefficients = data.frame(matrixColumn = seq_along(w), coefficientSharedScale = unname(w),
      isIntercept = c(TRUE, rep(FALSE, length(w) - 1L)), covariateId = c(NA, as.numeric(referenceColumns[-1]))))
}
# Saved models as the main run stores them, fitted loosely so that they stop short of their optimum.
fakeSavedModels <- function(train, loose) {
  x <- do.call(rbind, lapply(train, `[[`, "xMatrix"))
  y <- unlist(lapply(train, `[[`, "yLabels"))
  lambda <- 0.002
  pooled <- cyclopsFit(x, y, "laplace", 1L, loose, lambda)$w
  dualAvg <- cyclopsFit(x, y, "laplace", 1L, loose, 2 * lambda)$w
  sizes <- vapply(train, function(t) length(t$yLabels), numeric(1))
  weights <- sizes / sum(sizes)
  components <- Map(function(t, site, weight) {
    w <- cyclopsFit(t$xMatrix, t$yLabels, "laplace", 1L, loose, 0.004)$w
    savedComponent(w, 2 / (0.004 * length(t$yLabels))^2, weight, site)
  }, train, names(train), weights)
  average <- Reduce(`+`, Map(function(component, weight) weight * componentCoefficients(component), components, weights))
  list(PooledLasso = list(models = list(savedComponent(pooled, 2 / (lambda * length(y))^2)),
      config = list(intercept = TRUE), args = list()),
    DualAvg = list(models = list(savedComponent(dualAvg, 2 * lambda))),
    LocalAvgLasso = list(models = list(savedComponent(average)), components = unname(components)),
    LocalEnsembleLasso = list(models = unname(components)),
    BiggestSiteLasso = list(models = list(components[[which.max(sizes)]])))
}

test_that("converged refits meet the optimality conditions and the converged combinations follow the saved weights", {
  train <- list(A = simulatedDatabase(4000, 1), B = simulatedDatabase(2500, 2), C = simulatedDatabase(1500, 3))
  test <- simulatedDatabase(3000, 4)
  loose <- list(tolerance = 1e-3, maxIterations = 50L, threads = 1L)
  saved <- fakeSavedModels(train, loose)
  control <- list(tolerance = 1e-12, maxIterations = 20000L, threads = 1L)
  out <- referenceEvaluation(train, test, saved, control)
  fits <- out$referenceFits
  converged <- fits[startsWith(fits$fit, "converged"), ]
  expect_setequal(unique(converged$fit), c("convergedPooled", "convergedDualAvgPenalty", "convergedLocal", "convergedTarget"))
  expect_true(all(is.na(converged$error)))
  expect_true(all(converged$kktMax < 1e-6))
  expect_lt(fits$check[fits$fit == "averageCheck"], 1e-12)
  # Each local fit is refitted once and is no worse than the saved one on its own objective.
  expect_setequal(converged$sourceDatabase[converged$fit == "convergedLocal"], c("A", "B", "C"))
  savedLocal <- fits[fits$fit == "savedLocal", ]
  convergedLocal <- converged[converged$fit == "convergedLocal", ]
  expect_true(all(convergedLocal$objective <= savedLocal$objective[match(convergedLocal$sourceDatabase,
    savedLocal$sourceDatabase)] + 1e-12))
  # The converged pooled lasso is the plain lasso solution on the original columns.
  x <- do.call(rbind, lapply(train, `[[`, "xMatrix"))
  y <- unlist(lapply(train, `[[`, "yLabels"))
  plain <- cyclopsFit(x, y, "laplace", 1L, control, 0.002)$w
  coefficients <- out$referenceCoefficients
  reference <- coefficients$coefficient[coefficients$model == "convergedPooled"]
  expect_equal(reference, plain, tolerance = 1e-5)
  local <- coefficients[coefficients$model == "convergedLocal", ]
  weights <- vapply(saved$LocalAvgLasso$components, `[[`, numeric(1), "weight")
  sources <- vapply(saved$LocalAvgLasso$components, `[[`, "", "sourceClientId")
  expected <- Reduce(`+`, Map(function(site, weight) weight * local$coefficient[local$sourceDatabase == site], sources, weights))
  expect_equal(coefficients$coefficient[coefficients$model == "convergedLocalAvgLasso"], expected)

  scores <- out$referenceScores
  expect_identical(scores$model[[1]], "convergedPooled")
  expect_equal(scores$aucDifference[[1]], 0)
  expect_equal(scores$linearPredictorPearson[[1]], 1)
  expect_setequal(scores$model, c("convergedPooled", names(saved), "convergedDualAvgPenalty", "convergedLocalAvgLasso",
    "convergedBiggestSiteLasso", "convergedLocalEnsembleLasso", "convergedLocal", "TargetLasso", "convergedTargetLasso"))
  expect_setequal(scores$sourceDatabase[scores$model == "convergedLocal"], c("A", "B", "C"))
  expect_identical(scores$kind[scores$model == "PooledLasso"], "saved")
  # The converged largest-database model is the converged local fit of the largest database.
  biggest <- local$coefficient[local$sourceDatabase == "A"]
  expect_equal(scores$auc[scores$model == "convergedBiggestSiteLasso"],
    delongAuc(test$yLabels, cbind(stats::plogis(as.numeric(test$xMatrix %*% biggest))))$auc)
  expect_true(all(scores$top1Overlap >= 0 & scores$top1Overlap <= 1))

  # The target-database lasso: five tuned folds, each refitted to convergence at its own penalty.
  tuned <- fits[fits$fit == "tunedTarget", ]
  refitted <- fits[fits$fit == "convergedTarget", ]
  expect_identical(nrow(tuned), 5L)
  expect_equal(refitted$lambda, tuned$lambda)
  expect_true(all(refitted$objective <= tuned$objective + 1e-12))
  # A model has a curve on the framework's 50 percentiles, and an ICI, unless it gives fewer than five
  # distinct predictions; the reference has one.
  curves <- out$referenceCurves
  drawn <- paste(scores$model, scores$sourceDatabase) %in% unique(paste(curves$model, curves$sourceDatabase))
  expect_identical(is.finite(scores$calibrationIci), drawn)
  expect_true(drawn[[1]])
  expect_true(all(table(paste(curves$model, curves$sourceDatabase)) <= 50L))
  expect_true(all(curves$lower <= curves$observed & curves$observed <= curves$upper))
})

test_that("a model on other predictors or from an unknown database is refused", {
  train <- list(A = simulatedDatabase(1500, 5), B = simulatedDatabase(1500, 6))
  test <- simulatedDatabase(1000, 7)
  loose <- list(tolerance = 1e-3, maxIterations = 50L, threads = 1L)
  control <- list(tolerance = 1e-10, maxIterations = 20000L, threads = 1L)
  saved <- fakeSavedModels(train, loose)
  other <- saved
  other$DualAvg$models[[1]]$coefficients$covariateId[[5]] <- 999
  expect_error(referenceEvaluation(train, test, other, control), "does not use pooled lasso's predictors")
  stranger <- saved
  stranger$BiggestSiteLasso$models[[1]]$sourceClientId <- "Z"
  expect_error(referenceEvaluation(train, test, stranger, control), "not a training database")
})
