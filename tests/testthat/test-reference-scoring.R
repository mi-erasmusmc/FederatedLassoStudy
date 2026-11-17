test_that("saved models are scored against the exported converged pooled lasso", {
  set.seed(11)
  n <- 4000
  x <- cbind(1, age = stats::runif(n, 0.55, 0.85), female = stats::rbinom(n, 1, 0.5),
    b = stats::rbinom(n, 1, 0.2))
  w <- c(-5, 3, 0.4, 0.6)
  y <- stats::rbinom(n, 1, stats::plogis(as.numeric(x %*% w)))
  columns <- c("(Intercept)", "1002", "8532001", "101")
  component <- function(coefficients) list(coefficients = data.frame(matrixColumn = seq_along(coefficients),
    coefficientSharedScale = coefficients, isIntercept = c(TRUE, FALSE, FALSE, FALSE),
    covariateId = c(NA, 1002, 8532001, 101)))
  saved <- list(DualAvg = list(models = list(component(w)), config = list(intercept = TRUE)),
    ADAP2 = list(models = list(component(w + c(0.5, -1, 0, 0.3))), config = list(intercept = TRUE)))
  reference <- file.path(tempfile("reference-"), "export")
  dir.create(reference, recursive = TRUE)
  # Stored in a shuffled order, as the export may list covariates in any order.
  utils::write.csv(data.frame(task = "BleedWarfarin", featureSet = "ageSexPhenotypes", heldOutDatabase = "C",
    model = "convergedPooled", sourceDatabase = NA, covariateId = rev(columns), coefficient = rev(w)),
    file.path(reference, "referenceCoefficients.csv"), row.names = FALSE)
  local_mocked_bindings(
    scoredModels = function(directories) data.frame(directory = "d", task = "BleedWarfarin",
      featureSet = "ageSexPhenotypes", fold = 3L, heldOutDatabase = "C", method = c("DualAvg", "ADAP2"),
      stringsAsFactors = FALSE),
    studySites = function(...) list(BleedWarfarin = c(A = "a", B = "b", C = "c")))
  local_mocked_bindings(
    loadModel = function(directory, task, holdout, featureSet, method) saved[[method]],
    loadClientData = function(path, population) list(path = path),
    createClientMatrix = function(plpData, config) list(xMatrix = Matrix::Matrix(x, sparse = TRUE), yLabels = y),
    .package = "FederatedLearning")
  out <- tempfile("scored-")
  tables <- runReferenceScoring("data", c("A", "B", "C"), "d", dirname(reference), out)
  scores <- tables$referenceScores
  expect_identical(scores$model, c("convergedPooled", "DualAvg", "ADAP2"))
  expect_identical(scores$kind, c("converged", "saved", "saved"))
  expect_equal(scores$aucDifference[1:2], c(0, 0), tolerance = 1e-12)
  expect_equal(scores$linearPredictorPearson[2], 1, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(scores$aucDifference[3], 0)))
  expect_true(all(is.finite(scores$calibrationInTheLarge)))
  expect_true(file.exists(file.path(out, "export", "referenceScores.csv")))
  expect_setequal(unique(tables$scoredCoefficients$model), c("DualAvg", "ADAP2"))
  expect_true(all(c("convergedPooled", "DualAvg") %in% tables$referenceCurves$model))

  saved$DualAvg$models[[1]]$coefficients$covariateId[4] <- 999
  saved$ADAP2 <- NULL
  local_mocked_bindings(scoredModels = function(directories) data.frame(directory = "d", task = "BleedWarfarin",
    featureSet = "ageSexPhenotypes", fold = 3L, heldOutDatabase = "C", method = "DualAvg", stringsAsFactors = FALSE))
  expect_error(runReferenceScoring("data", c("A", "B", "C"), "d", dirname(reference), tempfile()), "same predictors")
})

test_that("converged models are rescored from their exported coefficients", {
  set.seed(12)
  n <- 3000
  x <- cbind(1, age = stats::runif(n, 0.55, 0.85), female = stats::rbinom(n, 1, 0.5),
    b = stats::rbinom(n, 1, 0.2))
  y <- stats::rbinom(n, 1, stats::plogis(as.numeric(x %*% c(-5, 3, 0.4, 0.6))))
  columns <- c("(Intercept)", "1002", "8532001", "101")
  local <- list(A = c(-4.8, 2.7, 0.3, 0.5), B = c(-5.3, 3.4, 0.5, 0.8))
  weights <- c(A = 0.75, B = 0.25)
  component <- function(coefficients, site = NULL, weight = NULL) c(list(coefficients = data.frame(
    matrixColumn = seq_along(coefficients), coefficientSharedScale = coefficients,
    isIntercept = c(TRUE, FALSE, FALSE, FALSE), covariateId = c(NA, 1002, 8532001, 101))),
    if (!is.null(site)) list(sourceClientId = site, weight = weight))
  pooled <- c(-5, 3, 0.4, 0.6)
  saved <- list(
    LocalEnsembleLasso = list(models = lapply(c("A", "B"), function(s) component(local[[s]] + 0.1, s, weights[[s]])),
      config = list(intercept = TRUE)),
    BiggestSiteLasso = list(models = list(component(local$A + 0.1, "A", 1)), config = list(intercept = TRUE)))
  average <- weights[["A"]] * local$A + weights[["B"]] * local$B
  rows <- function(model, w, source = NA) data.frame(task = "BleedWarfarin", featureSet = "ageSexPhenotypes",
    heldOutDatabase = "C", model = model, sourceDatabase = source, covariateId = columns, coefficient = w)
  reference <- file.path(tempfile("reference-"), "export")
  dir.create(reference, recursive = TRUE)
  utils::write.csv(rbind(rows("convergedPooled", pooled), rows("convergedLocalAvgLasso", average),
    rows("convergedLocal", local$A, "A"), rows("convergedLocal", local$B, "B")),
    file.path(reference, "referenceCoefficients.csv"), row.names = FALSE)
  local_mocked_bindings(
    scoredModels = function(directories) data.frame(directory = "d", task = "BleedWarfarin",
      featureSet = "ageSexPhenotypes", fold = 3L, heldOutDatabase = "C",
      method = c("LocalEnsembleLasso", "BiggestSiteLasso"), stringsAsFactors = FALSE),
    studySites = function(...) list(BleedWarfarin = c(A = "a", B = "b", C = "c")))
  local_mocked_bindings(
    loadModel = function(directory, task, holdout, featureSet, method) saved[[method]],
    loadClientData = function(path, population) list(path = path),
    createClientMatrix = function(plpData, config) list(xMatrix = Matrix::Matrix(x, sparse = TRUE), yLabels = y),
    .package = "FederatedLearning")
  tables <- runReferenceScoring("data", c("A", "B", "C"), "d", dirname(reference), tempfile("scored-"))
  scores <- tables$referenceScores
  expect_identical(scores$model, c("convergedPooled", "LocalEnsembleLasso", "BiggestSiteLasso",
    "convergedLocalAvgLasso", "convergedBiggestSiteLasso", "convergedLocalEnsembleLasso", "convergedLocal",
    "convergedLocal"))
  expect_identical(scores$sourceDatabase[7:8], c("A", "B"))
  expect_true(all(scores$kind[4:8] == "converged"))
  auc <- function(score) as.numeric(pROC::auc(pROC::roc(y, score, levels = c(0, 1), direction = "<", quiet = TRUE)))
  ensemble <- weights[["A"]] * stats::plogis(as.numeric(x %*% local$A)) +
    weights[["B"]] * stats::plogis(as.numeric(x %*% local$B))
  expect_equal(scores$auc[6], auc(ensemble), tolerance = 1e-10)
  expect_equal(scores$auc[4], auc(as.numeric(x %*% average)), tolerance = 1e-10)
  expect_equal(scores$auc[5], scores$auc[7], tolerance = 1e-12)
  expect_true(all(c("convergedLocalEnsembleLasso", "convergedLocal") %in% tables$referenceCurves$model))
  # Without scoreConverged only the reference and the saved models are scored.
  plain <- runReferenceScoring("data", c("A", "B", "C"), "d", dirname(reference), tempfile("scored-"),
    scoreConverged = FALSE)
  expect_identical(plain$referenceScores$model, c("convergedPooled", "LocalEnsembleLasso", "BiggestSiteLasso"))
})
