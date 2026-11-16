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
