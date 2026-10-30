fakeResults <- function() {
  methods <- c("PooledLasso", "DualAvg", "ODAL1", "ADAP1", "LocalAvgLasso", "LocalEnsembleLasso",
    "BiggestSiteLasso", "LocalSiteLasso", "LocalSiteLasso")
  labels <- c("default", "default", "default", "default", "coefficientAverage=sampleSize",
    "predictionAverage=sampleSize", "source=B", "source=A", "source=B")
  fold <- function(fold, heldOut, pooledAuc) {
    data.frame(task = "DvtTetraplegia", featureSet = "ageSexPhenotypes", fold = fold, clientId = heldOut,
      method = methods, configLabel = labels, error = NA, n = 100 * fold, outcomes = 10 * fold, p = 52,
      nonzeroPredictors = 40, auc = pooledAuc + c(0, 0.01, -0.02, 0.005, 0, 0.002, -0.01, -0.05, -0.01),
      calibrationInTheLarge = 0.1 * fold, calibrationSlope = 1, calibrationIci = 0.01,
      optimConvergence = c(NA, NA, if (fold == 1) 1 else 0, NA, NA, NA, NA, NA, NA),
      adapFinalConverged = c(NA, NA, NA, 1, NA, NA, NA, NA, NA),
      pooledObjectiveGap = c(0, 1e-4, if (fold == 1) 2 else 0.01, if (fold == 2) 0.06 else 0, 0.001, NA, 0.002, NA, NA),
      communicationBasis = "protocol", exchanges = 1, messages = 2, numbers = 3, bytes = 4,
      stringsAsFactors = FALSE)
  }
  rbind(fold(1L, "A", 0.70), fold(2L, "B", 0.60))
}

test_that("the report pairs each method with pooled lasso and flags failed fits", {
  results <- fakeResults()
  local_mocked_bindings(
    readComparisonResults = function(outputDirectory) {
      if (outputDirectory == "run1") results[results$fold == 1L, ] else results[results$fold == 2L, ]
    },
    .package = "FederatedLearning"
  )
  out <- tempfile()
  on.exit(unlink(out, recursive = TRUE), add = TRUE)
  evaluations <- writeStudyReport(c("run1", "run2"), out)
  expect_equal(evaluations$aucDifference[evaluations$method == "DualAvg"], c(0.01, 0.01))
  flagged <- evaluations[evaluations$flagged, c("fold", "method", "flagReason")]
  expect_identical(flagged$method, c("ODAL1", "ADAP1"))
  expect_identical(flagged$flagReason, c("check failed; training log loss", "training log loss"))
  expect_identical(evaluations$sourceDatabase[evaluations$method == "BiggestSiteLasso"], c("B", "B"))
  expect_identical(evaluations$retainedPredictors[[1]], 51)

  discrimination <- utils::read.csv(file.path(out, "discriminationSummary.csv"))
  odal <- discrimination[discrimination$method == "ODAL1" & discrimination$task == "overall", ]
  expect_equal(odal$meanAuc, 0.63)
  expect_equal(odal$meanAucDifference, -0.02)
  expect_identical(c(odal$pairs, odal$flagged, odal$unflaggedPairs), c(2L, 1L, 1L))
  expect_false("LocalSiteLasso" %in% discrimination$method)

  population <- utils::read.csv(file.path(out, "population.csv"))
  expect_identical(population$database, c("A", "B"))
  expect_equal(population$outcomePercent, c(10, 10))

  local <- utils::read.csv(file.path(out, "localBaselines.csv"))
  expect_equal(local$minSourceAuc, c(0.65, 0.55))
  expect_identical(local$bestSourceDatabase, c("B", "B"))
  expect_equal(local$coefficientAverageAuc, c(0.70, 0.60))
  expect_identical(local$largestSiteDatabase, c("B", "B"))
})

test_that("the report rejects an evaluation present in two directories", {
  local_mocked_bindings(
    readComparisonResults = function(outputDirectory) fakeResults(),
    .package = "FederatedLearning"
  )
  expect_error(writeStudyReport(c("run1", "run2"), tempfile()), "more than one results directory")
})

test_that("the report figures are drawn from the report tables", {
  skip_if_not_installed("ggplot2")
  results <- fakeResults()
  results$featureSet[results$fold == 2L] <- "ageSexConditions"
  results <- results[!(results$method == "ODAL1" & results$fold == 2L), ]
  local_mocked_bindings(
    readComparisonResults = function(outputDirectory) results,
    .package = "FederatedLearning"
  )
  out <- tempfile()
  on.exit(unlink(out, recursive = TRUE), add = TRUE)
  writeStudyReport("run1", out)
  paths <- writeStudyFigures(out)
  expect_setequal(basename(paths), c("discrimination.pdf", "discrimination.png", "localBaselines.pdf",
    "localBaselines.png"))
  expect_true(all(file.size(paths) > 0))

  predicted <- exp(seq(log(0.005), log(0.1), length.out = 20))
  curves <- merge(data.frame(task = "BleedWarfarin", featureSet = "ageSexPhenotypes", heldOutDatabase = c("A", "B"),
    stringsAsFactors = FALSE), data.frame(method = c("PooledLasso", "DualAvg", "ODAL1", "ADAP1"),
    stringsAsFactors = FALSE))
  curves <- merge(curves, data.frame(predicted = predicted))
  curves$observed <- curves$predicted
  curves$predicted[curves$method == "ODAL1" & curves$heldOutDatabase == "A" & curves$predicted > 0.09] <- 0.999
  curveDirectory <- file.path(out, "curves")
  dir.create(file.path(curveDirectory, "export"), recursive = TRUE)
  utils::write.csv(curves, file.path(curveDirectory, "export", "calibrationCurves.csv"), row.names = FALSE)
  paths <- writeStudyFigures(out, calibrationCurveDirectory = curveDirectory)
  expect_true(all(file.size(file.path(out, c("calibration.pdf", "calibration.png"))) > 0))
})
