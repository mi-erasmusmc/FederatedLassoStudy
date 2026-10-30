test_that("the calibration export collects every saved curve with its model identity", {
  curve <- function(n) data.frame(predicted = seq(0.1, 0.2, length.out = n), observed = 0.15,
    lower = 0.1, upper = 0.2)
  models <- list(
    PooledLasso = list(calibrationCurves = cbind(sourceClientId = NA, clientId = "C", curve(3))),
    DualAvg = list(calibrationCurves = cbind(clientId = "C", client = 1, curve(2))),
    LocalSiteLasso = list(calibrationCurves = cbind(sourceClientId = c("A", "A", "B", "B"), clientId = "C",
      curve(4))),
    ODAL1 = list(calibrationCurves = NULL)
  )
  local_mocked_bindings(
    readComparisonResults = function(outputDirectory) data.frame(task = "DvtTetraplegia",
      featureSet = "ageSexPhenotypes", fold = 1L, clientId = "C",
      method = c("PooledLasso", "DualAvg", "LocalSiteLasso", "LocalSiteLasso", "ODAL1", "ADAP1"),
      modelFile = c("p.rds", "d.rds", "l.rds", "l.rds", "o.rds", NA), error = c(NA, NA, NA, NA, NA, "failed")),
    loadModel = function(outputDirectory, task, holdout, feature, method) models[[method]],
    .package = "FederatedLearning"
  )
  out <- tempfile()
  on.exit(unlink(out, recursive = TRUE), add = TRUE)
  curves <- runCalibrationCurveExport("results/run1", out)
  expect_identical(as.vector(table(curves$method)[c("PooledLasso", "DualAvg", "LocalSiteLasso")]), c(3L, 2L, 4L))
  expect_identical(curves$sourceDatabase[curves$method == "LocalSiteLasso"], c("A", "A", "B", "B"))
  expect_true(all(is.na(curves$sourceDatabase[curves$method != "LocalSiteLasso"])))
  expect_identical(unique(curves$resultsDirectory), "run1")
  expect_identical(nrow(utils::read.csv(file.path(out, "export", "calibrationCurves.csv"))), 9L)
})
