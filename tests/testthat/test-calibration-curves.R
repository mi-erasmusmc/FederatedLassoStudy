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

test_that("the calibration figures draw one figure per task and predictor set", {
  skip_if_not_installed("ggplot2")
  curve <- function(task, featureSet, database, method, shift = 0) {
    predicted <- exp(seq(log(0.005), log(0.1), length.out = 20))
    data.frame(resultsDirectory = "run1", task = task, featureSet = featureSet, fold = 1L,
      heldOutDatabase = database, method = method, sourceDatabase = NA, predicted = predicted,
      observed = predicted * (1 + shift), lower = predicted * 0.8, upper = predicted * 1.2)
  }
  curves <- rbind(
    curve("DvtTetraplegia", "ageSexPhenotypes", "A", "PooledLasso"),
    curve("DvtTetraplegia", "ageSexPhenotypes", "A", "DualAvg", 0.1),
    curve("DvtTetraplegia", "ageSexPhenotypes", "B", "PooledLasso"),
    curve("DvtTetraplegia", "ageSexPhenotypes", "B", "DualAvg"),
    curve("DvtTetraplegia", "ageSexConditions", "A", "PooledLasso"),
    curve("DvtTetraplegia", "ageSexConditions", "A", "LocalSiteLasso"))
  curves$predicted[[30]] <- 0.999
  input <- tempfile()
  out <- tempfile()
  on.exit(unlink(c(input, out), recursive = TRUE), add = TRUE)
  dir.create(file.path(input, "export"), recursive = TRUE)
  utils::write.csv(curves, file.path(input, "export", "calibrationCurves.csv"), row.names = FALSE)
  paths <- writeCalibrationFigures(input, out)
  expect_setequal(basename(paths), paste0("calibration_DvtTetraplegia_",
    rep(c("ageSexPhenotypes", "ageSexConditions"), each = 2), c(".pdf", ".png")))
  expect_true(all(file.size(paths) > 0))
  plot <- calibrationFigure(curves[curves$featureSet == "ageSexPhenotypes", ])
  expect_identical(levels(plot$data$panel), c("Pooled lasso", "DualAvg"))
  offScale <- plot$layers[[5]]$data
  expect_identical(offScale$label, "1 point off scale")
})
