fakeModel <- function(...) {
  vectors <- list(...)
  if (is.null(names(vectors))) names(vectors) <- rep("", length(vectors))
  list(models = lapply(seq_along(vectors), function(i) {
    w <- vectors[[i]]
    source <- names(vectors)[[i]]
    list(sourceClientId = if (nzchar(source)) source else NULL,
      coefficients = data.frame(matrixColumn = seq_along(w), covariateId = c(NA, names(w)[-1]),
        isIntercept = c(TRUE, rep(FALSE, length(w) - 1L)), coefficientSharedScale = unname(w)))
  }))
}

test_that("coefficient agreement pairs each model with pooled lasso of the same fold", {
  models <- list(
    PooledLasso = fakeModel(c("(Intercept)" = -2, "1002" = 5, "8532001" = 0.2, "11" = 0.5, "12" = 0, "13" = 1)),
    DualAvg = fakeModel(c("(Intercept)" = -1, "1002" = 9, "8532001" = 0.2, "11" = 0.4, "12" = 0.3, "13" = 1)),
    # The feature filter dropped predictor 13 and kept 14.
    ADAP2 = fakeModel(c("(Intercept)" = -2, "1002" = 5, "8532001" = 0, "11" = 0, "12" = 0, "14" = 2)),
    LocalSiteLasso = fakeModel(A = c("(Intercept)" = 0, "1002" = 1, "8532001" = 0.2, "11" = 0.5, "12" = 0, "13" = 1),
      B = c("(Intercept)" = 0, "1002" = 1, "8532001" = 0, "11" = 0, "12" = 0, "13" = 0))
  )
  loaded <- list()
  local_mocked_bindings(
    readComparisonResults = function(outputDirectory) data.frame(task = "DvtTetraplegia",
      featureSet = "ageSexPhenotypes", fold = 1L, clientId = "C",
      method = c("PooledLasso", "DualAvg", "ADAP2", "LocalSiteLasso", "LocalSiteLasso", "LocalEnsembleLasso", "ODAL1"),
      modelFile = c("p.rds", "d.rds", "a.rds", "l.rds", "l.rds", "e.rds", NA),
      error = c(NA, NA, NA, NA, NA, NA, "failed")),
    loadModel = function(outputDirectory, task, holdout, feature, method) {
      loaded[[length(loaded) + 1L]] <<- method
      models[[method]]
    },
    .package = "FederatedLearning"
  )
  out <- tempfile()
  on.exit(unlink(out, recursive = TRUE), add = TRUE)
  agreement <- runCoefficientAgreement("run1", out)
  expect_setequal(unlist(loaded), names(models))
  expect_identical(agreement$method, c("DualAvg", "ADAP2", "LocalSiteLasso", "LocalSiteLasso"))
  expect_identical(agreement$sourceDatabase, c(NA, NA, "A", "B"))
  expect_identical(agreement$heldOutDatabase[[1]], "C")

  # Age and intercept are left out: four binary predictors remain.
  dual <- agreement[1, ]
  expect_identical(dual$sharedPredictors, 4L)
  expect_equal(dual$rmsDifference, sqrt((0.1^2 + 0.3^2) / 4))
  expect_equal(dual$correlation, stats::cor(c(0.2, 0.4, 0.3, 1), c(0.2, 0.5, 0, 1)))
  expect_identical(dual$maxAbsDifferenceCovariateId, "12")
  expect_identical(c(dual$nonzero, dual$pooledNonzero, dual$bothNonzero), c(4L, 3L, 3L))
  expect_equal(dual$supportJaccard, 3 / 4)

  adap <- agreement[2, ]
  expect_identical(c(adap$sharedPredictors, adap$onlyModelPredictors, adap$onlyPooledPredictors), c(3L, 1L, 1L))
  expect_identical(c(adap$nonzero, adap$bothNonzero), c(1L, 0L))
  expect_equal(adap$supportJaccard, 0)

  expect_equal(agreement$correlation[[3]], 1)
  expect_equal(agreement$rmsDifference[[3]], 0)
  # A local model with no nonzero predictors has no correlation.
  expect_true(is.na(agreement$correlation[[4]]))
  expect_equal(agreement$supportJaccard[[4]], 0)

  coefficients <- utils::read.csv(file.path(out, "export", "coefficients.csv"))
  expect_identical(nrow(coefficients), 5L * 6L)
  expect_true(all(c("(Intercept)", "1002") %in% coefficients$covariateId))
  expect_identical(nrow(utils::read.csv(file.path(out, "export", "coefficientAgreement.csv"))), 4L)
})

test_that("coefficient agreement rejects a shared model with several vectors", {
  local_mocked_bindings(
    readComparisonResults = function(outputDirectory) data.frame(task = "DvtTetraplegia",
      featureSet = "ageSexPhenotypes", fold = 1L, clientId = "C", method = c("PooledLasso", "DualAvg"),
      modelFile = c("p.rds", "d.rds")),
    loadModel = function(outputDirectory, task, holdout, feature, method) {
      w <- c("(Intercept)" = 0, "11" = 1)
      if (method == "DualAvg") fakeModel(A = w, B = w) else fakeModel(w)
    },
    .package = "FederatedLearning"
  )
  expect_error(runCoefficientAgreement("run1", tempfile()), "DualAvg has 2 coefficient vectors")
})
