test_that("the preprocessing audit runs every task with the study settings", {
  calls <- list()
  local_mocked_bindings(
    runComparison = function(taskSettings, sites, featureSettings, methodSettings,
                             outputDirectory, runSettings) {
      calls[[length(calls) + 1L]] <<- list(task = taskSettings$name, sites = sites,
        features = featureSettings, auditOnly = runSettings$auditOnly, hosts = runSettings$hosts)
      invisible(NULL)
    },
    .package = "FederatedLearning"
  )
  local_mocked_bindings(writePreprocessingExport = function(...) data.frame())
  runPreprocessingAudit("root", c("A", "B", "C"), "out",
    taskDirectories = c(BleedWarfarin = "warfarin", other = "ignored"))
  expect_identical(vapply(calls, `[[`, "", "task"), names(studyTasks()))
  expect_identical(calls[[1]]$sites, c(A = "root/DvtTetraplegia/A", B = "root/DvtTetraplegia/B",
    C = "root/DvtTetraplegia/C"))
  expect_identical(unname(calls[[3]]$sites), file.path("root", "warfarin", c("A", "B", "C")))
  expect_identical(calls[[1]]$features, studyFeatureSettings())
  expect_true(all(vapply(calls, `[[`, logical(1), "auditOnly")))
  expect_error(runPreprocessingAudit("root", c("A", "B"), "out", tasks = "Other"), "studyTasks")
  expect_error(runPreprocessingAudit("root", "A", "out"), "two distinct sites")
})

test_that("the export holds predictor ids and counts, not patient counts", {
  output <- tempfile("audit-")
  on.exit(unlink(output, recursive = TRUE), add = TRUE)
  directory <- file.path(output, "preprocessing")
  dir.create(directory, recursive = TRUE)
  audit <- function(antibiotic) {
    ids <- c("1002", "8532001", "1152049", if (antibiotic) "999001049")
    data.frame(covariateId = ids, trainingRows = 500, nonzeroRows = 3,
      removedRare = ids == "1152049", removedNearConstant = FALSE, retained = ids != "1152049")
  }
  saveRDS(list(audit = audit(TRUE)), file.path(directory, "DvtTetraplegia_fold2_ageSexPhenotypes.rds"))
  saveRDS(list(audit = audit(FALSE)), file.path(directory, "BleedWarfarin_fold1_ageSexPhenotypes.rds"))
  utils::write.csv(data.frame(task = c("DvtTetraplegia", "BleedWarfarin", "BleedTetraPlegia"),
    fold = c(2L, 1L, 1L), featureSet = "ageSexPhenotypes", candidatePredictors = c(4L, 3L, 9L),
    retainedPredictors = c(3L, 2L, 9L), removedRare = 1L, removedNearConstant = 0L,
    requiresRefit = TRUE, fingerprint = c("f1", "f2", "f3"), referenceStatus = "not_requested"),
    file.path(directory, "summary.csv"), row.names = FALSE)

  summary <- writePreprocessingExport(output, c("A", "B", "C"), c("DvtTetraplegia", "BleedWarfarin"))
  expect_identical(summary$heldOutSite, c("B", "A"))
  expect_identical(summary$trainingSites, c("A;C", "B;C"))
  expect_identical(summary$antibioticRetained, c(TRUE, NA))
  predictors <- utils::read.csv(file.path(output, "export", "predictors.csv"), stringsAsFactors = FALSE)
  expect_identical(names(predictors), c("task", "heldOutSite", "featureSet", "covariateId",
    "retained", "removedRare", "removedNearConstant"))
  expect_equal(nrow(predictors), 7L)
  exported <- utils::read.csv(file.path(output, "export", "preprocessingSummary.csv"))
  expect_false(any(c("trainingRows", "nonzeroRows", "mean") %in% names(exported)))
})
