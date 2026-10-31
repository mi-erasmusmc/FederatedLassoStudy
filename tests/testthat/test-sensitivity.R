test_that("the atrial fibrillation sensitivity set is the phenotype set without that phenotype", {
  features <- studyFeatureSettings(c("ageSexPhenotypes", "ageSexPhenotypesNoAtrialFibrillation"))
  main <- unclass(features$ageSexPhenotypes)
  noAf <- unclass(features$ageSexPhenotypesNoAtrialFibrillation)
  expect_false("excludeCovariateIds" %in% names(main))
  expect_identical(noAf$excludeCovariateIds, 1160049)
  expect_identical(noAf[names(main)], main)
  covariateRef <- data.frame(covariateId = c(1002, 8532001, 1160049, 1165049), analysisId = c(2, 1, 49, 49),
    covariateName = letters[1:4])
  selected <- FederatedLearning::filterCovariateRef(covariateRef, features$ageSexPhenotypesNoAtrialFibrillation)
  expect_identical(selected$covariateId, c(1002, 8532001, 1165049))
})

test_that("sensitivity settings add options only where they differ from the main analysis", {
  options <- function(...) {
    settings <- studyMethodSettings(...)
    stats::setNames(lapply(settings, `[[`, "options"), vapply(settings, `[[`, "", "method"))
  }
  main <- options()
  for (method in c("PooledLasso", "LocalAvgLasso", "LocalEnsembleLasso")) expect_length(main[[method]], 0L)
  expect_false(any(c("leadRank", "localAverageWeighting", "pooledTuning") %in% unlist(lapply(main, names))))

  lead <- options(c("ODAL1", "ADAP2", "DualAvg"), leadRank = 4L)
  expect_identical(lead$ODAL1$leadRank, 4L)
  expect_identical(lead$ADAP2$leadRank, 4L)
  expect_null(lead$DualAvg$leadRank)
  expect_identical(lead$ODAL1[setdiff(names(lead$ODAL1), "leadRank")], main$ODAL1)

  equal <- options(c("LocalAvgLasso", "LocalEnsembleLasso"), localWeighting = "equalClient")
  expect_identical(equal$LocalAvgLasso, list(localAverageWeighting = "equalClient"))
  expect_identical(equal$LocalEnsembleLasso, list(localEnsembleWeighting = "equalClient"))
  expect_identical(options("PooledLasso", pooledTuning = "rows")$PooledLasso, list(pooledTuning = "rows"))
  expect_identical(options("ADAP1", adapInitBackend = "localCyclops")$ADAP1$adapInitBackend, "localCyclops")
  expect_error(studyMethodSettings(leadRank = 0), "leadRank")
})

test_that("each sensitivity analysis runs into its own result directory", {
  calls <- list()
  local_mocked_bindings(runStudy = function(...) {
    calls[[length(calls) + 1L]] <<- list(...)
    invisible(NULL)
  })
  sites <- c("CCAE", "JMDC", "MDCD", "MDCR", "OPTUMEHR")
  runSensitivityAnalysis("noAtrialFibrillation", "data", sites)
  runSensitivityAnalysis("noAtrialFibrillation", "data", sites, holdouts = c("JMDC", "CCAE"))
  runSensitivityAnalysis("pooledRowTuning", "data", sites, featureSet = "ageSexConditions")
  runSensitivityAnalysis("leadSmallest", "data", sites)

  expect_identical(vapply(calls, `[[`, "", "outputDirectory"), c(
    "results/sensitivity/noAtrialFibrillation_ageSexPhenotypesNoAtrialFibrillation",
    "results/sensitivity/noAtrialFibrillation_ageSexPhenotypesNoAtrialFibrillation_JMDC_CCAE",
    "results/sensitivity/pooledRowTuning_ageSexConditions",
    "results/sensitivity/leadSmallest_ageSexPhenotypes"))
  expect_identical(calls[[1]]$tasks, "BleedWarfarin")
  expect_identical(calls[[1]]$methods, studyMethods())
  expect_identical(calls[[2]]$holdouts, c("JMDC", "CCAE"))
  expect_identical(calls[[3]]$pooledTuning, "rows")
  expect_identical(calls[[3]]$methods, "PooledLasso")
  expect_null(calls[[3]]$tasks)
  expect_identical(calls[[4]]$leadRank, 4L)
  expect_true(all(vapply(calls, `[[`, logical(1), "mirai")))

  expect_error(runSensitivityAnalysis("leadSmallest", "data", sites, featureSet = "ageSexConditions"), "applies to")
  expect_error(runSensitivityAnalysis("unknown", "data", sites), "must name one of")
})
