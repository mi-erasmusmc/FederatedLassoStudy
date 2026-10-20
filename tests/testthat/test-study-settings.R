test_that("task horizons and cohorts follow the protocol", {
  tasks <- studyTasks()
  expect_identical(names(tasks), c("DvtTetraplegia", "BleedTetraPlegia", "BleedWarfarin"))
  horizons <- vapply(tasks, function(t) t$populationSettings$riskWindowEnd, numeric(1))
  expect_equal(horizons, c(DvtTetraplegia = 365, BleedTetraPlegia = 365, BleedWarfarin = 30))
  expect_equal(vapply(tasks, function(t) t$populationSettings$riskWindowStart, numeric(1)),
    c(DvtTetraplegia = 1, BleedTetraPlegia = 1, BleedWarfarin = 1))
  expect_equal(vapply(tasks, `[[`, numeric(1), "targetId"),
    c(DvtTetraplegia = 15026, BleedTetraPlegia = 15026, BleedWarfarin = 13524))
  expect_equal(vapply(tasks, `[[`, numeric(1), "outcomeId"),
    c(DvtTetraplegia = 19047, BleedTetraPlegia = 19012, BleedWarfarin = 19012))
  expect_identical(vapply(tasks, `[[`, "", "name"), stats::setNames(names(tasks), names(tasks)))
  for (task in tasks) expect_s3_class(task, "taskSettings")
})

test_that("population settings match the production extraction", {
  population <- unclass(studyPopulationSettings(365L))
  expected <- list(binary = TRUE, includeAllOutcomes = TRUE, firstExposureOnly = FALSE,
    washoutPeriod = 0, removeSubjectsWithPriorOutcome = TRUE, priorOutcomeLookback = 99999L,
    requireTimeAtRisk = FALSE, minTimeAtRisk = 1L, riskWindowStart = 1L,
    startAnchor = "cohort start", riskWindowEnd = 365L, endAnchor = "cohort start",
    restrictTarToCohortEnd = FALSE)
  expect_equal(population[names(expected)], expected)
})

test_that("predictors are age, one sex indicator and phenotypes", {
  features <- studyFeatureSettings()
  expect_identical(names(features), "ageSexPhenotypes")
  covariateRef <- data.frame(covariateId = c(1002, 8507001, 8532001, 1152049, 4000102),
    analysisId = c(2, 1, 1, 49, 102))
  selected <- FederatedLearning::filterCovariateRef(covariateRef, features$ageSexPhenotypes)
  expect_identical(selected$covariateId, c(1002, 8532001, 1152049))
  expect_identical(features$ageSexPhenotypes$minFraction, 0.001)
})
