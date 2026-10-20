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

test_that("predictors are age, one sex indicator, phenotypes and any antibiotic", {
  features <- studyFeatureSettings()
  expect_identical(names(features), "ageSexPhenotypes")
  combination <- features$ageSexPhenotypes$combinations[[1]]
  expect_identical(combination$covariateId, 999001049)
  expect_identical(combination$components, c(1201049, 1202049, 1203049, 1204049, 1205049, 1206049,
    1207049, 1208049, 1209049, 1210049, 1211049, 1212049, 1213049, 1214049))
  covariateRef <- data.frame(covariateId = c(1002, 8507001, 8532001, 1152049, 1201049, 1214049, 1215049, 4000102),
    analysisId = c(2, 1, 1, 49, 49, 49, 49, 102), covariateName = letters[1:8])
  selected <- FederatedLearning::filterCovariateRef(covariateRef, features$ageSexPhenotypes)
  expect_identical(selected$covariateId, c(1002, 8532001, 1152049, 1215049, 999001049))
  expect_identical(features$ageSexPhenotypes$minFraction, 0.001)
})
