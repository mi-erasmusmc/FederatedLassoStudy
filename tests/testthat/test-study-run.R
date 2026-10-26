test_that("feature specifications are named and validated", {
  expect_identical(names(studyFeatureSettings()), "ageSexPhenotypes")
  all <- studyFeatureSettings(c("ageSexPhenotypes", "ageSexConditions", "defaultCovariates"))
  expect_identical(all$ageSexConditions$analysisIds, 102)
  expect_identical(all$ageSexConditions$covariateIds, c(1002, 8532001))
  expect_false(1 %in% all$defaultCovariates$analysisIds)
  expect_identical(all$defaultCovariates$covariateIds, 8532001)
  expect_error(studyFeatureSettings("wide"), "featureSets must name")
})

test_that("method settings carry the previous run's non-default options", {
  settings <- studyMethodSettings()
  expect_identical(vapply(settings, `[[`, "", "method"), studyMethods())
  options <- stats::setNames(lapply(settings, `[[`, "options"), studyMethods())
  expect_identical(options$DualAvg, list(k = 1L, etaClient = 1, etaServer = 1, epsilon = 1e-10,
    dualavgRounds = 10000L))
  expect_identical(options$ODAL2$odalInit, "ridgeFallback")
  expect_identical(options$ADAP2, list(epsilon = 1e-10))
  expect_identical(options$PooledLasso, list())
  fixed <- studyMethodSettings("DualAvg", dualAvgEpsilon = 1e-8, dualAvgLambda = 0.02)[[1]]$options
  expect_identical(fixed[c("epsilon", "dualavgTuneLambda", "dualavgLambda")],
    list(epsilon = 1e-8, dualavgTuneLambda = FALSE, dualavgLambda = 0.02))
  expect_identical(studyMethodSettings("DualAvg", dualAvgRounds = 5e4)[[1]]$options$dualavgRounds, 50000L)
  expect_error(studyMethodSettings("FedAvg"), "studyMethods")
})

test_that("runStudy runs each task on its own site directories", {
  calls <- list()
  local_mocked_bindings(
    runComparison = function(taskSettings, sites, featureSettings, methodSettings,
                             outputDirectory, holdouts, runSettings) {
      calls[[length(calls) + 1L]] <<- list(task = taskSettings$name, sites = sites, holdouts = holdouts,
        methods = vapply(methodSettings, `[[`, "", "method"), features = names(featureSettings),
        mirai = runSettings$mirai)
      invisible(NULL)
    },
    readComparisonResults = function(outputDirectory) data.frame(),
    .package = "FederatedLearning"
  )
  runStudy("root", c("A", "B", "C"), "out", tasks = c("DvtTetraplegia", "BleedWarfarin"),
    holdouts = "B", methods = c("PooledLasso", "DualAvg"), featureSets = "ageSexConditions",
    taskDirectories = c(BleedWarfarin = "warfarin"), mirai = TRUE)
  expect_identical(vapply(calls, `[[`, "", "task"), c("DvtTetraplegia", "BleedWarfarin"))
  expect_true(calls[[1]]$mirai)
  expect_identical(unname(calls[[1]]$sites), file.path("root", "DvtTetraplegiaAgeSexConditions", c("A", "B", "C")))
  expect_identical(unname(calls[[2]]$sites), file.path("root", "warfarin", c("A", "B", "C")))
  expect_identical(calls[[1]]$holdouts, "B")
  expect_identical(calls[[1]]$methods, c("PooledLasso", "DualAvg"))
  expect_identical(calls[[1]]$features, "ageSexConditions")
  default <- studySites("root", c("A", "B"), "BleedWarfarin", "defaultCovariates")
  expect_identical(unname(default$BleedWarfarin), file.path("root", "BleedWarfarinDefault", c("A", "B")))
  expect_identical(unname(studySites("root", c("A", "B"), "BleedWarfarin", "ageSexPhenotypes")[[1]]),
    file.path("root", "BleedWarfarin", c("A", "B")))
  expect_error(studySites("root", c("A", "B"), "BleedWarfarin", c("ageSexPhenotypes", "defaultCovariates")),
    "one extraction")
})

test_that("the threshold pilot refits at the tuned penalty without the outer database", {
  calls <- list()
  local_mocked_bindings(
    readComparisonResults = function(outputDirectory) data.frame(task = "DvtTetraplegia",
      method = c("PooledLasso", "DualAvg"), featureSet = "ageSexPhenotypes", clientId = "JMDC",
      selectedLambda = c(0.01, 0.0042)),
    .package = "FederatedLearning"
  )
  local_mocked_bindings(
    runStudy = function(dataRoot, siteIds, outputDirectory, ...) {
      calls[[length(calls) + 1L]] <<- c(list(sites = siteIds, output = outputDirectory), list(...))
      invisible(NULL)
    },
    writeThresholdPilotExport = function(...) data.frame()
  )
  sites <- c("MDCD", "MDCR", "CCAE", "JMDC", "OPTUMEHR")
  runDualAvgThresholdPilot("root", sites, "tuning", "out", innerHoldouts = c("MDCD", "CCAE"))
  expect_length(calls, 2L)
  expect_identical(calls[[1]]$sites, c("MDCD", "MDCR", "CCAE", "OPTUMEHR"))
  expect_identical(calls[[1]]$holdouts, c("MDCD", "CCAE"))
  expect_identical(vapply(calls, `[[`, numeric(1), "dualAvgEpsilon"), c(1e-10, 1e-8))
  expect_identical(vapply(calls, `[[`, numeric(1), "dualAvgLambda"), c(0.0042, 0.0042))
  expect_identical(vapply(calls, `[[`, "", "output"), file.path("out", c("epsilon1e-10", "epsilon1e-08")))
  expect_error(runDualAvgThresholdPilot("root", sites, "tuning", "out", innerHoldouts = "JMDC"),
    "not the outer held-out database")
  expect_error(runDualAvgThresholdPilot("root", sites, "tuning", "out", holdout = "MDCD"),
    "No tuned DualAvg penalty")
})

test_that("the rounds check refits each task's fold at its tuned penalty with more rounds", {
  calls <- list()
  local_mocked_bindings(
    readComparisonResults = function(outputDirectory) data.frame(task = c("DvtTetraplegia", "BleedWarfarin"),
      method = "DualAvg", featureSet = "ageSexPhenotypes", clientId = "JMDC", selectedLambda = c(4e-4, 7e-4)),
    .package = "FederatedLearning"
  )
  local_mocked_bindings(
    runStudy = function(dataRoot, siteIds, outputDirectory, ...) {
      calls[[length(calls) + 1L]] <<- c(list(sites = siteIds, output = outputDirectory), list(...))
      invisible(NULL)
    },
    writeRoundsCheckExport = function(...) data.frame()
  )
  sites <- c("MDCD", "MDCR", "CCAE", "JMDC", "OPTUMEHR")
  runDualAvgRoundsCheck("root", sites, "final", "out", tasks = c("DvtTetraplegia", "BleedWarfarin"))
  expect_identical(vapply(calls, `[[`, "", "tasks"), c("DvtTetraplegia", "BleedWarfarin"))
  expect_identical(calls[[1]]$sites, sites)
  expect_identical(calls[[1]]$holdouts, "JMDC")
  expect_identical(vapply(calls, `[[`, numeric(1), "dualAvgLambda"), c(4e-4, 7e-4))
  expect_identical(vapply(calls, `[[`, integer(1), "dualAvgRounds"), c(50000L, 50000L))
  expect_error(runDualAvgRoundsCheck("root", sites, "final", "out", tasks = "BleedTetraplegia"),
    "No tuned DualAvg penalty")
})
