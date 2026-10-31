#' Sensitivity analyses
#'
#' The sensitivity analyses requested in review, each a set of [runStudy()]
#' arguments that differs from the main analysis in one respect:
#' - `noAtrialFibrillation`: bleeding after warfarin without the atrial
#'   fibrillation phenotype, which the cohort's entry criterion makes nearly
#'   constant. Every workflow, phenotype predictors only.
#' - `leadSecondLargest`, `leadSmallest`: ODAL and ADAP with the second-largest
#'   or the smallest training database as lead site. Phenotype predictors only,
#'   as ODAL is not fitted with condition predictors.
#' - `equalWeights`: coefficient averaging and the probability ensemble with
#'   equal instead of sample-size weights.
#' - `adapFromLocalModels`: ADAP started from the local Cyclops models of
#'   coefficient averaging instead of separate glmnet fits.
#' - `pooledRowTuning`: pooled lasso with its penalty selected by
#'   outcome-stratified three-fold CV over the pooled training rows.
#'
#' @return A named list; each element holds `featureSets` (the predictor sets
#'   it applies to) and the [runStudy()] arguments it sets.
#' @export
studySensitivityAnalyses <- function() {
  phenotypes <- "ageSexPhenotypes"
  both <- c("ageSexPhenotypes", "ageSexConditions")
  surrogates <- c("ODAL1", "ODAL2", "ADAP1", "ADAP2")
  list(
    noAtrialFibrillation = list(featureSets = "ageSexPhenotypesNoAtrialFibrillation",
      tasks = "BleedWarfarin", methods = studyMethods()),
    leadSecondLargest = list(featureSets = phenotypes, methods = surrogates, leadRank = 2L),
    leadSmallest = list(featureSets = phenotypes, methods = surrogates, leadRank = 4L),
    equalWeights = list(featureSets = both, methods = c("LocalAvgLasso", "LocalEnsembleLasso"),
      localWeighting = "equalClient"),
    adapFromLocalModels = list(featureSets = both, methods = c("ADAP1", "ADAP2"),
      adapInitBackend = "localCyclops"),
    pooledRowTuning = list(featureSets = both, methods = "PooledLasso", pooledTuning = "rows")
  )
}

#' Run one sensitivity analysis
#'
#' Runs a [studySensitivityAnalyses()] entry for one predictor set into its own
#' result directory, `<outputRoot>/<analysis>_<featureSet>`, with the held-out
#' databases appended when only some are run. Use one R session per result
#' directory; rerunning the same call resumes it.
#'
#' @param analysis Name of a [studySensitivityAnalyses()] entry.
#' @param dataRoot,siteIds,taskDirectories,hosts,mirai As for [runStudy()].
#' @param featureSet Predictor set; defaults to the analysis's first.
#' @param holdouts Databases to hold out; all by default.
#' @param outputRoot Directory under which the result directory is created.
#' @return The results table, invisibly.
#' @export
runSensitivityAnalysis <- function(analysis, dataRoot, siteIds, featureSet = NULL, holdouts = siteIds,
                                   outputRoot = "results/sensitivity", taskDirectories = NULL,
                                   hosts = "localhost", mirai = TRUE) {
  analyses <- studySensitivityAnalyses()
  if (length(analysis) != 1L || !analysis %in% names(analyses)) {
    stop("analysis must name one of: ", paste(names(analyses), collapse = ", "), call. = FALSE)
  }
  settings <- analyses[[analysis]]
  if (is.null(featureSet)) featureSet <- settings$featureSets[[1]]
  if (length(featureSet) != 1L || !featureSet %in% settings$featureSets) {
    stop(analysis, " applies to: ", paste(settings$featureSets, collapse = ", "), call. = FALSE)
  }
  if (!length(holdouts) || !all(holdouts %in% siteIds)) {
    stop("holdouts must name siteIds", call. = FALSE)
  }
  name <- paste(analysis, featureSet, sep = "_")
  if (!setequal(holdouts, siteIds)) name <- paste(c(name, holdouts), collapse = "_")
  arguments <- c(settings[setdiff(names(settings), "featureSets")],
    list(dataRoot = dataRoot, siteIds = siteIds, outputDirectory = file.path(outputRoot, name),
      featureSets = featureSet, holdouts = holdouts, taskDirectories = taskDirectories,
      hosts = hosts, mirai = mirai))
  do.call(runStudy, arguments)
}
