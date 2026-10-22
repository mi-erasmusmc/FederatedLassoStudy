#' Audit the shared feature filter
#'
#' Builds the shared training-only feature filter for every task and held-out
#' database from site summaries, without fitting models or loading held-out
#' data. Each site's saved population settings must match the task's.
#'
#' The `preprocessing` directory holds per-predictor training counts and stays
#' private. The `export` directory holds only predictor counts and covariate
#' ids, for review before export:
#' - `preprocessingSummary.csv`: one row per task and held-out database with
#'   candidate, retained and removed predictors, and whether the any-antibiotic
#'   indicator is retained (`NA` when no training site has it, as in the
#'   wider specifications);
#' - `predictors.csv`: each candidate covariate id and why it was removed.
#'
#' @inheritParams runStudy
#' @param outputDirectory A new result directory.
#' @param tasks Names of [studyTasks()] to audit.
#' @param featureSets Names of [studyFeatureSettings()] matching the extraction.
#' @return The export summary, invisibly.
#' @export
runPreprocessingAudit <- function(dataRoot, siteIds, outputDirectory,
                                  tasks = names(studyTasks()), featureSets = "ageSexPhenotypes",
                                  taskDirectories = NULL, hosts = "localhost") {
  sites <- studySites(dataRoot, siteIds, tasks, taskDirectories)
  featureSettings <- studyFeatureSettings(featureSets)
  for (task in tasks) {
    FederatedLearning::runComparison(
      taskSettings = studyTasks()[[task]],
      sites = sites[[task]],
      featureSettings = featureSettings,
      methodSettings = FederatedLearning::createMethodSettings("PooledLasso"),
      outputDirectory = outputDirectory,
      runSettings = FederatedLearning::createRunSettings(hosts = hosts, auditOnly = TRUE)
    )
  }
  invisible(writePreprocessingExport(outputDirectory, siteIds, tasks))
}

# Writes the exportable audit tables: predictor ids and counts, no patient counts.
writePreprocessingExport <- function(outputDirectory, siteIds, tasks) {
  directory <- file.path(outputDirectory, "preprocessing")
  summary <- utils::read.csv(file.path(directory, "summary.csv"), stringsAsFactors = FALSE)
  summary <- summary[summary$task %in% tasks, , drop = FALSE]
  predictors <- lapply(seq_len(nrow(summary)), function(i) {
    row <- summary[i, ]
    stem <- paste(row$task, paste0("fold", row$fold), row$featureSet, sep = "_")
    audit <- readRDS(file.path(directory, paste0(stem, ".rds")))$audit
    data.frame(task = row$task, heldOutSite = siteIds[[row$fold]], featureSet = row$featureSet,
      covariateId = audit$covariateId, retained = audit$retained,
      removedRare = audit$removedRare, removedNearConstant = audit$removedNearConstant,
      stringsAsFactors = FALSE)
  })
  predictors <- do.call(rbind, predictors)
  antibiotic <- predictors[predictors$covariateId == format(anyAntibioticCovariateId, scientific = FALSE), ]
  key <- function(x) paste(x$task, x$heldOutSite, x$featureSet)
  summary$heldOutSite <- siteIds[summary$fold]
  exported <- data.frame(task = summary$task, heldOutSite = summary$heldOutSite,
    trainingSites = vapply(summary$fold, function(fold) paste(siteIds[-fold], collapse = ";"), ""),
    featureSet = summary$featureSet, candidatePredictors = summary$candidatePredictors,
    retainedPredictors = summary$retainedPredictors, removedRare = summary$removedRare,
    removedNearConstant = summary$removedNearConstant,
    antibioticRetained = antibiotic$retained[match(key(summary), key(antibiotic))],
    fingerprint = summary$fingerprint, stringsAsFactors = FALSE)
  exportDirectory <- file.path(outputDirectory, "export")
  dir.create(exportDirectory, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(exported, file.path(exportDirectory, "preprocessingSummary.csv"), row.names = FALSE)
  utils::write.csv(predictors, file.path(exportDirectory, "predictors.csv"), row.names = FALSE)
  exported
}
