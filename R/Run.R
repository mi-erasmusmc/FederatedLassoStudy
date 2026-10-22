#' Run the study comparison
#'
#' Fits each method on all but one database and evaluates it on the held-out
#' database, for every task and every database in `holdouts`. Results and model
#' artifacts are written to `outputDirectory`; completed combinations there
#' are reused.
#'
#' @param dataRoot Directory holding the extracted data as `<task>/<site>`.
#' @param siteIds Database names: the site directory names under each task.
#' @param outputDirectory Result directory.
#' @param tasks Names of [studyTasks()] to run.
#' @param holdouts Databases to hold out, one outer fold each.
#' @param methods Methods from [studyMethods()].
#' @param featureSets Names of [studyFeatureSettings()] matching the extraction.
#' @param dualAvgEpsilon,dualAvgLambda Passed to [studyMethodSettings()].
#' @param taskDirectories Task directory names under `dataRoot`, named by task,
#'   where they differ from the task names.
#' @param hosts Worker hosts, one per site or a single host for all sites.
#' @return The results table, invisibly.
#' @export
runStudy <- function(dataRoot, siteIds, outputDirectory, tasks = names(studyTasks()),
                     holdouts = siteIds, methods = studyMethods(),
                     featureSets = "ageSexPhenotypes", dualAvgEpsilon = 1e-10,
                     dualAvgLambda = NULL, taskDirectories = NULL, hosts = "localhost") {
  sites <- studySites(dataRoot, siteIds, tasks, taskDirectories)
  methodSettings <- studyMethodSettings(methods, dualAvgEpsilon, dualAvgLambda)
  featureSettings <- studyFeatureSettings(featureSets)
  for (task in tasks) {
    FederatedLearning::runComparison(
      taskSettings = studyTasks()[[task]],
      sites = sites[[task]],
      featureSettings = featureSettings,
      methodSettings = methodSettings,
      outputDirectory = outputDirectory,
      holdouts = holdouts,
      runSettings = FederatedLearning::createRunSettings(hosts = hosts)
    )
  }
  invisible(FederatedLearning::readComparisonResults(outputDirectory))
}

# Site data directories for each task, named by site.
studySites <- function(dataRoot, siteIds, tasks, taskDirectories = NULL) {
  if (!length(tasks) || !all(tasks %in% names(studyTasks()))) {
    stop("tasks must name studyTasks(): ", paste(names(studyTasks()), collapse = ", "), call. = FALSE)
  }
  if (!is.character(siteIds) || length(siteIds) < 2L || anyDuplicated(siteIds)) {
    stop("siteIds must name at least two distinct sites", call. = FALSE)
  }
  directories <- stats::setNames(tasks, tasks)
  override <- intersect(names(taskDirectories), tasks)
  directories[override] <- taskDirectories[override]
  lapply(directories, function(directory) {
    stats::setNames(file.path(dataRoot, directory, siteIds), siteIds)
  })
}
