#' Run the study comparison
#'
#' Fits each method on all but one database and evaluates it on the held-out
#' database, for every task and every database in `holdouts`. Results and model
#' artifacts are written to `outputDirectory`; completed combinations there
#' are reused.
#'
#' @param dataRoot Directory holding the extracted data as
#'   `<task><extraction>/<site>`, where the extraction suffix is empty for
#'   `ageSexPhenotypes`, `AgeSexConditions` for `ageSexConditions` and `Default`
#'   for `defaultCovariates`.
#' @param siteIds Database names: the site directory names under each task.
#' @param outputDirectory Result directory.
#' @param tasks Names of [studyTasks()] to run.
#' @param holdouts Databases to hold out, one outer fold each.
#' @param methods Methods from [studyMethods()].
#' @param featureSets Names of [studyFeatureSettings()] matching the extraction.
#' @param dualAvgEpsilon,dualAvgLambda,dualAvgRounds Passed to [studyMethodSettings()].
#' @param taskDirectories Task directory names under `dataRoot`, named by task,
#'   where they differ from that layout.
#' @param hosts Worker hosts, one per site or a single host for all sites.
#' @param mirai Use mirai workers instead of socket (PSOCK) workers. Results are
#'   the same; socket workers stall on each message of a few kilobytes or more,
#'   which makes DualAvg about five times slower at thousands of predictors.
#' @return The results table, invisibly.
#' @export
runStudy <- function(dataRoot, siteIds, outputDirectory, tasks = names(studyTasks()),
                     holdouts = siteIds, methods = studyMethods(),
                     featureSets = "ageSexPhenotypes", dualAvgEpsilon = 1e-10,
                     dualAvgLambda = NULL, dualAvgRounds = 10000L, taskDirectories = NULL,
                     hosts = "localhost", mirai = FALSE) {
  sites <- studySites(dataRoot, siteIds, tasks, featureSets, taskDirectories)
  methodSettings <- studyMethodSettings(methods, dualAvgEpsilon, dualAvgLambda, dualAvgRounds)
  featureSettings <- studyFeatureSettings(featureSets)
  for (task in tasks) {
    FederatedLearning::runComparison(
      taskSettings = studyTasks()[[task]],
      sites = sites[[task]],
      featureSettings = featureSettings,
      methodSettings = methodSettings,
      outputDirectory = outputDirectory,
      holdouts = holdouts,
      runSettings = FederatedLearning::createRunSettings(hosts = hosts, mirai = mirai)
    )
  }
  invisible(FederatedLearning::readComparisonResults(outputDirectory))
}

# Site data directories for each task, named by site.
studySites <- function(dataRoot, siteIds, tasks, featureSets, taskDirectories = NULL) {
  if (!length(tasks) || !all(tasks %in% names(studyTasks()))) {
    stop("tasks must name studyTasks(): ", paste(names(studyTasks()), collapse = ", "), call. = FALSE)
  }
  if (!is.character(siteIds) || length(siteIds) < 2L || anyDuplicated(siteIds)) {
    stop("siteIds must name at least two distinct sites", call. = FALSE)
  }
  suffix <- unique(extractionSuffixes[names(studyFeatureSettings(featureSets))])
  if (length(suffix) != 1L) stop("featureSets must come from one extraction", call. = FALSE)
  directories <- stats::setNames(paste0(tasks, suffix), tasks)
  override <- intersect(names(taskDirectories), tasks)
  directories[override] <- taskDirectories[override]
  lapply(directories, function(directory) {
    stats::setNames(file.path(dataRoot, directory, siteIds), siteIds)
  })
}

# Task directory suffix of each predictor specification's extraction.
extractionSuffixes <- c(ageSexPhenotypes = "", ageSexConditions = "AgeSexConditions",
  defaultCovariates = "Default")
