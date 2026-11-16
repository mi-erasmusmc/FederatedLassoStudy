# Refits DualAvg for every main-analysis evaluation with centred predictors, scores the new models
# against the converged pooled lasso and bundles the aggregate results for export. From the
# study's working directory (the one holding data/ and results/), in R:
#
#   source(system.file("extras", "RunCentredDualAvg.R", package = "FederatedLassoStudy"))
#
# DualAvg runs as in the main analysis (penalty tuned by leaving out each training database,
# at most 10,000 rounds, mirai workers), except that age is centred at its pooled training mean
# and binary predictors present in most training entries are coded the other way round. That
# leaves the objective unchanged and removes the near-collinearity with the intercept along
# which the main-run fits stopped short. It needs results/convergedCalibrated (from
# RunConvergedReference.R) for the reference coefficients.
#
# Each held-out database is one background Rscript process with its own log in
# results/dualAvgCentred/logs, `workers` at a time (default 2, about eight busy processors),
# warfarin first, since its folds took 3-10 hours each in the main analysis. Leave this session
# open: it reports each run as it finishes, then scores the models and writes
# dualAvgCentred.tar.gz. Sourcing it again resumes unfinished runs and refits fits that ended in
# an error. Set `sites`, `dataRoot` or `workers` before sourcing to change them; the site order
# must match the main analysis, as it numbers the folds.

library(FederatedLassoStudy)
local({
  pin <- sub(".*@", "", trimws(utils::packageDescription("FederatedLassoStudy")$Remotes))
  installed <- utils::packageDescription("FederatedLearning")$RemoteSha
  if (is.null(installed) || !identical(substr(installed, 1, 7), substr(pin, 1, 7))) {
    stop("FederatedLearning ", if (is.null(installed)) "was not installed from GitHub" else
      paste("is at", substr(installed, 1, 7)), "; this study needs ", substr(pin, 1, 7), ". Run: ",
      "remotes::install_github(\"mi-erasmusmc/FederatedLearning@", substr(pin, 1, 7), "\", upgrade = \"never\")",
      call. = FALSE)
  }
})
if (!exists("sites", inherits = FALSE)) sites <- c("MDCD", "MDCR", "CCAE", "JMDC", "OPTUMEHR")
if (!exists("dataRoot", inherits = FALSE)) dataRoot <- "data"
if (!exists("workers", inherits = FALSE)) workers <- 2L
referenceDirectory <- "results/convergedCalibrated"
if (!file.exists(file.path(referenceDirectory, "export", "referenceCoefficients.csv"))) {
  stop("Run RunConvergedReference.R first: ", referenceDirectory, "/export/referenceCoefficients.csv is missing",
    call. = FALSE)
}
outputRoot <- "results/dualAvgCentred"
logDirectory <- file.path(outputRoot, "logs")
dir.create(logDirectory, recursive = TRUE, showWarnings = FALSE)

grid <- expand.grid(holdout = sites, featureSet = c("ageSexPhenotypes", "ageSexConditions"),
  task = c("BleedWarfarin", "DvtTetraplegia", "BleedTetraPlegia"), stringsAsFactors = FALSE)
jobs <- lapply(seq_len(nrow(grid)), function(i) as.list(grid[i, c("task", "featureSet", "holdout")]))
names(jobs) <- paste(grid$task, grid$featureSet, grid$holdout, sep = "_")

jobFile <- function(name, extension) file.path(logDirectory, paste0(name, extension))

launch <- function(name) {
  job <- jobs[[name]]
  unlink(jobFile(name, c(".pid", ".status")))
  call <- bquote(runStudy(.(dataRoot), .(sites), .(file.path(outputRoot, name)), tasks = .(job$task),
    holdouts = .(job$holdout), methods = "DualAvg", featureSets = .(job$featureSet), dualAvgCenter = TRUE,
    mirai = TRUE, rerunErrors = TRUE))
  writeLines(c(
    sprintf("writeLines(as.character(Sys.getpid()), %s)", deparse1(jobFile(name, ".pid"))),
    sprintf("log <- file(%s, open = \"at\")", deparse1(jobFile(name, ".log"))),
    "sink(log)",
    "sink(log, type = \"message\")",
    "status <- tryCatch({",
    "  library(FederatedLassoStudy)",
    paste0("  results <- ", deparse1(call)),
    "  failed <- if (\"error\" %in% names(results)) !is.na(results$error) & nzchar(results$error) else FALSE",
    "  sprintf(\"ok: %d fits, %d with errors\", nrow(results), sum(failed))",
    "}, error = function(e) paste(\"failed:\", conditionMessage(e)))",
    "message(status)",
    sprintf("writeLines(status, %s)", deparse1(jobFile(name, ".status")))
  ), jobFile(name, ".R"))
  system2(file.path(R.home("bin"), "Rscript"), shQuote(jobFile(name, ".R")), stdout = NULL, stderr = NULL,
    wait = FALSE)
  Sys.time()
}

alive <- function(pid) {
  if (.Platform$OS.type == "windows") {
    tasks <- system2("tasklist", c("/FI", shQuote(paste("PID eq", pid)), "/NH"), stdout = TRUE)
    any(grepl(paste0("\\b", pid, "\\b"), tasks))
  } else {
    isTRUE(tools::pskill(pid, 0L))
  }
}

# The run's status line, or NULL while it is still running.
finished <- function(name, started) {
  if (file.exists(jobFile(name, ".status"))) return(readLines(jobFile(name, ".status"))[[1]])
  if (!file.exists(jobFile(name, ".pid"))) {
    return(if (difftime(Sys.time(), started, units = "mins") > 10) "failed: did not start" else NULL)
  }
  if (alive(as.integer(readLines(jobFile(name, ".pid"))[[1]]))) return(NULL)
  if (file.exists(jobFile(name, ".status"))) readLines(jobFile(name, ".status"))[[1]] else "failed: stopped, see log"
}

message(sprintf("%s  Starting %d DualAvg runs, %d at a time. Logs: %s", format(Sys.time(), "%H:%M"), length(jobs),
  workers, logDirectory))
queue <- names(jobs)
running <- list()
statuses <- data.frame(run = character(), status = character(), hours = numeric())
repeat {
  for (name in names(running)) {
    status <- finished(name, running[[name]])
    if (is.null(status)) next
    hours <- as.numeric(difftime(Sys.time(), running[[name]], units = "hours"))
    statuses[nrow(statuses) + 1L, ] <- list(name, status, round(hours, 2))
    running[[name]] <- NULL
    message(sprintf("%s  %d/%d %s: %s (%.1f h)", format(Sys.time(), "%H:%M"), nrow(statuses), length(jobs),
      name, status, hours))
  }
  while (length(running) < workers && length(queue)) {
    running[[queue[[1]]]] <- launch(queue[[1]])
    queue <- queue[-1]
  }
  if (!length(running)) break
  Sys.sleep(60)
}
utils::write.csv(statuses, file.path(logDirectory, "status.csv"), row.names = FALSE)

message(format(Sys.time(), "%H:%M"), "  All runs finished. Scoring against the converged pooled lasso.")
directories <- file.path(outputRoot, names(jobs))
directories <- directories[file.exists(file.path(directories, "comparison_results.csv"))]
scores <- runReferenceScoring(dataRoot, sites, directories, referenceDirectory, "results/dualAvgCentredScores")
bundle <- Sys.glob(c(outer(directories, c("comparison_results.csv", "run_manifest.csv"), file.path),
  "results/dualAvgCentredScores/export/*.csv", file.path(logDirectory, "status.csv")))
withCallingHandlers(utils::tar("dualAvgCentred.tar.gz", files = bundle, compression = "gzip", tar = "internal"),
  warning = function(w) if (grepl("not portable", conditionMessage(w))) invokeRestart("muffleWarning"))

results <- do.call(rbind, lapply(directories, FederatedLearning::readComparisonResults))
fits <- results[results$method == "DualAvg", , drop = FALSE]
failed <- if ("error" %in% names(fits)) !is.na(fits$error) & nzchar(fits$error) else rep(FALSE, nrow(fits))
dual <- scores$referenceScores[scores$referenceScores$model == "DualAvg", , drop = FALSE]
message(sprintf(paste0("\nWrote dualAvgCentred.tar.gz (%d files). %d of 30 DualAvg fits finished without errors; ",
  "%d ran all 10,000 rounds. AUC difference from the converged pooled lasso: median %+.4f, range %+.3f to %+.3f ",
  "(the main run's DualAvg: -0.010 to +0.023)."), length(bundle), sum(!failed),
  sum(fits$roundsCompleted[!failed] >= 10000, na.rm = TRUE), stats::median(dual$aucDifference),
  min(dual$aucDifference), max(dual$aucDifference)))
print(statuses, row.names = FALSE)
