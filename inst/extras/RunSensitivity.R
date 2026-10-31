# Runs every sensitivity analysis requested in review and bundles their
# aggregate results for export. From the study's working directory (the one
# holding data/ and results/), in R:
#
#   source(system.file("extras", "RunSensitivity.R", package = "FederatedLassoStudy"))
#
# Each run is a background Rscript process with its own log in
# results/sensitivity/logs. `workers` of them run at a time, longest first.
# Leave this session open: it reports each run as it finishes, then writes
# the coefficient and calibration exports and sensitivity.tar.gz. Sourcing it
# again resumes unfinished runs. Set `sites`, `dataRoot` or `workers` before
# sourcing to change them; the site order must match the main analysis, as it
# numbers the folds.

library(FederatedLassoStudy)
if (!exists("sites", inherits = FALSE)) sites <- c("MDCD", "MDCR", "CCAE", "JMDC", "OPTUMEHR")
if (!exists("dataRoot", inherits = FALSE)) dataRoot <- "data"
if (!exists("workers", inherits = FALSE)) workers <- 5L
outputRoot <- "results/sensitivity"
logDirectory <- file.path(outputRoot, "logs")
dir.create(logDirectory, recursive = TRUE, showWarnings = FALSE)

# The atrial fibrillation analysis runs one held-out database per process, as
# DualAvg took 3-10 hours per fold on this task in the main analysis.
jobs <- c(
  lapply(sites[c(1, 4, 3, 2, 5)], function(holdout) {
    list("noAtrialFibrillation", "ageSexPhenotypesNoAtrialFibrillation", holdout)
  }),
  list(list("leadSecondLargest", "ageSexPhenotypes", sites), list("leadSmallest", "ageSexPhenotypes", sites)),
  unlist(lapply(c("adapFromLocalModels", "pooledRowTuning", "equalWeights"), function(analysis) {
    lapply(c("ageSexConditions", "ageSexPhenotypes"), function(featureSet) list(analysis, featureSet, sites))
  }), recursive = FALSE))
# Named as runSensitivityAnalysis() names their result directories.
names(jobs) <- vapply(jobs, function(job) {
  name <- paste(job[[1]], job[[2]], sep = "_")
  if (length(job[[3]]) < length(sites)) name <- paste(c(name, job[[3]]), collapse = "_")
  name
}, "")

jobFile <- function(name, extension) file.path(logDirectory, paste0(name, extension))

launch <- function(name) {
  job <- jobs[[name]]
  unlink(jobFile(name, c(".pid", ".status")))
  call <- bquote(runSensitivityAnalysis(.(job[[1]]), .(dataRoot), .(sites), featureSet = .(job[[2]]),
    holdouts = .(job[[3]]), outputRoot = .(outputRoot)))
  writeLines(c(
    sprintf("writeLines(as.character(Sys.getpid()), %s)", deparse1(jobFile(name, ".pid"))),
    sprintf("log <- file(%s, open = \"at\")", deparse1(jobFile(name, ".log"))),
    "sink(log)",
    "sink(log, type = \"message\")",
    "status <- tryCatch({",
    "  library(FederatedLassoStudy)",
    paste0("  results <- ", deparse1(call)),
    "  errors <- if (\"error\" %in% names(results)) sum(!is.na(results$error) & nzchar(results$error)) else 0L",
    "  sprintf(\"ok: %d rows, %d with errors\", nrow(results), errors)",
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

message(sprintf("%s  Starting %d runs, %d at a time. Logs: %s", format(Sys.time(), "%H:%M"), length(jobs),
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

message(format(Sys.time(), "%H:%M"), "  All runs finished. Exporting coefficients and calibration curves.")
directories <- file.path(outputRoot, names(jobs))
directories <- directories[file.exists(file.path(directories, "comparison_results.csv"))]
try(runCoefficientAgreement(grep("/noAtrialFibrillation_", directories, value = TRUE),
  "results/sensitivityCoefficients"))
try(runCalibrationCurveExport(directories, "results/sensitivityCalibration"))
bundle <- Sys.glob(c(outer(directories, c("comparison_results.csv", "run_manifest.csv",
  "preprocessing/summary.csv"), file.path), "results/sensitivity*/export/*.csv", file.path(logDirectory, "status.csv")))
# Some paths exceed the old 100-byte tar limit; R and current tar programs read them.
withCallingHandlers(utils::tar("sensitivity.tar.gz", files = bundle, compression = "gzip", tar = "internal"),
  warning = function(w) if (grepl("not portable", conditionMessage(w))) invokeRestart("muffleWarning"))

failed <- statuses$status[!startsWith(statuses$status, "ok: ") | !grepl(" 0 with errors", statuses$status)]
message(sprintf("\nWrote sensitivity.tar.gz (%d files). %d of %d runs ended without errors.", length(bundle),
  nrow(statuses) - length(failed), length(jobs)))
print(statuses, row.names = FALSE)
