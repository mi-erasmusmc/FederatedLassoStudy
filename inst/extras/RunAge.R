# Runs runAgeAnalysis() on the main analysis (all tasks, both predictor sets,
# phenotype predictors first) and bundles its aggregate tables for export.
# From the study's working directory (the one holding data/ and results/), in R:
#
#   source(system.file("extras", "RunAge.R", package = "FederatedLassoStudy"))
#
# It reads the saved pooled lasso and DualAvg models and the site data, refits
# logistic models on the pooled training data of each evaluation, and writes
# results/ageAnalysis/export and age.tar.gz. The export tables are rewritten
# after each task and predictor set; if the run stops, sourcing the script
# again continues from there. Set `sites`, `dataRoot` or `threads` before
# sourcing to change them.

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
if (!exists("threads", inherits = FALSE)) threads <- 4L

started <- Sys.time()
tables <- runAgeAnalysis(dataRoot, sites, Sys.glob("results/final/*"), "results/ageAnalysis", threads = threads)
bundle <- Sys.glob("results/ageAnalysis/export/*.csv")
utils::tar("age.tar.gz", files = bundle, compression = "gzip", tar = "internal")
fits <- tables$ageFits
tight <- fits[fits$fit %in% c("tightPooled", "tightDualAvgPenalty"), ]
message(sprintf(paste0("\n%s  Wrote age.tar.gz (%d files) after %.1f h. %d of %d fits failed. ",
  "Largest optimality violation of the tight lasso refits %.1g (should be below 1e-6)."),
  format(Sys.time(), "%H:%M"), length(bundle), as.numeric(difftime(Sys.time(), started, units = "hours")),
  sum(!is.na(fits$error)), nrow(fits), max(tight$kktMax, na.rm = TRUE)))
