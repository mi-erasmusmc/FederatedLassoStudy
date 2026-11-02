# Runs runProductionChecks() on the main analysis and the sensitivity analyses
# and bundles its aggregate tables for export. From the study's working
# directory (the one holding data/ and results/), after RunSensitivity.R has
# finished, in R:
#
#   source(system.file("extras", "RunChecks.R", package = "FederatedLassoStudy"))
#
# It reads saved models and site data and fits nothing but the local logistic
# fits with which ODAL starts. It writes results/productionChecks/export and
# checks.tar.gz. Set `sites` or `dataRoot` before sourcing to change them.

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

# The main analysis first: its pooled lasso is the reference wherever it exists.
directories <- c(Sys.glob("results/final/*"), Sys.glob("results/sensitivity/*"))
directories <- directories[file.exists(file.path(directories, "comparison_results.csv"))]
message(format(Sys.time(), "%H:%M"), "  Checking the saved models of ", length(directories), " result directories.")
started <- Sys.time()
runProductionChecks(directories, dataRoot, sites, "results/productionChecks")
bundle <- Sys.glob("results/productionChecks/export/*.csv")
utils::tar("checks.tar.gz", files = bundle, compression = "gzip", tar = "internal")
message(sprintf("\n%s  Wrote checks.tar.gz (%d files) after %.1f h.", format(Sys.time(), "%H:%M"), length(bundle),
  as.numeric(difftime(Sys.time(), started, units = "hours"))))
