# Runs runPseudoSiteAnalysis() on the main analysis (phenotype predictors, all
# tasks and held-out databases, 100, 30 and 10 outcomes per pseudo-site) and
# bundles its aggregate tables for export. From the study's working directory
# (the one holding data/ and results/), after RunChecks.R has finished, in R:
#
#   source(system.file("extras", "RunPseudoSites.R", package = "FederatedLassoStudy"))
#
# It fits local lasso models at each pseudo-site, on `cores` processes, and
# writes results/pseudoSites/export and pseudoSites.tar.gz. Set `sites`,
# `dataRoot` or `cores` before sourcing to change them.

library(FederatedLassoStudy)
if (!exists("sites", inherits = FALSE)) sites <- c("MDCD", "MDCR", "CCAE", "JMDC", "OPTUMEHR")
if (!exists("dataRoot", inherits = FALSE)) dataRoot <- "data"
if (!exists("cores", inherits = FALSE)) cores <- 4L

started <- Sys.time()
tables <- runPseudoSiteAnalysis(dataRoot, sites, Sys.glob("results/final/*"), "results/pseudoSites", cores = cores)
bundle <- Sys.glob("results/pseudoSites/export/*.csv")
utils::tar("pseudoSites.tar.gz", files = bundle, compression = "gzip", tar = "internal")
checks <- tables$pseudoSites$check
message(sprintf("\n%s  Wrote pseudoSites.tar.gz (%d files) after %.1f h. Largest check difference %.2g (should be near 0).",
  format(Sys.time(), "%H:%M"), length(bundle), as.numeric(difftime(Sys.time(), started, units = "hours")),
  max(checks, na.rm = TRUE)))
