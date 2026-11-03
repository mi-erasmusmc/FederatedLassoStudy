# Runs the production checks and the pseudo-site analysis for the two
# tetraplegia tasks only, which need little memory, while the sensitivity
# analyses are still running. The warfarin parts follow with RunChecks.R and
# RunPseudoSites.R once the server is free. From the study's working
# directory, after putting the library holding this package version first on
# .libPaths(), in R:
#
#   source(system.file("extras", "RunEarlyChecks.R", package = "FederatedLassoStudy"))
#
# Writes results/productionChecksTetraplegia, results/pseudoSitesTetraplegia
# and early.tar.gz. Set `sites`, `dataRoot`, `cores` or `siteLabelDirectory`
# before sourcing to change them.

library(FederatedLassoStudy)
if (!exists("sites", inherits = FALSE)) sites <- c("MDCD", "MDCR", "CCAE", "JMDC", "OPTUMEHR")
if (!exists("dataRoot", inherits = FALSE)) dataRoot <- "data"
if (!exists("cores", inherits = FALSE)) cores <- 2L
if (!exists("siteLabelDirectory", inherits = FALSE)) siteLabelDirectory <- "private/siteLabels"
if (!dir.exists(siteLabelDirectory)) {
  message("No site labels in ", siteLabelDirectory, "; running random splits only.")
  siteLabelDirectory <- NULL
}
tasks <- c("DvtTetraplegia", "BleedTetraPlegia")
started <- Sys.time()
directories <- Sys.glob(file.path("results/final", paste0(tasks, "_*")))
message(format(Sys.time(), "%H:%M"), "  Checking the saved models of ", length(directories), " result directories.")
runProductionChecks(directories, dataRoot, sites, "results/productionChecksTetraplegia")
message(format(Sys.time(), "%H:%M"), "  Pseudo-site analysis.")
pseudo <- runPseudoSiteAnalysis(dataRoot, sites, Sys.glob("results/final/*"), "results/pseudoSitesTetraplegia",
  tasks = tasks, cores = cores, siteLabelDirectory = siteLabelDirectory)
bundle <- Sys.glob(c("results/productionChecksTetraplegia/export/*.csv", "results/pseudoSitesTetraplegia/export/*.csv"))
utils::tar("early.tar.gz", files = bundle, compression = "gzip", tar = "internal")
message(sprintf("\n%s  Wrote early.tar.gz (%d files) after %.1f h. Largest check difference %.2g (should be near 0).",
  format(Sys.time(), "%H:%M"), length(bundle), as.numeric(difftime(Sys.time(), started, units = "hours")),
  max(pseudo$pseudoSites$check, na.rm = TRUE)))
