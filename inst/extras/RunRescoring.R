# Rescores saved models against the converged pooled lasso with the framework's root-finding
# calibration-in-the-large, and bundles the aggregate results for export. From the study's working
# directory (the one holding data/ and results/), in R:
#
#   source(system.file("extras", "RunRescoring.R", package = "FederatedLassoStudy"))
#
# It scores (1) every saved model of the main analysis (results/final), which gives a finite
# calibration-in-the-large for the ODAL1 fit whose glm calculation ran off, together with the
# converged models of results/convergedCalibrated rebuilt from their exported coefficients, so that
# every calibration curve is drawn on the framework's full-range grid; (2) the ODAL2 fits with
# condition predictors of the sensitivity analysis (results/sensitivity/odal2Conditions_*), so that
# Table S16 is measured against the converged reference; and (3) the centred DualAvg fits of
# RunCentredDualAvg.R (results/dualAvgCentred), if present. Only the held-out databases' data are
# read; nothing is refitted. Calibration curves leave out the 10 most extreme entries at each end. It needs results/convergedCalibrated from
# RunConvergedReference.R and writes results/rescored and rescored.tar.gz. Set `sites` or
# `dataRoot` before sourcing to change them.

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
referenceDirectory <- "results/convergedCalibrated"
if (!file.exists(file.path(referenceDirectory, "export", "referenceCoefficients.csv"))) {
  stop("Run RunConvergedReference.R first: ", referenceDirectory, "/export/referenceCoefficients.csv is missing",
    call. = FALSE)
}
started <- Sys.time()
main <- runReferenceScoring(dataRoot, sites, Sys.glob("results/final/*"), referenceDirectory, "results/rescored/main")
odal2 <- Sys.glob("results/sensitivity/odal2Conditions_*")
odal2 <- odal2[file.exists(file.path(odal2, "comparison_results.csv"))]
conditions <- if (length(odal2)) {
  runReferenceScoring(dataRoot, sites, odal2, referenceDirectory, "results/rescored/odal2Conditions",
    featureSets = "ageSexConditions", scoreConverged = FALSE)
}
centred <- Sys.glob("results/dualAvgCentred/*")
centred <- centred[file.exists(file.path(centred, "comparison_results.csv"))]
centredScores <- if (length(centred)) {
  runReferenceScoring(dataRoot, sites, centred, referenceDirectory, "results/rescored/dualAvgCentred",
    scoreConverged = FALSE)
}
bundle <- Sys.glob("results/rescored/*/export/*.csv")
utils::tar("rescored.tar.gz", files = bundle, compression = "gzip", tar = "internal")

scores <- main$referenceScores
odal1 <- scores[scores$model == "ODAL1" & scores$task == "BleedWarfarin" & scores$featureSet == "ageSexPhenotypes" &
  scores$heldOutDatabase == "CCAE", ]
message(sprintf(paste0("\nWrote rescored.tar.gz (%d files) in %.1f h. Main analysis: %d models scored in %d evaluations; ",
  "calibration-in-the-large finite for %d of %d. ODAL1, warfarin, CCAE held out: %.2f (finite; NA before the fix). ",
  "ODAL2 with condition predictors: %d models scored. Converged models rescored: %d. ",
  "Centred DualAvg: %d models scored."), length(bundle),
  as.numeric(difftime(Sys.time(), started, units = "hours")), sum(scores$kind == "saved"),
  length(unique(paste(scores$task, scores$featureSet, scores$heldOutDatabase))),
  sum(is.finite(scores$calibrationInTheLarge)), nrow(scores), odal1$calibrationInTheLarge[1],
  if (is.null(conditions)) 0L else sum(conditions$referenceScores$kind == "saved"),
  sum(scores$kind == "converged" & scores$model != "convergedPooled"),
  if (is.null(centredScores)) 0L else sum(centredScores$referenceScores$kind == "saved")))
