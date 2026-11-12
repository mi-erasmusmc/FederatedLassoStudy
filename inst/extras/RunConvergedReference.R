# Runs runConvergedReference() on the main analysis (all tasks, both predictor sets, phenotype
# predictors first) and bundles its aggregate tables for export. From the study's working directory
# (the one holding data/ and results/), in R:
#
#   source(system.file("extras", "RunConvergedReference.R", package = "FederatedLassoStudy"))
#
# It reads the saved main-analysis models and the site data, refits pooled lasso, the pooled data at
# DualAvg's penalty, every local lasso and the target-database lasso to convergence, scores every
# model against the converged pooled lasso with calibration curves, and writes
# results/convergedModels/export and convergedModels.tar.gz. The export tables are rewritten after
# each task and predictor set; if the run stops, sourcing the script again continues from there. Set
# `sites`, `dataRoot` or `threads` before sourcing to change them.

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
tables <- runConvergedReference(dataRoot, sites, Sys.glob("results/final/*"), "results/convergedModels",
  threads = threads)
bundle <- Sys.glob("results/convergedModels/export/*.csv")
utils::tar("convergedModels.tar.gz", files = bundle, compression = "gzip", tar = "internal")
fits <- tables$referenceFits
refits <- fits[startsWith(fits$fit, "converged"), ]
checks <- fits$check[fits$fit == "averageCheck"]
# The saved models must score as the main run reported them, and the target-database lasso as the
# production checks did.
reported <- do.call(rbind, lapply(Sys.glob("results/final/*"), function(directory) {
  r <- FederatedLearning::readComparisonResults(directory)
  r[is.finite(r$auc), c("task", "featureSet", "clientId", "method", "auc", "calibrationSlope", "calibrationIci")]
}))
scores <- tables$referenceScores
key <- function(d, column) paste(d$task, d$featureSet, d[[column]], d$model)
saved <- scores[scores$kind == "saved", ]
index <- match(key(saved, "heldOutDatabase"), paste(reported$task, reported$featureSet, reported$clientId, reported$method))
gap <- function(a, b) max(abs(a - b), na.rm = TRUE)
targetGap <- NA_real_
if (file.exists("results/productionChecks/export/auc.csv")) {
  previous <- utils::read.csv("results/productionChecks/export/auc.csv", stringsAsFactors = FALSE)
  previous <- previous[previous$method == "TargetLasso" & startsWith(previous$resultsDirectory, paste0(previous$task, "_")), ]
  target <- scores[scores$model == "TargetLasso", ]
  targetGap <- gap(target$auc, previous$auc[match(paste(target$task, target$featureSet, target$heldOutDatabase),
    paste(previous$task, previous$featureSet, previous$heldOutDatabase))])
}
message(sprintf(paste0("\n%s  Wrote convergedModels.tar.gz (%d files) after %.1f h. %d of %d refits failed; ",
  "target-database lasso failed in %d evaluations (1 expected). ",
  "Largest optimality violation of the refits %.1g (should be below 1e-6). ",
  "Largest coefficient-average check %.1g (should be below 1e-8). ",
  "Saved models reproduce the reported AUCs to %.1g, slopes to %.1g and ICI to %.1g, and the target-database lasso ",
  "its earlier AUCs to %.1g (all should be below 1e-6)."),
  format(Sys.time(), "%H:%M"), length(bundle), as.numeric(difftime(Sys.time(), started, units = "hours")),
  sum(!is.na(refits$error)), nrow(refits), sum(fits$fit == "targetLasso"), max(refits$kktMax, na.rm = TRUE),
  max(checks, na.rm = TRUE), gap(saved$auc, reported$auc[index]), gap(saved$calibrationSlope, reported$calibrationSlope[index]),
  gap(saved$calibrationIci, reported$calibrationIci[index]), targetGap))
