#' Study predictor specifications
#'
#' Each specification needs its own production extraction:
#' - `ageSexPhenotypes`, the protocol's predictors: age (covariate 1002,
#'   divided by 100 when loaded), the female indicator (8532001) and binary
#'   PhenotypeLibrary phenotype covariates (analysis 49, days -365 through 0).
#'   The 14 antibiotic-family phenotypes (PhenotypeLibrary 1201-1214) are
#'   replaced by one indicator for any antibiotic, their logical OR, with
#'   covariate id 999001049.
#' - `ageSexPhenotypesNoAtrialFibrillation`, a sensitivity analysis for
#'   bleeding after warfarin, whose cohort requires atrial fibrillation: the
#'   protocol's predictors without the atrial fibrillation phenotype
#'   (PhenotypeLibrary 1160). It uses the phenotype extraction.
#' - `ageSexConditions`, a wider candidate: age, the female indicator and
#'   condition occurrence in the prior year (analysis 102).
#' - `defaultCovariates`, the widest candidate: FeatureExtraction's default
#'   covariate analyses with the female indicator in place of both sex
#'   indicators. It has age groups rather than age, and includes index year,
#'   race and continuous risk scores, which are not rescaled.
#'
#' A single sex indicator avoids collinearity with the intercept. The shared
#' training-only filter removes covariates present in fewer than 0.1% of
#' outer-training entries.
#'
#' @param featureSets Names of the specifications to return.
#' @return A named list of `featureSettings` objects.
#' @export
studyFeatureSettings <- function(featureSets = "ageSexPhenotypes") {
  phenotypes <- function(excludeCovariateIds = NULL) {
    FederatedLearning::createFeatureSettings(
      covariateIds = c(ageCovariateId, femaleCovariateId),
      analysisIds = 49,
      minFraction = 0.001,
      combinations = list(list(
        covariateId = anyAntibioticCovariateId,
        name = "Any antibiotic",
        components = antibioticCovariateIds
      )),
      excludeCovariateIds = excludeCovariateIds
    )
  }
  all <- list(
    ageSexPhenotypes = phenotypes(),
    ageSexPhenotypesNoAtrialFibrillation = phenotypes(atrialFibrillationCovariateId),
    ageSexConditions = FederatedLearning::createFeatureSettings(
      covariateIds = c(ageCovariateId, femaleCovariateId),
      analysisIds = 102,
      minFraction = 0.001
    ),
    defaultCovariates = FederatedLearning::createFeatureSettings(
      covariateIds = femaleCovariateId,
      analysisIds = defaultAnalysisIds,
      minFraction = 0.001
    )
  )
  if (!length(featureSets) || !all(featureSets %in% names(all))) {
    stop("featureSets must name: ", paste(names(all), collapse = ", "), call. = FALSE)
  }
  all[featureSets]
}

ageCovariateId <- 1002
femaleCovariateId <- 8532001

# PhenotypeLibrary antibiotic-family phenotypes 1201-1214 as cohort covariates.
antibioticCovariateIds <- (1201:1214) * 1000 + 49

# PhenotypeLibrary 1160, atrial fibrillation, as a cohort covariate.
atrialFibrillationCovariateId <- 1160049

# Covariate id of the combined antibiotic indicator, outside PhenotypeLibrary's range.
anyAntibioticCovariateId <- 999001049

# FeatureExtraction::createDefaultCovariateSettings() analyses other than gender (1).
defaultAnalysisIds <- c(3:7, 210, 212, 410, 412, 413, 502, 504, 602, 604, 702, 704, 710,
  712, 714, 716, 802, 804, 806, 808, 901:904)
