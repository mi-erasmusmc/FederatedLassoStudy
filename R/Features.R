#' Study predictor specification
#'
#' Age (covariate 1002, divided by 100 when loaded), the female indicator
#' (8532001) and binary PhenotypeLibrary phenotype covariates (analysis 49,
#' days -365 through 0). A single sex indicator avoids collinearity with the
#' intercept. The shared training-only filter removes covariates present in
#' fewer than 0.1% of outer-training entries.
#'
#' @return A named list of `featureSettings` objects.
#' @export
studyFeatureSettings <- function() {
  list(
    ageSexPhenotypes = FederatedLearning::createFeatureSettings(
      covariateIds = c(1002, 8532001),
      analysisIds = 49,
      minFraction = 0.001
    )
  )
}
