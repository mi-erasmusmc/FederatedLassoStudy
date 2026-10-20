#' Study prediction tasks
#'
#' The three tasks of the final comparison. Cohort ids are the study's ATLAS
#' definitions: tetraplegia target 15026, warfarin target 13524 (first warfarin
#' exposure with atrial fibrillation in days -30 to 0), bleeding outcome 19012
#' and DVT outcome 19047. Every population setting is stated explicitly,
#' because PatientLevelPrediction's defaults differ from the protocol.
#'
#' @return A named list of `taskSettings` objects.
#' @export
studyTasks <- function() {
  task <- function(name, targetId, outcomeId, riskWindowEnd) {
    FederatedLearning::createTaskSettings(
      name = name,
      targetId = targetId,
      outcomeId = outcomeId,
      populationSettings = studyPopulationSettings(riskWindowEnd)
    )
  }
  list(
    DvtTetraplegia = task("DvtTetraplegia", 15026, 19047, riskWindowEnd = 365L),
    BleedTetraPlegia = task("BleedTetraPlegia", 15026, 19012, riskWindowEnd = 365L),
    BleedWarfarin = task("BleedWarfarin", 13524, 19012, riskWindowEnd = 30L)
  )
}

#' Study population settings
#'
#' Outcome risk from day 1 to `riskWindowEnd` after cohort start, excluding
#' patients with any prior outcome. Full follow-up is not required: one day at
#' risk suffices.
#'
#' @param riskWindowEnd Last day of the risk window, relative to cohort start.
#' @return A PatientLevelPrediction `populationSettings` object.
#' @export
studyPopulationSettings <- function(riskWindowEnd) {
  PatientLevelPrediction::createStudyPopulationSettings(
    binary = TRUE,
    includeAllOutcomes = TRUE,
    firstExposureOnly = FALSE,
    washoutPeriod = 0,
    removeSubjectsWithPriorOutcome = TRUE,
    priorOutcomeLookback = 99999L,
    requireTimeAtRisk = FALSE,
    minTimeAtRisk = 1L,
    riskWindowStart = 1L,
    startAnchor = "cohort start",
    riskWindowEnd = riskWindowEnd,
    endAnchor = "cohort start",
    restrictTarToCohortEnd = FALSE
  )
}
