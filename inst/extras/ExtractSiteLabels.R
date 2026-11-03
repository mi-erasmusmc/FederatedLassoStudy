# Writes the natural-site label of every entry of the study's target cohorts,
# for the databases that have a usable site field (see ExploreSiteFields.R):
# provider network (care site of the last visit before index) in Optum EHR and
# state in CCAE and MDCR. MDCD and JMDC have none and stay whole. The files are
# patient-level and stay on this server: one per cohort and database,
# <outputDirectory>/<cohort id>/<database>.csv, with subjectId,
# cohortStartDate and label. runPseudoSiteAnalysis() reads them.
#
# From the directory the extraction ran in (holding extras/fetchTaskData.R and
# private/fetch_data_sources.yml), in R:
#
#   source("<path to this file>")
#
# Set `fetchScript`, `dataSourcesFile`, `cohortIds`, `fields` or
# `outputDirectory` before sourcing to change them. Reads the databases only.

if (!exists("fetchScript", inherits = FALSE)) fetchScript <- "extras/fetchTaskData.R"
if (!exists("dataSourcesFile", inherits = FALSE)) dataSourcesFile <- "private/fetch_data_sources.yml"
if (!exists("cohortIds", inherits = FALSE)) cohortIds <- c(tetraplegia = 15026L, warfarin = 13524L)
if (!exists("fields", inherits = FALSE)) fields <- c(OPTUMEHR = "lastVisitCareSite", CCAE = "state", MDCR = "state")
if (!exists("outputDirectory", inherits = FALSE)) outputDirectory <- "private/siteLabels"

source(fetchScript, local = (helpers <- new.env()))
config <- yaml::read_yaml(dataSourcesFile)

labelSql <- list(
  state = "SELECT CAST(c.person_id AS VARCHAR(50)) AS subject_id, c.cohort_start_date,
      CAST(l.state AS VARCHAR(255)) AS label FROM cohort c
    LEFT JOIN @cdm.person p ON p.person_id = c.person_id LEFT JOIN @cdm.location l ON l.location_id = p.location_id",
  lastVisitCareSite = "SELECT CAST(c.person_id AS VARCHAR(50)) AS subject_id, c.cohort_start_date,
      CAST(v.care_site_id AS VARCHAR(255)) AS label FROM cohort c
    LEFT JOIN (SELECT vo.person_id, c2.cohort_start_date, vo.care_site_id,
        ROW_NUMBER() OVER (PARTITION BY vo.person_id, c2.cohort_start_date ORDER BY vo.visit_start_date DESC) AS rn
      FROM @cdm.visit_occurrence vo JOIN cohort c2 ON vo.person_id = c2.person_id
      WHERE vo.visit_start_date <= c2.cohort_start_date AND vo.care_site_id IS NOT NULL) v
    ON v.person_id = c.person_id AND v.cohort_start_date = c.cohort_start_date AND v.rn = 1")

for (name in names(fields)) {
  dataSource <- config$dataSources[[name]]
  if (is.null(dataSource)) stop("No data source named ", name, " in ", dataSourcesFile, call. = FALSE)
  connection <- DatabaseConnector::connect(helpers$makeConnectionDetails(dataSource, config$connectionProfiles))
  for (cohort in names(cohortIds)) {
    labels <- DatabaseConnector::renderTranslateQuerySql(connection,
      "WITH cohort AS (SELECT subject_id AS person_id, cohort_start_date FROM @cohort_schema.@cohort_table
        WHERE cohort_definition_id = @cohort_id) @labelled;",
      labelled = labelSql[[fields[[name]]]], cdm = dataSource$cdmDatabaseSchema,
      cohort_schema = dataSource$cohortDatabaseSchema, cohort_table = dataSource$cohortTable,
      cohort_id = cohortIds[[cohort]], tempEmulationSchema = dataSource$tempEmulationSchema,
      snakeCaseToCamelCase = TRUE)
    directory <- file.path(outputDirectory, cohortIds[[cohort]])
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    labels$cohortStartDate <- format(as.Date(labels$cohortStartDate))
    utils::write.csv(labels, file.path(directory, paste0(name, ".csv")), row.names = FALSE)
    message(sprintf("%s  %s %s: %d entries, %d labels, %.1f%% unlabelled", format(Sys.time(), "%H:%M"), name,
      cohort, nrow(labels), length(unique(stats::na.omit(labels$label))), 100 * mean(is.na(labels$label))))
  }
  DatabaseConnector::disconnect(connection)
}
message("\nWrote the site labels to ", normalizePath(outputDirectory), ". They are patient-level: do not upload them.")
