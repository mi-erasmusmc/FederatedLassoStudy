# Which CDM fields could split each database into natural sites (regions,
# care sites, providers, insurers)? For the entries of the study's target
# cohorts, counts how many have a value of each candidate field, how many
# distinct values there are, and how entries spread over them. Only counts are
# written, never the values themselves. Reads the databases only; it can run
# while study sessions are running and needs no package installation.
#
# From the directory the extraction ran in (holding extras/fetchTaskData.R and
# private/fetch_data_sources.yml), in R:
#
#   source("<path to this file>")
#
# Set `fetchScript`, `dataSourcesFile`, `cohortIds` or `outputFile` before
# sourcing to change them. Writes siteFields.csv.

if (!exists("fetchScript", inherits = FALSE)) fetchScript <- "extras/fetchTaskData.R"
if (!exists("dataSourcesFile", inherits = FALSE)) dataSourcesFile <- "private/fetch_data_sources.yml"
if (!exists("cohortIds", inherits = FALSE)) cohortIds <- c(tetraplegia = 15026L, warfarin = 13524L)
if (!exists("outputFile", inherits = FALSE)) outputFile <- "siteFields.csv"

# The extraction's own connection helpers; the script's guard keeps it from running.
source(fetchScript, local = (helpers <- new.env()))
config <- yaml::read_yaml(dataSourcesFile)

# One row per cohort entry with the candidate label (NULL when absent).
labelSql <- list(
  personLocation = "SELECT c.person_id, CAST(p.location_id AS VARCHAR(255)) AS label FROM cohort c
    LEFT JOIN @cdm.person p ON p.person_id = c.person_id",
  state = "SELECT c.person_id, CAST(l.state AS VARCHAR(255)) AS label FROM cohort c
    LEFT JOIN @cdm.person p ON p.person_id = c.person_id LEFT JOIN @cdm.location l ON l.location_id = p.location_id",
  county = "SELECT c.person_id, CAST(l.county AS VARCHAR(255)) AS label FROM cohort c
    LEFT JOIN @cdm.person p ON p.person_id = c.person_id LEFT JOIN @cdm.location l ON l.location_id = p.location_id",
  personCareSite = "SELECT c.person_id, CAST(p.care_site_id AS VARCHAR(255)) AS label FROM cohort c
    LEFT JOIN @cdm.person p ON p.person_id = c.person_id",
  personProvider = "SELECT c.person_id, CAST(p.provider_id AS VARCHAR(255)) AS label FROM cohort c
    LEFT JOIN @cdm.person p ON p.person_id = c.person_id",
  lastVisitCareSite = "SELECT c.person_id, CAST(v.care_site_id AS VARCHAR(255)) AS label FROM cohort c
    LEFT JOIN (SELECT vo.person_id, c2.cohort_start_date, vo.care_site_id,
        ROW_NUMBER() OVER (PARTITION BY vo.person_id, c2.cohort_start_date ORDER BY vo.visit_start_date DESC) AS rn
      FROM @cdm.visit_occurrence vo JOIN cohort c2 ON vo.person_id = c2.person_id
      WHERE vo.visit_start_date <= c2.cohort_start_date AND vo.care_site_id IS NOT NULL) v
    ON v.person_id = c.person_id AND v.cohort_start_date = c.cohort_start_date AND v.rn = 1",
  payer = "SELECT c.person_id, MIN(CAST(pp.payer_source_value AS VARCHAR(255))) AS label FROM cohort c
    LEFT JOIN @cdm.payer_plan_period pp ON pp.person_id = c.person_id
      AND pp.payer_plan_period_start_date <= c.cohort_start_date AND pp.payer_plan_period_end_date >= c.cohort_start_date
    GROUP BY c.person_id, c.cohort_start_date",
  plan = "SELECT c.person_id, MIN(CAST(pp.plan_source_value AS VARCHAR(255))) AS label FROM cohort c
    LEFT JOIN @cdm.payer_plan_period pp ON pp.person_id = c.person_id
      AND pp.payer_plan_period_start_date <= c.cohort_start_date AND pp.payer_plan_period_end_date >= c.cohort_start_date
    GROUP BY c.person_id, c.cohort_start_date",
  sponsor = "SELECT c.person_id, MIN(CAST(pp.sponsor_source_value AS VARCHAR(255))) AS label FROM cohort c
    LEFT JOIN @cdm.payer_plan_period pp ON pp.person_id = c.person_id
      AND pp.payer_plan_period_start_date <= c.cohort_start_date AND pp.payer_plan_period_end_date >= c.cohort_start_date
    GROUP BY c.person_id, c.cohort_start_date")

summarise <- function(connection, dataSource, cohortId, field) {
  # The label query goes in as text, so that SqlRender sees every parameter in one template.
  cohort <- paste("WITH cohort AS (SELECT subject_id AS person_id, cohort_start_date FROM @cohort_schema.@cohort_table
    WHERE cohort_definition_id = @cohort_id), labelled AS (", labelSql[[field]], ")")
  sql <- paste(cohort, "SELECT label_n FROM (SELECT label, COUNT(*) AS label_n FROM labelled
    WHERE label IS NOT NULL GROUP BY label) t")
  totalSql <- paste(cohort, "SELECT COUNT(*) AS entries FROM labelled")
  run <- function(query) {
    DatabaseConnector::renderTranslateQuerySql(connection, query, cdm = dataSource$cdmDatabaseSchema, cohort_schema = dataSource$cohortDatabaseSchema,
      cohort_table = dataSource$cohortTable, cohort_id = cohortId, tempEmulationSchema = dataSource$tempEmulationSchema,
      snakeCaseToCamelCase = TRUE)
  }
  tryCatch({
    entries <- run(totalSql)$entries
    sizes <- sort(run(sql)$labelN, decreasing = TRUE)
    data.frame(entries = entries, labelled = sum(sizes), labels = length(sizes),
      labelsWith100 = sum(sizes >= 100), labelsWith1000 = sum(sizes >= 1000), labelsWith10000 = sum(sizes >= 10000),
      largestShare = if (length(sizes)) round(sizes[[1]] / sum(sizes), 3) else NA_real_,
      medianEntries = if (length(sizes)) stats::median(sizes) else NA_real_, error = NA_character_)
  }, error = function(e) data.frame(entries = NA, labelled = NA, labels = NA, labelsWith100 = NA,
    labelsWith1000 = NA, labelsWith10000 = NA, largestShare = NA, medianEntries = NA,
    error = substr(conditionMessage(e), 1, 200)))
}

rows <- list()
for (name in names(config$dataSources)) {
  dataSource <- config$dataSources[[name]]
  message(format(Sys.time(), "%H:%M"), "  ", name)
  connection <- DatabaseConnector::connect(helpers$makeConnectionDetails(dataSource, config$connectionProfiles))
  for (cohort in names(cohortIds)) for (field in names(labelSql)) {
    rows[[length(rows) + 1L]] <- cbind(database = name, cohort = cohort, field = field,
      summarise(connection, dataSource, cohortIds[[cohort]], field))
  }
  DatabaseConnector::disconnect(connection)
}
result <- do.call(rbind, rows)
utils::write.csv(result, outputFile, row.names = FALSE)
print(result[, c("database", "cohort", "field", "entries", "labelled", "labels", "labelsWith1000", "largestShare")],
  row.names = FALSE)
message("\nWrote ", outputFile, ". Upload it to OneDrive natural-sites/.")
