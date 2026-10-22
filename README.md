# FederatedLassoStudy

Study package for comparing federated, local and pooled logistic prediction
workflows across five OMOP CDM databases, with leave-one-database-out external
validation. Models are fitted with the
[FederatedLearning](https://github.com/mi-erasmusmc/FederatedLearning)
framework, pinned to commit `5b38a6b` (after release `v0.3.0`).

This package defines the study: its prediction tasks, predictors and method
settings. Reusable methods live in the framework. Database connection
details, credentials and patient-level data are never part of this
repository.

## Prediction Tasks

| Task | Target cohort | Outcome cohort | Risk window (days) |
|---|---|---|---|
| `DvtTetraplegia` | Tetraplegia (15026) | Deep vein thrombosis (19047) | 1-365 |
| `BleedTetraPlegia` | Tetraplegia (15026) | Bleeding (19012) | 1-365 |
| `BleedWarfarin` | Warfarin with atrial fibrillation (13524) | Bleeding (19012) | 1-30 |

Patients with a prior outcome are excluded, and one day at risk suffices.
`studyTasks()` returns these definitions. A test checks every horizon against
the protocol.

## Predictors

`studyFeatureSettings()` defines the following predictors:
- age;
- the female indicator;
- binary PhenotypeLibrary phenotype covariates (analysis 49, days -365
  through 0);
- one indicator for any antibiotic (covariate `999001049`), replacing the 14
  antibiotic-family phenotypes 1201-1214 by their logical OR. A family that no
  patient in a database has counts as zero there.

A shared training-only filter removes covariates present in fewer than 0.1% of
outer-training entries.

## Preprocessing Audit

`runPreprocessingAudit()` builds the shared training-only feature filter for
every task and held-out database without fitting models:

```r
FederatedLassoStudy::runPreprocessingAudit(
  dataRoot = "data",
  siteIds = c("siteA", "siteB", "siteC", "siteD", "siteE"),
  outputDirectory = "results/preprocessingAudit"
)
```

The data must be laid out as `<dataRoot>/<task><extraction>/<site>`, with no
suffix for the phenotype extraction and `AgeSexConditions` or `Default` for the
wider ones; `taskDirectories`
maps tasks to differently named directories. The `preprocessing` directory
holds per-predictor training counts and stays private. The `export` directory
holds only predictor counts and covariate ids.

## Predictor Specifications

`studyFeatureSettings()` names three specifications, each matching its own
extraction: the protocol's `ageSexPhenotypes`, and two wider candidates,
`ageSexConditions` (age, female and conditions in the prior year) and
`defaultCovariates` (FeatureExtraction's defaults with a single sex
indicator). Pass `featureSets` to `runPreprocessingAudit()` or `runStudy()`.

## Running

`runStudy()` runs the comparison with `studyMethods()` and
`studyMethodSettings()`. Pilots before the final run:

```r
# Matched pooled lasso and DualAvg penalty selection on one outer fold.
runTuningPilot("data", sites, "results/tuningPilot",
  task = "DvtTetraplegia", holdout = "JMDC")
# DualAvg at that penalty with stopping thresholds 1e-10 and 1e-8, on the
# outer-training databases only.
runDualAvgThresholdPilot("data", sites, "results/tuningPilot", "results/thresholdPilot",
  task = "DvtTetraplegia", holdout = "JMDC")
```

Each writes an `export` directory of aggregate summaries for review.

## Status

Pending from the final-run plan: the stopping threshold, the predictor
specification, a pinned `renv.lock` and reporting scripts.

## Installation

```r
remotes::install_github("mi-erasmusmc/FederatedLassoStudy")
```

## License

Apache License 2.0.
