# FederatedLassoStudy

Study package for comparing federated, local and pooled logistic prediction
workflows across five OMOP CDM databases, with leave-one-database-out external
validation. Models are fitted with the
[FederatedLearning](https://github.com/mi-erasmusmc/FederatedLearning)
framework, pinned to release `v0.2.0`.

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
  through 0).

A shared training-only filter removes covariates present in fewer than 0.1% of
outer-training entries.

## Status

Pending changes from the final-run plan:
- Combine the 14 antibiotic phenotypes into one indicator.
- Final method settings: matched pooled/DualAvg database-level tuning, log
  loss for penalty selection, three-fold local lasso CV.
- External calibration outputs.
- Run and reporting scripts, with a pinned `renv.lock`.

## Installation

```r
remotes::install_github("mi-erasmusmc/FederatedLassoStudy")
```

## License

Apache License 2.0.
