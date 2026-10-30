#' Study report figures
#'
#' Draws the main manuscript figures from the tables of [writeStudyReport()]:
#' - `discrimination`: each workflow's AUC difference from pooled lasso in
#'   every task and held-out database, with the mean, by predictor set.
#'   Flagged fits are drawn hollow;
#' - `localBaselines`: per task and held-out database, the AUCs of the
#'   single-source models beside coefficient averaging, probability
#'   ensembling and pooled lasso, by predictor set.
#'
#' Each figure is written as PDF and PNG. Requires ggplot2.
#'
#' @param reportDirectory Output of [writeStudyReport()].
#' @param outputDirectory Directory for the figures; defaults to `reportDirectory`.
#' @return The figure paths, invisibly.
#' @export
writeStudyFigures <- function(reportDirectory, outputDirectory = reportDirectory) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("writeStudyFigures() needs the ggplot2 package", call. = FALSE)
  }
  evaluations <- utils::read.csv(file.path(reportDirectory, "evaluations.csv"), stringsAsFactors = FALSE)
  local <- utils::read.csv(file.path(reportDirectory, "localBaselines.csv"), stringsAsFactors = FALSE)
  dir.create(outputDirectory, recursive = TRUE, showWarnings = FALSE)
  figures <- list(
    discrimination = list(plot = discriminationFigure(evaluations), width = 7.5, height = 5),
    localBaselines = list(plot = localBaselineFigure(local, evaluations), width = 7.5, height = 6.5)
  )
  paths <- character()
  for (name in names(figures)) {
    for (extension in c("pdf", "png")) {
      path <- file.path(outputDirectory, paste0(name, ".", extension))
      ggplot2::ggsave(path, figures[[name]]$plot, width = figures[[name]]$width,
        height = figures[[name]]$height, dpi = 300, bg = "white")
      paths <- c(paths, path)
    }
  }
  invisible(paths)
}

# Validated categorical slots 1-3 (all pairs, light surface); shapes carry
# the same identity for print and colour-vision deficiency.
figureTaskColors <- c("#2a78d6", "#eb6834", "#1baf7a")
figureTaskShapes <- c(16, 17, 15)
figureTaskShapesHollow <- c(1, 2, 0)
figureInk <- c(primary = "#0b0b0b", secondary = "#52514e", muted = "#8a8984", grid = "#e6e5e1")

figureLabels <- list(
  task = c(DvtTetraplegia = "DVT, tetraplegia", BleedTetraPlegia = "Bleeding, tetraplegia",
    BleedWarfarin = "Bleeding, warfarin"),
  featureSet = c(ageSexPhenotypes = "Phenotype predictors", ageSexConditions = "Condition predictors"),
  method = c(DualAvg = "DualAvg", ODAL2 = "ODAL2", ADAP2 = "ADAP2", ODAL1 = "ODAL1", ADAP1 = "ADAP1",
    LocalAvgLasso = "Coefficient average", LocalEnsembleLasso = "Probability ensemble",
    BiggestSiteLasso = "Largest database")
)

figureTheme <- function() {
  ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = figureInk[["primary"]]),
      axis.text = ggplot2::element_text(colour = figureInk[["secondary"]]),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(colour = figureInk[["grid"]], linewidth = 0.3),
      strip.text = ggplot2::element_text(face = "bold", hjust = 0),
      legend.position = "top", legend.justification = "left",
      legend.title = ggplot2::element_blank(),
      panel.spacing = ggplot2::unit(12, "pt"),
      plot.margin = ggplot2::margin(5, 16, 5, 5)
    )
}

discriminationFigure <- function(evaluations) {
  .data <- ggplot2::.data
  methods <- names(figureLabels$method)
  federated <- c("DualAvg", "ODAL2", "ADAP2", "ODAL1", "ADAP1")
  groups <- c("Federated", "Local-model baselines")
  rows <- evaluations[evaluations$method %in% methods & is.finite(evaluations$aucDifference), , drop = FALSE]
  rows$group <- factor(ifelse(rows$method %in% federated, groups[[1]], groups[[2]]), levels = groups)
  rows$method <- factor(figureLabels$method[rows$method], levels = rev(figureLabels$method))
  # Each task keeps its colour and shape when others are absent.
  rows$task <- droplevels(factor(figureLabels$task[rows$task], levels = figureLabels$task))
  taskSlot <- match(levels(rows$task), figureLabels$task)
  rows$featureSet <- factor(figureLabels$featureSet[rows$featureSet], levels = figureLabels$featureSet)
  rows$shape <- interaction(rows$task, rows$flagged)
  means <- stats::aggregate(aucDifference ~ method + featureSet + group, rows, mean)
  # Methods fitted for one predictor set only are labelled in the other panel.
  fitted <- unique(rows[c("method", "featureSet", "group")])
  missing <- merge(unique(rows[c("method", "group")]), data.frame(featureSet = levels(rows$featureSet)))
  missing <- missing[!paste(missing$method, missing$featureSet) %in% paste(fitted$method, fitted$featureSet), ]
  missing$featureSet <- factor(missing$featureSet, levels = levels(rows$featureSet))
  shapes <- stats::setNames(c(figureTaskShapes[taskSlot], figureTaskShapesHollow[taskSlot]),
    c(paste(levels(rows$task), "FALSE", sep = "."), paste(levels(rows$task), "TRUE", sep = ".")))
  signed <- function(x) ifelse(x > 0, sprintf("+%.3f", x), sprintf("%.3f", x))
  ggplot2::ggplot(rows, ggplot2::aes(x = .data$aucDifference, y = .data$method)) +
    ggplot2::geom_vline(xintercept = 0, colour = figureInk[["muted"]], linewidth = 0.4) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$task, shape = .data$shape), size = 1.8,
      position = ggplot2::position_jitter(height = 0.18, width = 0, seed = 1), stroke = 0.6) +
    ggplot2::geom_point(data = means, shape = 124, size = 5, colour = figureInk[["primary"]]) +
    ggplot2::geom_text(data = missing, ggplot2::aes(x = 0, label = "Not fitted"), size = 2.6,
      colour = figureInk[["secondary"]], fontface = "italic") +
    ggplot2::facet_grid(group ~ featureSet, scales = "free_y", space = "free_y", switch = "y") +
    ggplot2::scale_colour_manual(values = stats::setNames(figureTaskColors[taskSlot], levels(rows$task))) +
    ggplot2::scale_shape_manual(values = shapes, guide = "none") +
    ggplot2::scale_x_continuous(labels = signed) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(shape = figureTaskShapes[taskSlot], size = 2.2))) +
    ggplot2::labs(x = "AUC difference from pooled lasso in the held-out database", y = NULL,
      caption = "Each marker is one task and held-out database; the black tick is the mean. Hollow: flagged fit.") +
    figureTheme() +
    ggplot2::theme(strip.placement = "outside", strip.text.y.left = ggplot2::element_text(angle = 90),
      plot.caption = ggplot2::element_text(colour = figureInk[["secondary"]], hjust = 0))
}

localBaselineFigure <- function(local, evaluations) {
  .data <- ggplot2::.data
  local$task <- factor(figureLabels$task[local$task], levels = figureLabels$task)
  local$featureSet <- factor(figureLabels$featureSet[local$featureSet], levels = figureLabels$featureSet)
  local$row <- factor(local$heldOutDatabase, levels = rev(sort(unique(local$heldOutDatabase))))
  sources <- evaluations[evaluations$method == "LocalSiteLasso" & evaluations$available, , drop = FALSE]
  sources$task <- factor(figureLabels$task[sources$task], levels = figureLabels$task)
  sources$featureSet <- factor(figureLabels$featureSet[sources$featureSet], levels = figureLabels$featureSet)
  sources$row <- factor(sources$heldOutDatabase, levels = levels(local$row))
  aggregates <- rbind(
    data.frame(local[c("task", "featureSet", "row")], model = "Coefficient average", auc = local$coefficientAverageAuc),
    data.frame(local[c("task", "featureSet", "row")], model = "Probability ensemble", auc = local$probabilityEnsembleAuc),
    data.frame(local[c("task", "featureSet", "row")], model = "Pooled lasso", auc = local$pooledAuc))
  aggregates$model <- factor(aggregates$model, levels = c("Coefficient average", "Probability ensemble", "Pooled lasso"))
  # Coefficient average above and probability ensemble below the row, so
  # near-equal AUCs stay visible.
  aggregates$nudge <- c(`Coefficient average` = 0.2, `Probability ensemble` = -0.2, `Pooled lasso` = 0)[
    as.character(aggregates$model)]
  ggplot2::ggplot(local, ggplot2::aes(y = .data$row)) +
    ggplot2::geom_segment(ggplot2::aes(x = .data$minSourceAuc, xend = .data$maxSourceAuc, yend = .data$row),
      colour = figureInk[["grid"]], linewidth = 2.2, lineend = "round") +
    ggplot2::geom_point(data = sources, ggplot2::aes(x = .data$auc), colour = figureInk[["muted"]], size = 1.1) +
    ggplot2::geom_point(data = aggregates, ggplot2::aes(x = .data$auc, y = as.numeric(.data$row) + .data$nudge,
      colour = .data$model, shape = .data$model), size = 2, stroke = 0.7) +
    ggplot2::facet_grid(task ~ featureSet, scales = "free_y", space = "free_y", switch = "y") +
    ggplot2::scale_colour_manual(values = c(`Coefficient average` = figureTaskColors[[1]],
      `Probability ensemble` = figureTaskColors[[2]], `Pooled lasso` = figureInk[["primary"]])) +
    ggplot2::scale_shape_manual(values = c(`Coefficient average` = 16, `Probability ensemble` = 17,
      `Pooled lasso` = 3)) +
    ggplot2::labs(x = "AUC in the held-out database", y = NULL,
      caption = "Grey: models fitted in one training database (range and each model).") +
    figureTheme() +
    ggplot2::theme(strip.placement = "outside", strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 0),
      plot.caption = ggplot2::element_text(colour = figureInk[["secondary"]], hjust = 0))
}
