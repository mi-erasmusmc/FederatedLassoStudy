#' Study report figures
#'
#' Draws the main manuscript figures from the tables of [writeStudyReport()]:
#' - `discrimination`: each workflow's AUC difference from pooled lasso in
#'   every task and held-out database, with the mean, by predictor set.
#'   Flagged fits are drawn hollow;
#' - `localBaselines`: per task and held-out database, the AUCs of the
#'   single-source models beside coefficient averaging, probability
#'   ensembling and pooled lasso, by predictor set;
#' - `calibration`, with `calibrationCurveDirectory`: (A) each workflow's
#'   calibration-in-the-large and slope in every task and held-out database,
#'   and (B) the curves of pooled lasso, DualAvg and ODAL1 for bleeding after
#'   warfarin, both with phenotype predictors.
#'
#' Each figure is written as PDF and PNG. Requires ggplot2.
#'
#' @param reportDirectory Output of [writeStudyReport()].
#' @param outputDirectory Directory for the figures; defaults to `reportDirectory`.
#' @param calibrationCurveDirectory Output of [runCalibrationCurveExport()], or
#'   `NULL` to skip the calibration figure.
#' @return The figure paths, invisibly.
#' @export
writeStudyFigures <- function(reportDirectory, outputDirectory = reportDirectory,
                              calibrationCurveDirectory = NULL) {
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
  if (!is.null(calibrationCurveDirectory)) {
    curves <- utils::read.csv(file.path(calibrationCurveDirectory, "export", "calibrationCurves.csv"),
      stringsAsFactors = FALSE)
    figures$calibration <- list(plot = calibrationSummaryFigure(evaluations, curves), width = 7.5, height = 9)
  }
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

#' Calibration curve figures
#'
#' Draws the curves of [runCalibrationCurveExport()], one figure per task and
#' predictor set: a row per held-out database and a column per workflow, each
#' curve with its pointwise 95% interval over pooled lasso and the diagonal of
#' ideal calibration. Both axes are logarithmic and shared within a figure;
#' a panel notes the curve points beyond the axes. Single-database models are
#' left out. Each figure is written as PDF and PNG. Requires ggplot2.
#'
#' @param calibrationCurveDirectory Output of [runCalibrationCurveExport()].
#' @param outputDirectory Directory for the figures.
#' @return The figure paths, invisibly.
#' @export
writeCalibrationFigures <- function(calibrationCurveDirectory, outputDirectory) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("writeCalibrationFigures() needs the ggplot2 package", call. = FALSE)
  }
  curves <- utils::read.csv(file.path(calibrationCurveDirectory, "export", "calibrationCurves.csv"),
    stringsAsFactors = FALSE)
  curves <- curves[curves$method %in% c("PooledLasso", names(figureLabels$method)), , drop = FALSE]
  dir.create(outputDirectory, recursive = TRUE, showWarnings = FALSE)
  paths <- character()
  groups <- unique(curves[c("task", "featureSet")])
  for (i in seq_len(nrow(groups))) {
    rows <- curves[curves$task == groups$task[[i]] & curves$featureSet == groups$featureSet[[i]], , drop = FALSE]
    plot <- calibrationFigure(rows)
    for (extension in c("pdf", "png")) {
      path <- file.path(outputDirectory,
        paste0("calibration_", groups$task[[i]], "_", groups$featureSet[[i]], ".", extension))
      ggplot2::ggsave(path, plot, width = 10, height = 6.5, dpi = 300, bg = "white")
      paths <- c(paths, path)
    }
  }
  invisible(paths)
}

# Main calibration figure, phenotype predictors: the per-evaluation measures
# show that held-out databases differ more than workflows, and the warfarin
# curves show DualAvg against ODAL1 beside pooled lasso.
calibrationSummaryFigure <- function(evaluations, curves, featureSet = "ageSexPhenotypes",
                                     curveTask = "BleedWarfarin") {
  .data <- ggplot2::.data
  highlighted <- c(PooledLasso = "Pooled lasso", DualAvg = "DualAvg", ODAL1 = "ODAL1")
  groups <- c(highlighted, other = "Other workflows")
  colours <- stats::setNames(c(figureInk[["primary"]], figureTaskColors[1:2], figureInk[["muted"]]), groups)
  groupOf <- function(method) {
    factor(ifelse(method %in% names(highlighted), highlighted[method], groups[["other"]]), levels = groups)
  }
  methods <- c("PooledLasso", names(figureLabels$method))
  rows <- evaluations[evaluations$featureSet == featureSet & evaluations$method %in% methods &
    is.finite(evaluations$calibrationInTheLarge), , drop = FALSE]
  rows$group <- groupOf(rows$method)
  rows$task <- factor(figureLabels$task[rows$task], levels = figureLabels$task)
  rows$database <- factor(rows$heldOutDatabase, levels = rev(sort(unique(rows$heldOutDatabase))))
  measures <- c("Calibration-in-the-large", "Calibration slope")
  long <- rbind(
    data.frame(rows[c("task", "database", "group", "flagged")], measure = measures[[1]],
      value = rows$calibrationInTheLarge, ideal = 0),
    data.frame(rows[c("task", "database", "group", "flagged")], measure = measures[[2]],
      value = rows$calibrationSlope, ideal = 1))
  long$measure <- factor(long$measure, levels = measures)
  long$shape <- paste(long$group, long$flagged)
  # Highlighted workflows sit above and below the row, other workflows on it.
  long$y <- as.numeric(long$database) + c(0, 0.22, -0.22, 0)[as.integer(long$group)]
  long <- long[order(long$group != groups[["other"]]), ]
  shapeValues <- stats::setNames(c(3, 3, 16, 1, 17, 2, 1, 1), paste(rep(groups, each = 2), c(FALSE, TRUE)))
  ideal <- unique(long[c("measure", "ideal")])
  levels <- levels(long$database)
  measuresPlot <- ggplot2::ggplot(long, ggplot2::aes(x = .data$value, y = .data$y)) +
    ggplot2::geom_vline(data = ideal, ggplot2::aes(xintercept = .data$ideal), colour = figureInk[["muted"]],
      linewidth = 0.4) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$group, shape = .data$shape), size = 1.7, stroke = 0.6) +
    ggplot2::facet_grid(task ~ measure, scales = "free_x", switch = "y") +
    ggplot2::scale_y_continuous(breaks = seq_along(levels), labels = levels) +
    ggplot2::scale_colour_manual(values = colours, drop = FALSE) +
    ggplot2::scale_shape_manual(values = shapeValues, guide = "none") +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(shape = c(3, 16, 17, 1), size = 2.2))) +
    ggplot2::labs(x = NULL, y = NULL, tag = "A",
      caption = paste0("Held-out database per row. Ideal: 0 for calibration-in-the-large, 1 for the slope.\n",
        "Negative calibration-in-the-large means predicted risk above observed. Hollow triangle: flagged ODAL1 fit.")) +
    figureTheme() +
    ggplot2::theme(strip.placement = "outside", strip.text.y.left = ggplot2::element_text(angle = 90),
      plot.caption = ggplot2::element_text(colour = figureInk[["secondary"]], hjust = 0))

  curve <- curves[curves$featureSet == featureSet & curves$task == curveTask &
    curves$method %in% names(highlighted), , drop = FALSE]
  curve$group <- groupOf(curve$method)
  curve$database <- factor(curve$heldOutDatabase, levels = sort(unique(curve$heldOutDatabase)))
  # The axes cover every curve point below 50% risk; points above are noted.
  shown <- curve$predicted <= 0.5 & curve$observed <= 0.5
  limits <- c(min(curve$predicted, curve$observed) / 1.15,
    max(curve$predicted[shown], curve$observed[shown]) * 1.15)
  outside <- !shown
  offScale <- stats::aggregate(list(points = outside), curve[c("database", "group")], sum)
  offScale <- offScale[offScale$points > 0, , drop = FALSE]
  offScale$label <- sprintf("%s: %d curve %s\nat %s predicted risk", offScale$group, offScale$points,
    ifelse(offScale$points == 1, "point", "points"),
    vapply(seq_len(nrow(offScale)), function(i) {
      at <- curve$predicted[outside & curve$database == offScale$database[[i]] & curve$group == offScale$group[[i]]]
      at <- unique(formatC(100 * range(at), format = "f", digits = 1))
      paste0(paste(at, collapse = "-"), "%")
    }, character(1)))
  breaks <- c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5)
  percent <- function(x) paste0(formatC(100 * x, format = "fg", digits = 2), "%")
  curvePlot <- ggplot2::ggplot(curve, ggplot2::aes(x = .data$predicted, y = .data$observed, colour = .data$group)) +
    ggplot2::geom_abline(colour = figureInk[["muted"]], linewidth = 0.3, linetype = "22") +
    ggplot2::geom_line(linewidth = 0.55) +
    ggplot2::geom_text(data = offScale, ggplot2::aes(label = .data$label), x = -Inf, y = Inf, hjust = -0.04,
      vjust = 1.4, size = 2.1, colour = figureInk[["secondary"]], inherit.aes = FALSE) +
    ggplot2::facet_grid(. ~ database) +
    ggplot2::scale_colour_manual(values = colours, guide = "none") +
    ggplot2::scale_x_log10(breaks = breaks, labels = percent) +
    ggplot2::scale_y_log10(breaks = breaks, labels = percent) +
    ggplot2::coord_cartesian(xlim = limits, ylim = limits, expand = FALSE) +
    ggplot2::labs(x = "Predicted risk", y = "Observed risk (smoothed)", tag = "B",
      title = paste0(figureLabels$task[[curveTask]], ": calibration curves by held-out database"),
      caption = paste("Curves span the 1st to 99th percentile of predicted risk in the held-out database.",
        "Dashed: ideal calibration. Logarithmic axes.")) +
    figureTheme() +
    ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA, colour = figureInk[["grid"]]),
      panel.grid.major.y = ggplot2::element_line(colour = figureInk[["grid"]], linewidth = 0.3),
      panel.spacing = ggplot2::unit(6, "pt"), plot.title = ggplot2::element_text(size = 9),
      strip.text = ggplot2::element_text(face = "bold", hjust = 0.5),
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
      plot.caption = ggplot2::element_text(colour = figureInk[["secondary"]], hjust = 0))
  # Building the grobs needs a device; a null one leaves no Rplots.pdf behind.
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grobTree(
    grid::grobTree(ggplot2::ggplotGrob(measuresPlot), vp = grid::viewport(y = 1, height = 0.68, just = "top")),
    grid::grobTree(ggplot2::ggplotGrob(curvePlot), vp = grid::viewport(y = 0, height = 0.32, just = "bottom")))
}

calibrationFigure <- function(rows) {
  .data <- ggplot2::.data
  labels <- c(PooledLasso = "Pooled lasso", figureLabels$method)
  labels <- labels[names(labels) %in% rows$method]
  rows$panel <- factor(labels[rows$method], levels = labels)
  rows$database <- factor(rows$heldOutDatabase, levels = sort(unique(rows$heldOutDatabase)))
  rows$series <- "Workflow"
  # Pooled lasso is drawn in every panel as the reference.
  pooled <- rows[rows$method == "PooledLasso", setdiff(names(rows), "panel"), drop = FALSE]
  pooled <- merge(pooled, data.frame(panel = factor(labels, levels = labels)))
  pooled$series <- "Pooled lasso"
  rows$series[rows$method == "PooledLasso"] <- "Pooled lasso"
  rows$fill <- rows$series
  # Axes cover nearly all curve points; the rest are counted in their panel.
  range <- stats::quantile(c(rows$predicted, rows$observed), c(0.005, 0.995), names = FALSE)
  limits <- c(range[[1]] / 1.3, min(1, range[[2]] * 1.3))
  outside <- rows$predicted < limits[[1]] | rows$predicted > limits[[2]] |
    rows$observed < limits[[1]] | rows$observed > limits[[2]]
  offScale <- stats::aggregate(list(points = outside), rows[c("panel", "database")], sum)
  offScale <- offScale[offScale$points > 0, , drop = FALSE]
  offScale$label <- sprintf("%d %s off scale", offScale$points, ifelse(offScale$points == 1, "point", "points"))
  breaks <- c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5)
  percent <- function(x) paste0(formatC(100 * x, format = "fg", digits = 2), "%")
  ggplot2::ggplot(rows, ggplot2::aes(x = .data$predicted, y = .data$observed)) +
    ggplot2::geom_abline(colour = figureInk[["muted"]], linewidth = 0.3, linetype = "22") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper, fill = .data$fill), alpha = 0.18) +
    ggplot2::geom_line(data = pooled, ggplot2::aes(colour = .data$series), linewidth = 0.35) +
    ggplot2::geom_line(ggplot2::aes(colour = .data$series), linewidth = 0.55) +
    ggplot2::geom_text(data = offScale, ggplot2::aes(label = .data$label), x = -Inf, y = Inf,
      hjust = -0.05, vjust = 1.4, size = 2.1, colour = figureInk[["secondary"]], inherit.aes = FALSE) +
    ggplot2::facet_grid(database ~ panel, switch = "y") +
    ggplot2::scale_colour_manual(values = c(Workflow = figureTaskColors[[1]], `Pooled lasso` = figureInk[["primary"]]),
      breaks = c("Workflow", "Pooled lasso")) +
    ggplot2::scale_fill_manual(values = c(Workflow = figureTaskColors[[1]], `Pooled lasso` = figureInk[["muted"]]),
      guide = "none") +
    ggplot2::scale_x_log10(breaks = breaks, labels = percent) +
    ggplot2::scale_y_log10(breaks = breaks, labels = percent) +
    ggplot2::coord_cartesian(xlim = limits, ylim = limits, expand = FALSE) +
    ggplot2::labs(x = "Predicted risk", y = "Observed risk (smoothed)",
      title = paste0(figureLabels$task[[rows$task[[1]]]], ", ", tolower(figureLabels$featureSet[[rows$featureSet[[1]]]])),
      caption = paste("Each curve spans the 1st to 99th percentile of predicted risk in the held-out database,",
        "with its pointwise 95% interval. Dashed: ideal calibration. Logarithmic axes.")) +
    figureTheme() +
    ggplot2::theme(text = ggplot2::element_text(size = 7.5, colour = figureInk[["primary"]]),
      strip.placement = "outside", strip.text.y.left = ggplot2::element_text(angle = 90),
      strip.text.x = ggplot2::element_text(hjust = 0.5),
      panel.border = ggplot2::element_rect(fill = NA, colour = figureInk[["grid"]]),
      panel.grid.major.y = ggplot2::element_line(colour = figureInk[["grid"]], linewidth = 0.3),
      panel.spacing = ggplot2::unit(4, "pt"),
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
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
