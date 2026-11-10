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
#' With `productionChecksDirectory`, the main-analysis rows of its `auc` table
#' add:
#' - `fidelity`: each workflow's AUC difference from pooled lasso against (A)
#'   the correlation of its linear predictor with pooled lasso's, with the
#'   median difference of unflagged fits within bands of that correlation, and
#'   (B) the share of pooled lasso's top 1% of predicted risk it also places
#'   there. The lasso developed in the held-out database is left out;
#' - `validity`: per task and held-out database, pooled lasso's observed AUC
#'   beside its model-based AUC (expected from case mix alone) and the AUC of
#'   the lasso developed in the held-out database;
#' - `cost`: per workflow, the median bytes and exchanges of one evaluation
#'   (DualAvg with and without penalty tuning) beside median top 1% overlap
#'   and AUC difference from pooled lasso, with interquartile ranges;
#' - `aucIntervals_<featureSet>`: every paired AUC difference from pooled
#'   lasso with its DeLong 95% interval, one panel per workflow.
#'
#' Each figure is written as PDF and PNG. Requires ggplot2.
#'
#' @param reportDirectory Output of [writeStudyReport()].
#' @param outputDirectory Directory for the figures; defaults to `reportDirectory`.
#' @param calibrationCurveDirectory Output of [runCalibrationCurveExport()], or
#'   `NULL` to skip the calibration figure.
#' @param productionChecksDirectory Output of [runProductionChecks()] (holding
#'   `export/auc.csv`), or `NULL` to skip the figures that need it.
#' @return The figure paths, invisibly.
#' @export
writeStudyFigures <- function(reportDirectory, outputDirectory = reportDirectory,
                              calibrationCurveDirectory = NULL, productionChecksDirectory = NULL) {
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
  if (!is.null(productionChecksDirectory)) {
    auc <- mainAnalysisChecks(utils::read.csv(file.path(productionChecksDirectory, "export", "auc.csv"),
      stringsAsFactors = FALSE), evaluations)
    communication <- utils::read.csv(file.path(reportDirectory, "communication.csv"), stringsAsFactors = FALSE)
    figures$fidelity <- list(plot = fidelityFigure(auc), width = 7.5, height = 4)
    figures$validity <- list(plot = validityFigure(auc), width = 7.5, height = 6.5)
    figures$cost <- list(plot = costFigure(communication, auc), width = 7.5, height = 3.4)
    for (featureSet in intersect(names(figureLabels$featureSet), auc$featureSet)) {
      figures[[paste0("aucIntervals_", featureSet)]] <- list(plot = aucIntervalFigure(auc, featureSet),
        width = 7.5, height = 8.5)
    }
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
    BiggestSiteLasso = "Largest database"),
  database = c(OPTUMEHR = "Optum EHR"),
  target = "Target-database lasso"
)

# Held-out databases as a factor of display names, ordered by database code
# (reversed for rows read from the top).
databaseFactor <- function(database, reverse = FALSE) {
  label <- function(x) {
    out <- figureLabels$database[x]
    ifelse(is.na(out), x, out)
  }
  codes <- sort(unique(database), decreasing = reverse)
  factor(unname(label(database)), levels = unname(label(codes)))
}

# The main-analysis rows of the production checks, with the report's flags.
mainAnalysisChecks <- function(auc, evaluations) {
  auc <- auc[startsWith(auc$resultsDirectory, paste0(auc$task, "_")), , drop = FALSE]
  key <- function(x) paste(x$task, x$featureSet, x$heldOutDatabase, x$method)
  auc$flagged <- key(auc) %in% key(evaluations[evaluations$flagged, , drop = FALSE])
  auc
}

signedNumber <- function(x, digits = 3) {
  ifelse(x > 0, sprintf("+%.*f", digits, x), sprintf("%.*f", digits, x))
}

# Places plots side by side with their panels aligned on the heights of the first.
alignedPlots <- function(plots) {
  grobs <- lapply(plots, ggplot2::ggplotGrob)
  do.call(cbind, c(grobs, size = "first"))
}

# A figure assembled from several plots. ggplot2 measures text when a plot is
# converted to a grob, so the conversion waits until the figure is drawn, on
# the device that draws it; otherwise labels measured with one device's fonts
# are clipped on another's.
deferredFigure <- function(build) {
  grid::gTree(build = build, cl = "deferredFigure")
}

#' @exportS3Method grid::makeContent
makeContent.deferredFigure <- function(x) {
  grid::setChildren(x, grid::gList(x$build()))
}

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
      plot.margin = ggplot2::margin(5, 16, 5, 5),
      plot.caption.position = "plot"
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
  background <- c(federated = "Other federated workflows", local = "Local-model baselines")
  groups <- c(highlighted, background)
  colours <- stats::setNames(c(figureInk[["primary"]], figureTaskColors[1:2], rep(figureInk[["muted"]], 2)), groups)
  localMethods <- c("LocalAvgLasso", "LocalEnsembleLasso", "BiggestSiteLasso")
  groupOf <- function(method) {
    factor(ifelse(method %in% names(highlighted), highlighted[method],
      ifelse(method %in% localMethods, background[["local"]], background[["federated"]])), levels = groups)
  }
  methods <- c("PooledLasso", names(figureLabels$method))
  rows <- evaluations[evaluations$featureSet == featureSet & evaluations$method %in% methods &
    is.finite(evaluations$calibrationInTheLarge), , drop = FALSE]
  rows$group <- groupOf(rows$method)
  rows$task <- factor(figureLabels$task[rows$task], levels = figureLabels$task)
  rows$database <- databaseFactor(rows$heldOutDatabase, reverse = TRUE)
  measures <- c("Calibration-in-the-large", "Calibration slope")
  long <- rbind(
    data.frame(rows[c("task", "database", "group", "flagged")], measure = measures[[1]],
      value = rows$calibrationInTheLarge, ideal = 0),
    data.frame(rows[c("task", "database", "group", "flagged")], measure = measures[[2]],
      value = rows$calibrationSlope, ideal = 1))
  long$measure <- factor(long$measure, levels = measures)
  long$shape <- paste(long$group, long$flagged)
  # Highlighted workflows sit above and below the row, the others on it.
  long$y <- as.numeric(long$database) + c(0, 0.22, -0.22, 0, 0)[as.integer(long$group)]
  long <- long[order(!long$group %in% background), ]
  # Filled markers, hollow when flagged; the local-model baselines are never flagged.
  shapeValues <- stats::setNames(c(3, 3, 16, 1, 17, 2, 15, 0, 4, 4), paste(rep(groups, each = 2), c(FALSE, TRUE)))
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
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2,
      override.aes = list(shape = c(3, 16, 17, 15, 4), size = 2.2))) +
    ggplot2::labs(x = NULL, y = NULL, tag = "A",
      caption = paste0("Held-out database per row. Ideal: 0 for calibration-in-the-large, 1 for the slope.\n",
        "Negative calibration-in-the-large means predicted risk above observed. Hollow: flagged fit.\n",
        "Other federated workflows: ODAL2, ADAP2, ADAP1. ",
        "Local-model baselines: coefficient average, probability ensemble, largest database.")) +
    figureTheme() +
    ggplot2::theme(strip.placement = "outside", strip.text.y.left = ggplot2::element_text(angle = 90),
      plot.caption = ggplot2::element_text(colour = figureInk[["secondary"]], hjust = 0))

  curve <- curves[curves$featureSet == featureSet & curves$task == curveTask &
    curves$method %in% names(highlighted), , drop = FALSE]
  curve$group <- groupOf(curve$method)
  curve$database <- databaseFactor(curve$heldOutDatabase)
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
  deferredFigure(function() grid::grobTree(
    grid::grobTree(ggplot2::ggplotGrob(measuresPlot), vp = grid::viewport(y = 1, height = 0.68, just = "top")),
    grid::grobTree(ggplot2::ggplotGrob(curvePlot), vp = grid::viewport(y = 0, height = 0.32, just = "bottom"))))
}

calibrationFigure <- function(rows) {
  .data <- ggplot2::.data
  labels <- c(PooledLasso = "Pooled lasso", figureLabels$method)
  labels <- labels[names(labels) %in% rows$method]
  rows$panel <- factor(labels[rows$method], levels = labels)
  rows$database <- databaseFactor(rows$heldOutDatabase)
  rows$series <- "Workflow"
  # Pooled lasso is drawn in every panel as the reference.
  pooled <- rows[rows$method == "PooledLasso", setdiff(names(rows), "panel"), drop = FALSE]
  pooled <- merge(pooled, data.frame(panel = factor(labels, levels = labels)))
  pooled$series <- "Pooled lasso"
  rows$series[rows$method == "PooledLasso"] <- "Pooled lasso"
  rows$fill <- rows$series
  # Axes cover nearly all curve points; the rest are counted in their panel,
  # by the side of the plotting area they lie beyond.
  range <- stats::quantile(c(rows$predicted, rows$observed), c(0.005, 0.995), names = FALSE)
  limits <- c(range[[1]] / 1.3, min(1, range[[2]] * 1.3))
  sides <- data.frame(above = rows$observed > limits[[2]], below = rows$observed < limits[[1]],
    right = rows$predicted > limits[[2]], left = rows$predicted < limits[[1]])
  offScale <- stats::aggregate(sides, rows[c("panel", "database")], sum)
  offScale <- offScale[rowSums(offScale[names(sides)]) > 0, , drop = FALSE]
  offScale$label <- vapply(seq_len(nrow(offScale)), function(i) {
    counts <- unlist(offScale[i, names(sides)])
    paste("Off scale:", paste(counts[counts > 0], names(counts)[counts > 0], collapse = ", "))
  }, character(1))
  breaks <- c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5)
  percent <- function(x) paste0(formatC(100 * x, format = "fg", digits = 2), "%")
  # Pooled lasso is drawn last and thinner, so it stays visible where a
  # workflow's curve lies on top of it.
  ggplot2::ggplot(rows, ggplot2::aes(x = .data$predicted, y = .data$observed)) +
    ggplot2::geom_abline(colour = figureInk[["muted"]], linewidth = 0.3, linetype = "22") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper, fill = .data$fill), alpha = 0.18) +
    ggplot2::geom_line(ggplot2::aes(colour = .data$series), linewidth = 0.8) +
    ggplot2::geom_line(data = pooled, ggplot2::aes(colour = .data$series), linewidth = 0.3) +
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
  local$row <- databaseFactor(local$heldOutDatabase, reverse = TRUE)
  sources <- evaluations[evaluations$method == "LocalSiteLasso" & evaluations$available, , drop = FALSE]
  sources$task <- factor(figureLabels$task[sources$task], levels = figureLabels$task)
  sources$featureSet <- factor(figureLabels$featureSet[sources$featureSet], levels = figureLabels$featureSet)
  sources$row <- factor(databaseFactor(sources$heldOutDatabase), levels = levels(local$row))
  models <- c("Coefficient average", "Probability ensemble", "Pooled lasso", "Largest training database")
  aggregates <- rbind(
    data.frame(local[c("task", "featureSet", "row")], model = models[[1]], auc = local$coefficientAverageAuc),
    data.frame(local[c("task", "featureSet", "row")], model = models[[2]], auc = local$probabilityEnsembleAuc),
    data.frame(local[c("task", "featureSet", "row")], model = models[[3]], auc = local$pooledAuc),
    data.frame(local[c("task", "featureSet", "row")], model = models[[4]], auc = local$largestSiteAuc))
  aggregates <- aggregates[is.finite(aggregates$auc), , drop = FALSE]
  aggregates$model <- factor(aggregates$model, levels = models)
  # Coefficient average above and probability ensemble below the row, so
  # near-equal AUCs stay visible; the largest database's model is ringed on the row.
  aggregates$nudge <- c(0.2, -0.2, 0, 0)[as.integer(aggregates$model)]
  ggplot2::ggplot(local, ggplot2::aes(y = .data$row)) +
    ggplot2::geom_segment(ggplot2::aes(x = .data$minSourceAuc, xend = .data$maxSourceAuc, yend = .data$row),
      colour = figureInk[["grid"]], linewidth = 2.2, lineend = "round") +
    ggplot2::geom_point(data = sources, ggplot2::aes(x = .data$auc), colour = figureInk[["muted"]], size = 1.1) +
    ggplot2::geom_point(data = aggregates, ggplot2::aes(x = .data$auc, y = as.numeric(.data$row) + .data$nudge,
      colour = .data$model, shape = .data$model, size = .data$model), stroke = 0.7) +
    ggplot2::facet_grid(task ~ featureSet, scales = "free_y", space = "free_y", switch = "y") +
    ggplot2::scale_colour_manual(values = stats::setNames(c(figureTaskColors[1:2], figureInk[["primary"]],
      figureInk[["secondary"]]), models)) +
    ggplot2::scale_shape_manual(values = stats::setNames(c(16, 17, 3, 1), models)) +
    ggplot2::scale_size_manual(values = stats::setNames(c(2, 2, 2, 2.8), models)) +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2)) +
    ggplot2::labs(x = "AUC in the held-out database", y = NULL,
      caption = paste("Grey: models fitted in one training database (range and each model).",
        "Ring: the model of the largest training database.")) +
    figureTheme() +
    ggplot2::theme(strip.placement = "outside", strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 0),
      plot.caption = ggplot2::element_text(colour = figureInk[["secondary"]], hjust = 0))
}

# A caption under a combined figure, in the style of the single-plot captions.
withCaption <- function(grob, caption) {
  lines <- length(strsplit(caption, "\n", fixed = TRUE)[[1]])
  space <- grid::unit(lines * 9 + 8, "pt")
  grid::grobTree(
    grid::grobTree(grob, vp = grid::viewport(y = 1, height = grid::unit(1, "npc") - space, just = "top")),
    grid::textGrob(caption, x = grid::unit(5, "pt"), y = grid::unit(4, "pt"), just = c("left", "bottom"),
      gp = grid::gpar(fontsize = 7.2, col = figureInk[["secondary"]], lineheight = 1.1)))
}

# The eight workflows in three families, since only three categorical
# colours stay distinguishable in a scatter.
workflowFamily <- function(method) {
  families <- c("DualAvg, ADAP2, ODAL2", "ODAL1, ADAP1", "Local-model baselines")
  slot <- c(DualAvg = 1L, ADAP2 = 1L, ODAL2 = 1L, ODAL1 = 2L, ADAP1 = 2L, LocalAvgLasso = 3L,
    LocalEnsembleLasso = 3L, BiggestSiteLasso = 3L)[method]
  factor(families[slot], levels = families)
}

# Fidelity against performance: does a workflow that reproduces pooled lasso
# more closely also do better in the held-out database?
fidelityFigure <- function(auc) {
  .data <- ggplot2::.data
  rows <- auc[auc$method %in% names(figureLabels$method) & is.finite(auc$aucDifference) &
    is.finite(auc$linearPredictorPearson), , drop = FALSE]
  rows$family <- workflowFamily(rows$method)
  rows$shape <- paste(as.integer(rows$family), rows$flagged)
  shapes <- c(`1 FALSE` = 16, `1 TRUE` = 1, `2 FALSE` = 17, `2 TRUE` = 2, `3 FALSE` = 15, `3 TRUE` = 0)
  # Closeness is -log10(1 - correlation), so each step to the right is one more 9.
  rows$closeness <- -log10(pmax(1 - rows$linearPredictorPearson, 1e-6))
  cuts <- c(0.90, 0.95, 0.98, 0.995)
  unflagged <- rows[!rows$flagged, , drop = FALSE]
  band <- cut(unflagged$linearPredictorPearson, c(-Inf, cuts, Inf))
  edges <- -log10(1 - cuts)
  medians <- data.frame(median = as.numeric(tapply(unflagged$aucDifference, band, stats::median)),
    fits = as.numeric(table(band)), from = c(min(unflagged$closeness), edges),
    to = c(edges, max(unflagged$closeness)))
  medians <- medians[medians$fits > 0, , drop = FALSE]
  limits <- range(rows$aucDifference)
  breaks <- 1:5
  points <- list(
    ggplot2::geom_hline(yintercept = 0, colour = figureInk[["muted"]], linewidth = 0.4),
    ggplot2::geom_point(ggplot2::aes(colour = .data$family, shape = .data$shape), size = 1.5, stroke = 0.55,
      alpha = 0.85),
    ggplot2::scale_colour_manual(values = stats::setNames(figureTaskColors, levels(rows$family)), drop = FALSE),
    ggplot2::scale_shape_manual(values = shapes, guide = "none"),
    ggplot2::scale_y_continuous(labels = signedNumber),
    ggplot2::coord_cartesian(ylim = limits),
    figureTheme(),
    ggplot2::theme(panel.grid.major.y = ggplot2::element_line(colour = figureInk[["grid"]], linewidth = 0.3),
      plot.tag = ggplot2::element_text(face = "bold")))
  closeness <- ggplot2::ggplot(rows, ggplot2::aes(x = .data$closeness, y = .data$aucDifference)) +
    points +
    ggplot2::geom_segment(data = medians, ggplot2::aes(x = .data$from, xend = .data$to, y = .data$median,
      yend = .data$median), colour = figureInk[["primary"]], linewidth = 0.9, inherit.aes = FALSE) +
    ggplot2::scale_x_continuous(breaks = breaks, labels = as.character(1 - 10^-breaks)) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(shape = c(16, 17, 15), size = 2.2,
      alpha = 1))) +
    ggplot2::labs(x = "Correlation with pooled lasso's linear predictor", y = "AUC difference from pooled lasso",
      tag = "A")
  overlap <- ggplot2::ggplot(rows[is.finite(rows$top1Overlap), , drop = FALSE],
    ggplot2::aes(x = .data$top1Overlap, y = .data$aucDifference)) +
    points +
    ggplot2::scale_x_continuous(labels = function(x) paste0(round(100 * x), "%")) +
    ggplot2::labs(x = "Overlap with pooled lasso's top 1% of predicted risk", y = NULL,
      tag = "B") +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), legend.position = "none")
  deferredFigure(function() withCaption(alignedPlots(list(closeness, overlap)), paste0(
    "Each marker is one workflow in one task, held-out database and predictor set; hollow: flagged fit. ",
    "Black lines (A): median difference of the unflagged fits\nwith correlation up to 0.90, 0.95, 0.98, 0.995 ",
    "and above. Local-model baselines: coefficient average, probability ensemble, largest database. ",
    "The axis in A is\nlogarithmic in 1 minus the correlation.")))
}

# Case mix against validity: pooled lasso's model-based AUC (what its risks
# would give if they were correct in the held-out database) beside its
# observed AUC, and the AUC of a lasso developed in that database.
validityFigure <- function(auc) {
  .data <- ggplot2::.data
  models <- c("Pooled lasso, expected from case mix", "Pooled lasso, observed", figureLabels$target)
  keys <- c("task", "featureSet", "heldOutDatabase")
  pooled <- auc[auc$method == "PooledLasso", , drop = FALSE]
  target <- auc[auc$method == "TargetLasso" & is.finite(auc$auc), , drop = FALSE]
  levels <- levels(databaseFactor(pooled$heldOutDatabase, reverse = TRUE))
  labelled <- function(x) {
    x$task <- factor(figureLabels$task[x$task], levels = figureLabels$task)
    x$featureSet <- factor(figureLabels$featureSet[x$featureSet], levels = figureLabels$featureSet)
    x$row <- match(as.character(databaseFactor(x$heldOutDatabase)), levels)
    x
  }
  points <- labelled(rbind(
    data.frame(pooled[keys], model = models[[1]], auc = pooled$modelBasedAuc),
    data.frame(pooled[keys], model = models[[2]], auc = pooled$auc),
    data.frame(target[keys], model = models[[3]], auc = target$auc)))
  points$model <- factor(points$model, levels = models)
  # The target-database lasso sits below the row, so near-equal AUCs stay visible.
  points$y <- points$row + c(0, 0, -0.24)[as.integer(points$model)]
  gaps <- labelled(data.frame(pooled[keys], expected = pooled$modelBasedAuc, observed = pooled$auc))
  ggplot2::ggplot(points, ggplot2::aes(x = .data$auc, y = .data$y)) +
    ggplot2::geom_segment(data = gaps, ggplot2::aes(x = .data$expected, xend = .data$observed, y = .data$row,
      yend = .data$row), colour = figureInk[["muted"]], linewidth = 0.6, inherit.aes = FALSE) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$model, shape = .data$model), size = 2, stroke = 0.7) +
    ggplot2::facet_grid(task ~ featureSet, switch = "y") +
    ggplot2::scale_y_continuous(breaks = seq_along(levels), labels = levels) +
    ggplot2::scale_colour_manual(values = stats::setNames(c(figureTaskColors[[1]], figureInk[["primary"]],
      figureTaskColors[[2]]), models)) +
    ggplot2::scale_shape_manual(values = stats::setNames(c(1, 3, 17), models)) +
    ggplot2::labs(x = "AUC in the held-out database", y = NULL,
      caption = paste0("Expected: the AUC pooled lasso would reach if its predicted risks were the true risks ",
        "in the held-out database (model-based AUC),\nwhich depends on case mix only. Grey line: from expected to ",
        "observed. Target-database lasso: developed in the held-out database by\nfive-fold cross-validation ",
        "on pooled lasso's retained predictors, scored out of fold.")) +
    figureTheme() +
    ggplot2::theme(strip.placement = "outside", strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 0),
      plot.caption = ggplot2::element_text(colour = figureInk[["secondary"]], hjust = 0))
}

# Paired AUC differences from pooled lasso with their DeLong 95% intervals,
# one panel per workflow, rows grouped by task.
aucIntervalFigure <- function(auc, featureSet, limits = c(-0.11, 0.085)) {
  .data <- ggplot2::.data
  labels <- c(figureLabels$method, TargetLasso = figureLabels$target)
  rows <- auc[auc$featureSet == featureSet & auc$method %in% names(labels) & is.finite(auc$aucDifference), ,
    drop = FALSE]
  rows$panel <- factor(labels[rows$method], levels = labels[names(labels) %in% rows$method])
  databases <- sort(unique(auc$heldOutDatabase))
  step <- length(databases) + 1
  rows$y <- -((match(rows$task, names(figureLabels$task)) - 1) * step + match(rows$heldOutDatabase, databases))
  breaks <- -(rep((seq_along(figureLabels$task) - 1) * step, each = length(databases)) + seq_along(databases))
  rows$task <- factor(figureLabels$task[rows$task], levels = figureLabels$task)
  rows$shape <- paste(as.integer(rows$task), rows$flagged)
  shapes <- stats::setNames(c(figureTaskShapes, figureTaskShapesHollow),
    paste(seq_along(figureTaskShapes), rep(c(FALSE, TRUE), each = length(figureTaskShapes))))
  off <- rows$aucDifference < limits[[1]] | rows$aucDifference > limits[[2]]
  beyond <- rows[off, , drop = FALSE]
  beyond$direction <- sign(beyond$aucDifference)
  beyond$edge <- ifelse(beyond$direction < 0, limits[[1]], limits[[2]])
  beyond$hjust <- ifelse(beyond$direction < 0, 0, 1)
  shown <- rows[!off, , drop = FALSE]
  ggplot2::ggplot(shown, ggplot2::aes(y = .data$y, colour = .data$task)) +
    ggplot2::geom_vline(xintercept = 0, colour = figureInk[["muted"]], linewidth = 0.4) +
    ggplot2::geom_segment(ggplot2::aes(x = .data$differenceLower, xend = .data$differenceUpper, yend = .data$y),
      linewidth = 0.45) +
    ggplot2::geom_point(ggplot2::aes(x = .data$aucDifference, shape = .data$shape), size = 1.4, stroke = 0.5,
      fill = "white") +
    ggplot2::geom_segment(data = beyond, ggplot2::aes(x = .data$edge - .data$direction * 0.014, xend = .data$edge,
      yend = .data$y), linewidth = 0.45, arrow = grid::arrow(length = grid::unit(3, "pt"), type = "closed"),
      show.legend = FALSE) +
    ggplot2::geom_text(data = beyond, ggplot2::aes(x = .data$edge - .data$direction * 0.017,
      label = signedNumber(.data$aucDifference, 2), hjust = .data$hjust), size = 2.2,
      colour = figureInk[["secondary"]], show.legend = FALSE) +
    ggplot2::facet_wrap(~ panel, ncol = 3) +
    ggplot2::scale_y_continuous(breaks = breaks,
      labels = rep(levels(databaseFactor(databases)), length(figureLabels$task))) +
    ggplot2::scale_x_continuous(labels = function(x) signedNumber(x, 2)) +
    ggplot2::scale_colour_manual(values = stats::setNames(figureTaskColors, figureLabels$task), drop = FALSE) +
    ggplot2::scale_shape_manual(values = shapes, guide = "none") +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(shape = figureTaskShapes, size = 2.2,
      linetype = 0))) +
    ggplot2::coord_cartesian(xlim = limits) +
    ggplot2::labs(x = "AUC difference from pooled lasso in the held-out database, with 95% interval", y = NULL,
      title = figureLabels$featureSet[[featureSet]],
      caption = paste0("Paired DeLong intervals. Rows are held-out databases, grouped by task. Hollow: flagged fit. ",
        "Arrow: estimate beyond the axis, with its value.\nA missing row is a fit that failed.")) +
    figureTheme() +
    ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA, colour = figureInk[["grid"]]),
      axis.text.y = ggplot2::element_text(size = 6.5), plot.title = ggplot2::element_text(face = "bold"),
      plot.caption = ggplot2::element_text(colour = figureInk[["secondary"]], hjust = 0))
}

# Communication cost against agreement and performance, one row per workflow.
costFigure <- function(communication, auc) {
  .data <- ggplot2::.data
  methods <- figureLabels$method
  sets <- figureLabels$featureSet
  positioned <- function(x) {
    x$row <- factor(methods[x$method], levels = rev(methods))
    x$featureSet <- factor(sets[x$featureSet], levels = sets)
    x$y <- as.numeric(x$row) + c(0.16, -0.16)[as.integer(x$featureSet)]
    x
  }
  zero <- function(x) ifelse(is.na(x), 0, x)
  cost <- positioned(communication[communication$method %in% names(methods), , drop = FALSE])
  cost$bytes <- cost$medianBytes + zero(cost$medianTuningBytes)
  cost$exchanges <- cost$medianExchanges + zero(cost$medianTuningExchanges)
  tuned <- cost[cost$bytes > cost$medianBytes, , drop = FALSE]
  rows <- auc[auc$method %in% names(methods), , drop = FALSE]
  summarised <- function(column) {
    groups <- split(rows[[column]], list(rows$method, rows$featureSet), drop = TRUE, sep = "|")
    groups <- lapply(groups, function(v) v[is.finite(v)])
    parts <- do.call(rbind, strsplit(names(groups), "|", fixed = TRUE))
    positioned(data.frame(method = parts[, 1], featureSet = parts[, 2],
      median = vapply(groups, stats::median, numeric(1)),
      lower = vapply(groups, stats::quantile, numeric(1), 0.25, names = FALSE),
      upper = vapply(groups, stats::quantile, numeric(1), 0.75, names = FALSE), stringsAsFactors = FALSE))
  }
  colours <- stats::setNames(figureTaskColors[1:2], sets)
  shapes <- stats::setNames(c(16, 17), sets)
  levels <- levels(cost$row)
  panel <- function(data, title, first = FALSE) {
    ggplot2::ggplot(data, ggplot2::aes(y = .data$y, colour = .data$featureSet, shape = .data$featureSet)) +
      ggplot2::scale_y_continuous(breaks = seq_along(levels), labels = levels, limits = c(0.6, length(levels) + 0.4)) +
      ggplot2::scale_colour_manual(values = colours, drop = FALSE) +
      ggplot2::scale_shape_manual(values = shapes, drop = FALSE) +
      ggplot2::labs(x = NULL, y = NULL, title = title) +
      figureTheme() +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 8, face = "bold"),
        panel.grid.major.y = ggplot2::element_line(colour = figureInk[["grid"]], linewidth = 0.3),
        axis.text.y = if (first) ggplot2::element_text(colour = figureInk[["secondary"]]) else ggplot2::element_blank(),
        legend.position = if (first) "top" else "none", plot.margin = ggplot2::margin(5, 8, 5, 2))
  }
  logLabels <- function(units) function(x) {
    out <- character(length(x))
    for (i in seq_along(x)) {
      if (is.na(x[[i]])) next
      unit <- max(which(units <= x[[i]] * (1 + 1e-9)), 1)
      out[[i]] <- paste0(format(x[[i]] / units[[unit]], big.mark = ",", scientific = FALSE, trim = TRUE),
        names(units)[[unit]])
    }
    out
  }
  bytes <- panel(cost, "Data exchanged", first = TRUE) +
    ggplot2::geom_segment(data = tuned, ggplot2::aes(x = .data$medianBytes, xend = .data$bytes, yend = .data$y),
      linewidth = 0.4, show.legend = FALSE) +
    ggplot2::geom_point(data = tuned, ggplot2::aes(x = .data$medianBytes), shape = 1, size = 1.8, stroke = 0.6,
      show.legend = FALSE) +
    ggplot2::geom_point(ggplot2::aes(x = .data$bytes), size = 1.8) +
    ggplot2::scale_x_log10(breaks = 10^seq(2, 10, by = 2),
      labels = logLabels(c(` B` = 1, ` kB` = 1e3, ` MB` = 1e6, ` GB` = 1e9)))
  exchanges <- panel(cost, "Exchanges") +
    ggplot2::geom_segment(data = tuned, ggplot2::aes(x = .data$medianExchanges, xend = .data$exchanges,
      yend = .data$y), linewidth = 0.4) +
    ggplot2::geom_point(data = tuned, ggplot2::aes(x = .data$medianExchanges), shape = 1, size = 1.8,
      stroke = 0.6) +
    ggplot2::geom_point(ggplot2::aes(x = .data$exchanges), size = 1.8) +
    ggplot2::scale_x_log10(breaks = 10^seq(0, 6, by = 2), labels = logLabels(c(1, ` thousand` = 1e3)))
  interval <- function(data, title) {
    panel(data, title) +
      ggplot2::geom_segment(ggplot2::aes(x = .data$lower, xend = .data$upper, yend = .data$y), linewidth = 0.45) +
      ggplot2::geom_point(ggplot2::aes(x = .data$median), size = 1.8)
  }
  overlap <- interval(summarised("top1Overlap"), "Top 1% overlap") +
    ggplot2::scale_x_continuous(labels = function(x) paste0(round(100 * x), "%"))
  difference <- interval(summarised("aucDifference"), "AUC difference") +
    ggplot2::geom_vline(xintercept = 0, colour = figureInk[["muted"]], linewidth = 0.4) +
    ggplot2::scale_x_continuous(breaks = seq(-0.1, 0.1, by = 0.01), labels = function(x) signedNumber(x, 2))
  deferredFigure(function() withCaption(alignedPlots(list(bytes, exchanges, overlap, difference)), paste0(
    "Medians over the 15 evaluations of each predictor set; bars: interquartile range. Overlap and AUC difference ",
    "are with pooled lasso,\nthe reference, which pools patient-level data. Data: serialized messages in both ",
    "directions, as counted by the federated protocols or modeled for the\nlocal-model baselines (one model upload ",
    "per database). Hollow: DualAvg without its federated penalty tuning; the other workflows tune within one ",
    "database.")))
}
