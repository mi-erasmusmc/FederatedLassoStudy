simulatedLogistic <- function(n = 3000, seed = 3) {
  set.seed(seed)
  age <- runif(n, 0.4, 0.9)
  x <- cbind(1, age, Matrix::rsparsematrix(n, 6, density = 0.2, rand.x = function(k) rep(1, k)))
  y <- stats::rbinom(n, 1, stats::plogis(-3 + 1.5 * age + as.numeric(x[, 3:8] %*% c(0.8, -0.5, 0.4, 0, 0, 0.3))))
  list(x = methods::as(x, "CsparseMatrix"), y = y)
}
control <- list(tolerance = 1e-12, maxIterations = 20000L, threads = 1L)

test_that("the lasso refit on centred age is the uncentred solution, and meets the optimality conditions", {
  d <- simulatedLogistic()
  lambda <- 0.002
  plain <- cyclopsFit(d$x, d$y, "laplace", 1L, control, lambda)$w
  center <- mean(d$x[, 2])
  centred <- d$x
  centred[, 2] <- d$x[, 2] - center
  w <- cyclopsFit(centred, d$y, "laplace", 1L, control, lambda)$w
  w[1] <- w[1] - w[2] * center
  expect_equal(w, plain, tolerance = 1e-5)
  diagnostics <- lassoDiagnostics(d$x, d$y, w, lambda)
  expect_lt(diagnostics$kktMax, 1e-6)
  expect_equal(diagnostics$objective, meanLogLoss(d$x, d$y, w) + lambda * sum(abs(w[-1])))
  # A perturbed solution has a larger objective and violates the conditions.
  worse <- lassoDiagnostics(d$x, d$y, w + c(0.05, -0.1, rep(0, 6)), lambda)
  expect_gt(worse$objective, diagnostics$objective)
  expect_gt(worse$kktMax, 1e-4)
})

test_that("flipping a prevalent binary predictor and centring age give the same lasso and ridge solutions", {
  d <- simulatedLogistic()
  d$x[, 3] <- 1 - d$x[, 3]
  d$x <- Matrix::drop0(d$x)
  expect_true(binaryColumn(d$x, 3L) && !binaryColumn(d$x, 2L) && mean(d$x[, 3]) > 0.5)
  center <- mean(d$x[, 2])
  design <- reparametrize(d$x, 2L, center, 3L)
  expect_equal(as.numeric(design[, 3]), 1 - as.numeric(d$x[, 3]))
  expect_equal(mean(design[, 2]), 0)
  for (prior in list("laplace", list("normal", 0.5))) {
    plain <- cyclopsFit(d$x, d$y, prior, 1L, control, 0.002)$w
    w <- restoreCoefficients(cyclopsFit(design, d$y, prior, 1L, control, 0.002)$w, 2L, center, 3L)
    expect_equal(w, plain, tolerance = 1e-5)
  }
})

test_that("an unpenalized fit and its age standard error match glm", {
  d <- simulatedLogistic()
  fit <- cyclopsFit(d$x, d$y, "none", integer(), control)
  reference <- stats::glm(d$y ~ as.matrix(d$x) - 1, family = stats::binomial())
  expect_equal(fit$w, unname(stats::coef(reference)), tolerance = 1e-6)
  se <- ageStandardError(d$x, d$y, fit$w, rep(0, ncol(d$x)), 2L)
  expect_equal(se, unname(sqrt(diag(stats::vcov(reference)))[2]), tolerance = 1e-5)
  expect_equal(meanLogLoss(d$x, d$y, fit$w), -as.numeric(stats::logLik(reference)) / length(d$y))
})

test_that("a weak normal prior leaves age nearly unchanged and keeps it unpenalized", {
  d <- simulatedLogistic()
  free <- 1:2
  normal <- cyclopsFit(d$x, d$y, list("normal", 10), free, control)$w
  none <- cyclopsFit(d$x, d$y, "none", integer(), control)$w
  expect_equal(normal[2], none[2], tolerance = 0.02)
  strong <- cyclopsFit(d$x, d$y, list("normal", 1e-4), free, control)$w
  expect_true(all(abs(strong[3:8]) < 0.01))
  expect_gt(abs(strong[2]), 0.5)
})

test_that("database summaries suppress means over too few entries", {
  pooled <- list(models = list(list(coefficients = data.frame(matrixColumn = 1:2, coefficientSharedScale = 0,
    isIntercept = c(TRUE, FALSE), covariateId = c(NA, 1002)))))
  matrices <- list(A = list(xMatrix = cbind(1, c(0.5, 0.6, 0.7, 0.8)), yLabels = c(0, 0, 0, 1)))
  out <- databaseAges(matrices, pooled, "T", 2L)
  expect_equal(out$ageMean, 0.65)
  expect_true(is.na(out$ageMeanOutcome))
  expect_equal(out$ageMeanNoOutcome, 0.6)
})
