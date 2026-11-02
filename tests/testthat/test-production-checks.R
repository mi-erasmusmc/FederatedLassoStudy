test_that("DeLong AUCs and covariance match the definition", {
  set.seed(1)
  y <- rbinom(60, 1, 0.3)
  scores <- cbind(a = y + rnorm(60), b = round(y + rnorm(60, sd = 2)), c = 0)
  delong <- delongAuc(y, scores)
  x <- scores[y == 1, , drop = FALSE]
  z <- scores[y == 0, , drop = FALSE]
  psi <- lapply(seq_len(ncol(scores)), function(k) outer(x[, k], z[, k], function(a, b) (a > b) + 0.5 * (a == b)))
  expect_equal(delong$auc, vapply(psi, mean, numeric(1)), ignore_attr = TRUE)
  v10 <- sapply(psi, rowMeans)
  v01 <- sapply(psi, colMeans)
  expect_equal(delong$covariance, cov(v10) / nrow(x) + cov(v01) / nrow(z), ignore_attr = TRUE)
  expect_equal(unname(delong$auc[3]), 0.5)
  expect_equal(unname(delong$covariance[3, 3]), 0)
  expect_error(delongAuc(c(1, 0, 0), cbind(1:3)), "at least two")
})

test_that("prevalence suppresses small counts of patients with or without the covariate", {
  counts <- list(present = c("(Intercept)" = 100, "1160049" = 97, "1165049" = 40), n = 100)
  prevalence <- covariatePrevalence(counts, c(1160049, 1165049, 99), "T", "S", 5L)
  expect_identical(prevalence$suppressed, c(TRUE, FALSE, FALSE))
  expect_identical(prevalence$withCovariate, c(NA, 40, NA))
  expect_identical(prevalence$percent, c(NA, 40, NA))
})
