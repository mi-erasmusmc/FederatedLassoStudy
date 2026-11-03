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

test_that("the model-based AUC matches its pairwise definition", {
  set.seed(2)
  p <- c(round(runif(40), 1), 0.5, 0.5)
  pairs <- expand.grid(i = seq_along(p), j = seq_along(p))
  pairs <- pairs[pairs$i != pairs$j, ]
  weight <- p[pairs$i] * (1 - p[pairs$j])
  score <- (p[pairs$i] > p[pairs$j]) + 0.5 * (p[pairs$i] == p[pairs$j])
  expect_equal(modelBasedAuc(p), sum(weight * score) / sum(weight))
  expect_equal(modelBasedAuc(rep(0.2, 10)), 0.5)
})

test_that("top-risk sets have the requested size and break ties by position", {
  expect_identical(topRisk(c(0.1, 0.9, 0.5, 0.9), 0.5), c(2L, 4L))
  expect_identical(topRisk(c(0.1, 0.2), 0.01), 2L)
})

test_that("site labels are matched to the population by subject and index date", {
  directory <- tempfile()
  dir.create(file.path(directory, "13524"), recursive = TRUE)
  write.csv(data.frame(subjectId = c("20000000001", "20000000002", "20000000003"),
    cohortStartDate = c("2010-01-01", "2011-01-01", "2012-01-01"), label = c("TX", NA, "TX")),
    file.path(directory, "13524", "CCAE.csv"), row.names = FALSE)
  population <- data.frame(subjectId = c(20000000003, 20000000002, 20000000001), targetId = 13524,
    cohortStartDate = as.Date(c("2012-01-01", "2011-01-01", "2010-01-01")))
  expect_identical(siteLabels(directory, population, "CCAE"), c(1L, 2L, 1L))
  expect_null(siteLabels(directory, population, "MDCD"))
  population$cohortStartDate <- population$cohortStartDate + 1
  expect_error(siteLabels(directory, population, "CCAE"), "have a site label")
})
