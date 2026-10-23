library(testthat)
library(Rssa)
source(system.file("extdata", "common.test.methods.R", package = "Rssa"))
context("Tensor")

test_that("n-unfolding is correct", {
  tens <- c(1, 2) %o% c(0, 3) %o% c(1, 0.5)
  unfolds <- list()
  unfolds[[1]] <- matrix(c(0, 0, 3, 6, 0, 0, 1.5, 3), nrow = 2)
  unfolds[[2]] <- matrix(c(0, 3, 0, 6, 0, 1.5, 0, 3), nrow = 2)
  unfolds[[3]] <- matrix(c(0, 0, 0, 0, 3, 1.5, 6, 3), nrow = 2)
  for (m in 1:3) {
    expect_true(identical(n.unfold(tens, m), unfolds[[m]]))
  }
})

test_that("n-folding is correct", {
  tens <- c(1, 2) %o% c(0, 3) %o% c(1, 0.5)
  modes <- c(2, 2, 2)
  unfolds <- list()
  unfolds[[1]] <- matrix(c(0, 0, 3, 6, 0, 0, 1.5, 3), nrow = 2)
  unfolds[[2]] <- matrix(c(0, 3, 0, 6, 0, 1.5, 0, 3), nrow = 2)
  unfolds[[3]] <- matrix(c(0, 0, 0, 0, 3, 1.5, 6, 3), nrow = 2)
  for (m in 1:3) {
    expect_true(identical(n.fold(unfolds[[m]], m, modes), tens))
  }
})

test_that("reunfolding is correct", {
  set.seed(1)
  tens <- rnorm(3) %o% rnorm(4) %o% rnorm(5)
  modes <- c(3, 4, 5)
  for (m1 in 1:3) {
    mat <- n.unfold(tens, m1)
    for (m2 in 1:3) {
      expect_true(identical(.reunfold(mat, modes, m1, m2), n.unfold(tens, m2)))
    }
  }
})

test_that("Same-mode matrices are correctly pre-multiplying in n-mode product list variant", {
  set.seed(1)
  tens <- rnorm(3) %o% rnorm(4) %o% rnorm(5)
  mats <- list(
    rnorm(6) %o% rnorm(3), # 1st mode
    rnorm(3) %o% rnorm(4), # 2nd mode
    rnorm(4) %o% rnorm(6), # 1st mode
    rnorm(6) %o% rnorm(4), # 1st mode
    rnorm(2) %o% rnorm(3)  # 2nd mode
  )
  
  mats.pm <- list(
    mats[[4]] %*% mats[[3]] %*% mats[[1]],
    mats[[5]] %*% mats[[2]]
  )
  prod1 <- n.mode.prod.list(tens, mats, c(1, 2, 1, 1, 2))
  prod2 <- n.mode.prod.list(tens, mats.pm, c(1, 2))

  expect_true(mean((prod1 - prod2)^2) < sqrt(.Machine$double.eps))
})

test_that("Reconstruction of full-ranked HOSVD equals to initial tensor", {
  set.seed(1)
  F <- rnorm(10)
  tens <- new.trajtens(F, c(3, 4, 5))
  unfolds <- new.unfolds(F, c(3, 4, 5))
  svd.methods <- c("svd", "eigen", "primme")
  
  for (method in svd.methods) {
    dec <- hosvd(tens, unfolds, ranks = c(3, 4, 5), svd.method = method, calc.est = TRUE)
    
    expect_true(dec$fnorm_resid < sqrt(.Machine$double.eps))
  }
})

test_that("Reconstruction of full-ranked complex HOSVD equals to initial tensor", {
  set.seed(1)
  F <- rnorm(10) + 1i * rnorm(10)
  tens <- new.trajtens(F, c(3, 4, 5))
  unfolds <- new.unfolds(F, c(3, 4, 5))
  svd.methods <- c("svd", "eigen", "primme")
  
  for (method in svd.methods) {
    dec <- hosvd(tens, unfolds, ranks = c(3, 4, 5), svd.method = method, calc.est = TRUE)
    
    expect_true(dec$fnorm_resid < sqrt(.Machine$double.eps))
  }
})


test_that("Reconstruction of non-full-ranked HOSVD equals to initial tensor", {
  F <- exp(0.5 + 1:10 / 12) + 1
  tens <- new.trajtens(F, c(3, 4, 5))
  unfolds <- new.unfolds(F, c(3, 4, 5))
  svd.methods <- c("svd", "eigen", "primme")
  
  for (method in svd.methods) {
    dec <- hosvd(tens, unfolds, ranks = 2, svd.method = method, calc.est = TRUE)
    
    expect_true(dec$fnorm_resid < sqrt(.Machine$double.eps))
  }
})

test_that("Reconstruction of non-full-ranked complex HOSVD equals to initial tensor", {
  F <- exp(0.1 * 1:10 + 2i * pi * 1:10 / 2) + 1
  tens <- new.trajtens(F, c(3, 4, 5))
  unfolds <- new.unfolds(F, c(3, 4, 5))
  svd.methods <- c("svd", "eigen", "primme")
  
  for (method in svd.methods) {
    dec <- hosvd(tens, unfolds, ranks = 2, svd.method = method, calc.est = TRUE)
    
    expect_true(dec$fnorm_resid < sqrt(.Machine$double.eps))
  }
})
