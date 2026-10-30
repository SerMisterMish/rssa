library(testthat)
library(Rssa)
source(system.file("extdata", "common.test.methods.R", package = "Rssa"))
context("Tensor SSA")

test_that("3D-TSSA reconstructs finite-rank series exactly for all tensor decompositions and svd.methods", {
  N <- 13; xs <- 1:N; L <- c(4, 5)
  
  tens.decomps <- c("hutd", "hutd+hooi", "hosvd", "hooi")
  svd.methods <- c("svd", "eigen", "primme")
  
  s1 <- sin(2 * pi * xs / 3)
  r1 <- 2
  s2 <- exp(xs / N + 2i * pi * xs / 3)
  r2 <- 1
  
  for (tens.decomp in tens.decomps) {
    for (svd.method in svd.methods) {
      tssa1 <- ssa(s1, L, svd.method = svd.method, tensor.decomposition = tens.decomp)
      rec1 <- reconstruct(tssa1, groups = list(1:r1))
      expect_equal(rec1$F1, s1, tolerance = sqrt(.Machine$double.eps))
      
      tssa2 <- ssa(s2, L, kind = "tssa", svd.method = svd.method, tensor.decomposition = tens.decomp)
      rec2 <- reconstruct(tssa2, groups = list(1:r2))
      expect_equal(rec2$F1, s2, tolerance = sqrt(.Machine$double.eps))
    }
  }
})

test_that("4D-TSSA reconstructs finite-rank series exactly for all tensor decompositions and svd.methods", {
  N <- 19; xs <- 1:N; L <- c(4, 5, 6)
  
  tens.decomps <- c("hutd", "hutd+hooi", "hosvd", "hooi")
  svd.methods <- c("svd", "eigen", "primme")
  
  s1 <- sin(2 * pi * xs / 3)
  r1 <- 2
  s2 <- exp(xs / N + 2i * pi * xs / 3)
  r2 <- 1
  
  for (tens.decomp in tens.decomps) {
    for (svd.method in svd.methods) {
      tssa1 <- ssa(s1, L, svd.method = svd.method, tensor.decomposition = tens.decomp)
      rec1 <- reconstruct(tssa1, groups = list(1:r1))
      expect_equal(rec1$F1, s1, tolerance = sqrt(.Machine$double.eps))
      
      tssa2 <- ssa(s2, L, kind = "tssa", svd.method = svd.method, tensor.decomposition = tens.decomp)
      rec2 <- reconstruct(tssa2, groups = list(1:r2))
      expect_equal(rec2$F1, s2, tolerance = sqrt(.Machine$double.eps))
    }
  }
})
