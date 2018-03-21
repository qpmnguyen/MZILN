library(MZILN)
context("Testing MZILN's internal features")

test_that("Testing that the output of each function what expected",{
  data('test.main')
  data("test.covariates")
  nonzero <- which(test.main[1,] != 0)
  expect_output(log.trans.generation())
})
