library(testthat)
library(ZILN)
context("Testing data cleaning features of ZILN")


test_check("ZILN")

test_that("The function outputs a list"{
  data <- read.csv('tests/rep1.csv')
  data.main <- data[,8:267]
  data.covar <- data[,4:7]
  expect_output(features.confirm(data.main,data.covar), )
})

features.confirm(data.main, data.covar)
