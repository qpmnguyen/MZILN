library(MZILN)
context("Testing MZILN cleaning features")

test_that("Function generates the correct messages when it needs to",{
  dat.nozero <- c(1,2,3,4,5)
  covar <- "testing"
  data("test.covariates")
  data("test.main")
  expect_error(features.confirm(test.main,dat.nozero),
               "Your data must be in data frame format and your covariates must be in vector format!")
  expect_error(features.confirm(test.main,test.covariates,covariates = covar),
               "Your specific covariates does not fit covariates listed in your covariate data frame!")
  test.main[,1] <- rep(0,nrow(test.main))
  test.covariates <- test.covariates[-c(1,2,3,4),]
  expect_warning(features.confirm(test.main, test.covariates),
                 "The number of subjects is inconsistent between data frames")
})



