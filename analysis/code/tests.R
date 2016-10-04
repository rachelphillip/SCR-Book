## Tests for code in analysis/code/. Run all tests using test_file("tests.R")

## Loading all packages (and install any that are not yet installed).
source("packages.R")

## Tests for example chapter.
context("Example")
test_that("Example chapter tests",
{
    ## Load example fit.
    load("../objects/example.RData")
    ## These are the parameter estimates I expect.
    true.estimates <- c(9.2081067, 0.1606297, 27.6651518)
    ## These are the parameter estimates extracted from the model fit.
    estimates <- predict(example.fit)[, 2]
    ## We expect them to be equal (testing "equality" with a some tolerance for numerical error).
    expect_equal(estimates, expected = true.estimates, tolerance = 0.0001)
    ## We expect model fit to be an S3 class.
    expect_s3_class(example.fit, "secr")
})

## Tests for detector types chapter.
context("Detector types")
test_that("Detector types chapter test",
{
    ## Tests in here.
})
