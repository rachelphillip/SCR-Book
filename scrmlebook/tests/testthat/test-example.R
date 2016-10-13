## Tests for example chapter.
context("Example")
test_that("Example chapter tests",
{
    ## Load example fit.
    data("example-data")
    data("example-fits")
    ## These are the parameter estimates I expect.
    true.estimates <- c(9.2081067, 0.1606297, 27.6651518)
    ## These are the parameter estimates extracted from the model fit.
    estimates <- predict(example.fit)[, 2]
    ## We expect them to be equal (testing "equality" with a some tolerance for numerical error).
    expect_equal(estimates, expected = true.estimates, tolerance = 0.0001)
    ## We expect model fit to be an S3 class.
    expect_s3_class(example.fit, "secr")
})
