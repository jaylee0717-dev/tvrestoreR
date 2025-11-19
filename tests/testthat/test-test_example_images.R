library(testthat)

test_that("generate_example_images returns correct structure and sizes", {
  ex <- generate_example_images(
    size = c(50, 80),
    pattern = "checkerboard",
    noise_sigma = 0.05,
    missing_fraction = 0.1,
    mask_type = "random_pixels",
    seed = 123
  )

  # Expect it returns a list with named components
  expect_true(is.list(ex))
  expect_named(ex, c("original", "corrupted", "mask"))

  # Check sizes
  expect_equal(dim(ex$original), c(50, 80))
  expect_equal(dim(ex$corrupted), c(50, 80))
  expect_equal(dim(ex$mask), c(50, 80))

  # Check that mask is logical
  expect_true(is.logical(ex$mask))

  # Check that corrupted has NAs exactly where mask is FALSE
  expect_equal(which(!ex$mask, arr.ind = TRUE),
               which(is.na(ex$corrupted), arr.ind = TRUE))

  # Check that original and corrupted differ when noise_sigma > 0
  expect_false(all.equal(ex$original, ex$corrupted, check.attributes = FALSE))

  # Check missing fraction approximately correct
  missing_frac_observed <- sum(!ex$mask) / (50 * 80)
  expect_lt(abs(missing_frac_observed - 0.10), 0.05)
})

test_that("generate_example_images errors on invalid pattern or mask_type", {
  expect_error(generate_example_images(pattern = "unknown_pattern"),
               "Unknown pattern type")
  expect_error(generate_example_images(missing_fraction = 0.2, mask_type = "unknown_mask"),
               "Unknown mask_type")
})
