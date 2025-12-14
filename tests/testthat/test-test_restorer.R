test_that("tv_restore input validation works", {
  # Test non-matrix input
  expect_error(tv_restore(img = 1:10, lambda = 1), "must be a 2D numeric matrix")

  # Test non-square matrix
  rect_mat <- matrix(0, nrow = 10, ncol = 20)
  expect_error(tv_restore(img = rect_mat, lambda = 1), "must be a square matrix")

  # Test mismatched mask dimensions
  sq_mat <- matrix(0, 10, 10)
  bad_mask <- matrix(TRUE, 5, 5)
  expect_error(tv_restore(sq_mat, mask = bad_mask, lambda = 1, task = "inpainting"),
               "mask dimensions must match img")
})

test_that("tv_restore runs Primal-Dual Denoising correctly", {
  # Setup small synthetic data
  set.seed(123)
  img <- matrix(rep(0:1, each = 8), 4, 4) # 4x4 matrix
  img_noisy <- img + matrix(rnorm(16, 0, 0.1), 4, 4)

  # Run Denoising
  res <- tv_restore(img_noisy, lambda = 1, method = "primal_dual", task = "denoising")

  # Checks
  expect_true(is.matrix(res))
  expect_equal(dim(res), c(4, 4))
  expect_true(all(res >= 0 & res <= 1)) # Check clamping
  expect_false(any(is.na(res))) # Denoising should not return NAs
})

test_that("tv_restore auto-detects masks for Inpainting", {
  # Setup image with NAs
  img <- matrix(0.5, 10, 10)
  img[1:5, 1:5] <- NA

  # Should run without error and fill NAs
  expect_message(
    res <- tv_restore(img, lambda = 1, task = "inpainting", verbose = TRUE),
    "Auto-detected mask"
  )

  expect_false(any(is.na(res)))
})

test_that("tv_restore runs Split-Bregman solver", {
  img <- matrix(runif(100), 10, 10)
  # Run Split-Bregman
  res <- tv_restore(img, lambda = 1, method = "split_bregman", task = "denoising")

  expect_true(is.matrix(res))
  expect_equal(dim(res), c(10, 10))
  expect_true(all(res >= 0 & res <= 1))
})

test_that("tv_restore handles batch processing (List of Images)", {
  # Create list of 2 images
  img1 <- matrix(0, 8, 8)
  img2 <- matrix(1, 8, 8)
  img_list <- list(img1, img2)

  # Run batch
  res_list <- tv_restore(img_list, lambda = 1, task = "denoising", verbose = FALSE)

  expect_type(res_list, "list")
  expect_length(res_list, 2)
  expect_equal(dim(res_list[[1]]), c(8, 8))
  expect_equal(dim(res_list[[2]]), c(8, 8))
})


test_that("find_lambda validation works", {
  expect_error(find_lambda(img = 1:5, ground_truth = 1:5), "must be a matrix")

  # Fail if inpainting requested without mask/NAs
  clean_img <- matrix(0, 10, 10)
  expect_error(find_lambda(clean_img, ground_truth = clean_img, task = "inpainting"),
               "Mask required")
})

test_that("find_lambda executes grid search and returns valid structure", {
  # Setup synthetic problem
  set.seed(42)
  n <- 16
  truth <- matrix(0, n, n)
  truth[5:10, 5:10] <- 1

  # Add noise
  noisy <- truth + matrix(rnorm(n*n, 0, 0.2), n, n)

  # Define small grid
  l_grid <- c(0.1, 1, 10)

  # Run tuning (suppress messages for clean test output)
  res <- suppressMessages(find_lambda(
    img = noisy,
    ground_truth = truth,
    lambda_grid = l_grid,
    task = "denoising",
    method = "primal_dual",
    verbose = FALSE
  ))

  # Check Structure
  expect_type(res, "list")
  expect_named(res, c("best_lambda", "restored", "errors", "lambda_grid"))

  # Check Logic
  expect_true(res$best_lambda %in% l_grid)
  expect_true(is.matrix(res$restored))
  expect_equal(dim(res$restored), c(n, n))
  expect_length(res$errors, length(l_grid))

  # Check that errors are numeric and non-negative
  expect_true(all(res$errors >= 0))
})
