#' tv_restore
#'
#' @param img
#' @param mask
#' @param lambda
#' @param method
#' @param task
#' @param u0 Initialized image. Defaults to an entirely black image.
#' @param max_iter Maximum iteration. Default 10000.
#' @param tol Tolerance on the 2-norm difference of primal and dual variables. Default 1e-5.
#' @param verbose Controls if detailed printouts on the looped variables will be provided. Default False.
#'
#' @returns
#' @export
#'
#' @examples
tv_restore <- function(img, mask = NULL, lambda,
                       method = "primal_dual", task = "ROF_inpainting", u0 = NULL,
                       max_iter = 10000, tol = 1e-4, verbose = FALSE, D= NULL,
                       # Task-specific parameters here
                       Fmat = NULL, ...){
  ###### IF img IS LIST OF IMAGES
  if (is.list(img)) {
    if (verbose) message(sprintf("Processing list of %d images...", length(img)))

    # Pre-compute D and F once (using size of first frame)
    n <- nrow(img[[1]])
    D_pre <- compute_D(n)

    # 1. Capture the dots and remove 'D' if it exists there to prevent collision
    args <- list(...)
    args$D <- NULL

    # Apply recursively
    res <- lapply(img, function(frame) {
      do.call(tv_restore, c(list(
        img = frame,
        mask = mask,
        lambda = lambda,
        task = task,
        method = method,
        u0 = u0,
        max_iter = max_iter,
        tol = tol,
        verbose = FALSE,
        D = D_pre
      ), args))
    })
    if (verbose) message("Completed.")
    return(res)
  }


  ###### IF IMG IS ONE IMAGE
  ## Input Checks
  # Match method string
  method <- match.arg(method, choices = c("primal_dual", "split_bregman"))
  task <- match.arg(task, choices = c("denoising", "inpainting", "ROF_inpainting"))


  if (!is.matrix(img)) stop("`img` must be a 2D numeric matrix.")
  if (nrow(img) != ncol(img)) stop("`img` must be a square matrix.")
  n <- nrow(img)

  # Check Mask
  if (task == "denoising") {
    if (!is.null(mask) && verbose) message("Note: Mask ignored for 'denoising' task.")
    mask <- matrix(1, nrow = nrow(img), ncol = ncol(img))
  } else {
    if (is.null(mask)) {
      if (any(is.na(img))) {
        mask <- !is.na(img)
        if (verbose) message("Auto-detected mask from NA values.")
      } else {
        stop(sprintf("Task '%s' requires a 'mask' argument.", task))
      }
    }
  }
  if (!all(dim(mask) == c(n, n))) {
    stop("mask dimensions must match img.")
  }

  # Initialize u0
  if (is.null(u0)){
    u0 <- as.numeric(ifelse(is.na(img), 0.5, img))
  }

  # Replace NAs in image with 0, vectorize inputs for efficiency
  im_vec <- as.numeric(ifelse(is.na(img), 0, img))   # length n^2
  mask_vec <- as.numeric(as.vector(mask))   # as.vector keeps column-major order


  ## Pre-Compute Large Operators if not provided
  args <- list(...)
  D_op <- if(is.null(args$D)) compute_D(n) else args$D
  F_op <- if(is.null(args$Fmat) && method == "primal_dual") compute_F(n) else args$Fmat


  #### Call corresponding restore function
  if (verbose) message(sprintf("Running %s using %s...", task, method))
  if (method == "primal_dual") {
    # Solve saddle point
    res_vec <- find_primaldual_saddle_point(im_vec = im_vec, mask_vec = mask_vec,
                      task = task, lmda = lambda, u0 = u0,
                      tol = tol, max_iter = max_iter, verbose = verbose,
                      D = D_op, Fmat = F_op, ...)
  }
  else if (method == "split_bregman") {
    res_vec <- find_split_bregman(im_vec, mask_vec, task, lambda, u0 = u0,
                                  tol = tol, max_iter = max_iter, verbose = verbose,
                                  D = D_op, ...)
  }
  # Return as matrix, strictly capped between 0 and 1
  res <- matrix(res_vec, nrow = n, ncol = n)
  res[res < 0] <- 0; res[res > 1] <- 1;
  if (task == "denoising")
    res[!mask] <- NA
  return(res)
}


#' Tune Regularization Parameter (Lambda)
#'
#' Performs a grid search to find the optimal lambda.
#'
#' @inheritParams tv_restore
#' @param lambda_grid Vector of lambda values to test.
#' @param ground_truth Ground truth matrix for calculating true MSE.
#' @export
find_lambda <- function(img, mask = NULL, ground_truth,
                        lambda_grid = 10^seq(-2, 4, length.out = 12),
                        task = c("denoising", "inpainting", "ROF_inpainting"),
                        method = c("primal_dual", "split_bregman"),
                        verbose = TRUE,
                        ...) {

  ##  Validation
  task <- match.arg(task)
  method <- match.arg(method)
  if (!is.matrix(img)) stop("'img' must be a matrix.")

  # Auto-detect mask if needed
  if (task != "denoising" && is.null(mask)) {
    if (any(is.na(img))) mask <- !is.na(img)
    else stop("Mask required for inpainting tasks.")
  }
  if (is.null(mask)) mask <- matrix(TRUE, nrow=nrow(img), ncol=ncol(img))


  # --- 2. Pre-computation ---
  n <- nrow(img)
  if (verbose) message("Pre-computing operators...")
  D_pre <- compute_D(n)
  F_pre <- if (method == "primal_dual") compute_F(n) else NULL

  # Setup storage
  errors <- numeric(length(lambda_grid))
  names(errors) <- lambda_grid

  best_res <- NULL
  min_err <- Inf
  best_lam <- NA

  # Setup Error Metric Data
  observed_idx <- which(mask)
  obs_data <- img[observed_idx]
  # Replace NA in input for solver safety
  img[is.na(img)] <- 0

  # Grid Search Loop
  if (verbose) message(sprintf("Starting grid search for '%s' using '%s'...", task, method))

  for (i in seq_along(lambda_grid)) {
    lam <- lambda_grid[i]
    if (verbose) message(sprintf("  [%d/%d] Testing lambda = %.4f...", i, length(lambda_grid), lam))

    # Run tv_restore
    res_mat <- tryCatch({
      tv_restore(img = img, mask = mask, lambda = lam,
                 task = task, method = method,
                 verbose = FALSE,
                 D = D_pre, Fmat = F_pre,
                 ...)
    }, error = function(e) {
      warning(sprintf("Solver failed for lambda=%.4f: %s", lam, e$message))
      return(NULL)
    })

    if (is.null(res_mat)) { errors[i] <- Inf; next }

    # Calculate Error

    err <- mean((res_mat - ground_truth)^2, na.rm = TRUE)



    errors[i] <- err
    if (verbose) message(sprintf("    Error: %.6e", err))
    if (err < min_err) {
      min_err <- err
      best_lam <- lam
      best_res <- res_mat
    }
    rm(res_mat)
    gc(verbose = FALSE)
  }


  if (verbose) message(sprintf("Search complete. Best lambda: %.4f", best_lam))
  return(list(
    best_lambda = best_lam,
    restored = best_res,
    errors = errors,
    lambda_grid = lambda_grid
  ))
}
