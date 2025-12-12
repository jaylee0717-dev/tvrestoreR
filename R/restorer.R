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
                       max_iter = 10000, tol = 1e-4, verbose = FALSE,
                       # Task-specific parameters here
                       ...){
  ### Input Checks
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
