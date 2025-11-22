#' Title
#'
#' @param img
#' @param mask
#' @param lambda
#' @param method
#' @param tau
#' @param sigma
#' @param theta
#'
#' @returns
#' @export
#'
#' @examplesa
tv_restore <- function(img, mask = NULL, lambda,
                       method = "primal_dual", task = "ROF_inpainting", u0 = NULL,
                       # primal-dual params
                       tau, sigma, theta,
                       max_iter = 300, tol = 1e-4, verbose = FALSE){
  ### Input Checks
  if (!is.matrix(img)) stop("`img` must be a 2D numeric matrix.")
  nr <- nrow(img); nc <- ncol(img)
  if (is.null(mask)) {
    mask <- !is.na(img)
  }
  if (!all(dim(mask) == c(nr, nc))) {
    stop("mask dimensions must match img.")
  }
  if (is.null(u0)){
    u0 <- as.numeric(ifelse(is.na(img), 0, img))
  }
  ###################MORE TO ADD LATER ###################

  # Match method string
  method <- match.arg(method, choices = c("primal_dual", "split_bregman"))
  task <- match.arg(task, choices = c("denoising", "inpainting", "ROF_inpainting"))

  # Call corresponding restore function
  if (method == "primal_dual") {
    # Vectorize
    im_vec <- as.numeric(ifelse(is.na(img), 0, img))   # length n^2
    mask_vec <- as.numeric(as.vector(mask))   # as.vector keeps column-major order

    # Solve saddle point
    res_vec <- find_saddle_point(im_vec = im_vec, mask_vec = mask_vec,
                      task = task,
                      lmda = lambda, u0 = u0, tol = tol, max_iter = max_iter,
                      verbose = verbose)
    # Return as matrix
    res <- matrix(res_vec, nrow = nr, ncol = nc)   # column-major fill
    return(res)
  }
}
