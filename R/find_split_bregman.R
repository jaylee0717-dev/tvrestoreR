#' Title
#'
#' @param im_vec
#' @param mask_vec
#' @param task
#' @param lmda
#' @param u0
#' @param mu
#' @param tol
#' @param max_iter
#' @param verbose
#' @param D
#'
#' @returns
#' @export
#'
#' @examples
find_split_bregman <- function(im_vec, mask_vec, task,
                               lmda, u0, mu = 1.0,
                               tol = 1e-4, max_iter = 500, verbose = FALSE,
                               D = NULL) {

  vecnorm2 <- function(x) sqrt(sum(x^2))
  n2 <- length(im_vec)
  n <- as.integer(sqrt(n2))

  # Allocate D
  if (is.null(D)) {
    D <- compute_D(n)
  }

  # Compute Laplacian (D^T D) once
  Laplacian <- crossprod(D)

  # Build Linear System according to task (with stability padding)
  if (task == "inpainting") {
    LHS <- mu * Laplacian + Diagonal(n=length(im_vec), x=1e-10)
  } else { # Soft constraints (Denoising / ROF_inpainting)
    LHS <- Diagonal(x = mask_vec) + mu * Laplacian
  }
  solver_chol <- Cholesky(LHS)

  ## Initialization
  u <- as.numeric(u0)
  d <- numeric(nrow(D))
  b <- numeric(nrow(D))


  count <- 0L
  ## Main Loop
  while (TRUE) {
    u_old <- u

    # update u
    rhs_data <- mask_vec * im_vec
    rhs_reg  <- mu * as.numeric(crossprod(D, d - b))
    RHS <- rhs_data + rhs_reg
    # Solve Au = b
    u <- as.numeric(solve(solver_chol, RHS))
    if (task == "inpainting") {
      # Reset known pixels to original values
      u[mask_vec == 1] <- im_vec[mask_vec == 1]
    }

    # update d
    # d = shrink(Du + b, lambda / mu)
    Du <- as.numeric(D %*% u)
    d <- shrink(Du + b, tau = (lmda / mu))

    # update b
    b <- b + (Du - d)

    count <- count + 1L

    ## Check Convergence
    if (count > max_iter) break

    rel_change <- vecnorm2(u - u_old) / (vecnorm2(u) + 1e-10)
    if (rel_change < tol) break


    ## Output Control
    if (verbose && count %% 50 == 0) {
      cat(sprintf("Iter %d | Rel Change: %.2e\n", count, rel_change))
    }
  }

  return(u)
}
