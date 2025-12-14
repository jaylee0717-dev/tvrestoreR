#' Split Bregman Solver for TV Restoration
#'
#' @description
#' Solves the Total Variation problem using the Split Bregman iteration technique.
#'
#' @param im_vec Numeric vector. Vectorized input image.
#' @param mask_vec Numeric vector. Vectorized mask (0/1).
#' @param task Character. Task type. Note: \code{"inpainting"} enforces hard constraints via high penalty weights.
#' @param lmda Numeric. Regularization parameter $\lambda$.
#' @param u0 Numeric vector. Initial guess.
#' @param mu Numeric. Penalty parameter for the split variable.
#' @param tol Numeric. Convergence tolerance.
#' @param max_iter Integer. Maximum iterations.
#' @param verbose Logical. Print progress if TRUE.
#' @param D Sparse matrix. The discrete gradient operator.
#'
#'
#' @return A numeric vector representing the restored image.
#'
#' @importFrom Matrix Diagonal Cholesky crossprod solve
#' @keywords internal
#' @noRd
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

  # --- 1. DETERMINE WEIGHTS & THRESHOLDS ---
  if (task == "inpainting") {
    fidelity_weight <- 10000.0 # Huge penalty to anchor known pixels
    shrink_tau <- 1.0 / mu

  } else {
    fidelity_weight <- 1.0
    shrink_tau <- 1.0 / (lmda * mu)
  }

  # Compute Laplacian (D^T D) once
  Laplacian <- crossprod(D)

  # Build LHS Matrix
  Fidelity <- Diagonal(x = (mask_vec * fidelity_weight))
  LHS <- Fidelity + (mu * Laplacian) + Diagonal(n = n2, x = 1e-10)

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
    rhs_data <- (mask_vec * fidelity_weight * im_vec)
    rhs_reg  <- mu * as.numeric(crossprod(D, d - b))
    # Solve Au = b
    u <- as.numeric(solve(solver_chol, rhs_data + rhs_reg))
    if (task == "inpainting") {
      # Reset known pixels to original values
      u[mask_vec == 1] <- im_vec[mask_vec == 1]
    }

    # update d
    Du <- as.numeric(D %*% u)
    d <- shrink(Du + b, tau = shrink_tau)

    # update b
    b <- b + (Du - d)

    count <- count + 1L

    ## Check Convergence
    if (count > max_iter) break

    rel_change <- vecnorm2(u - u_old) / (vecnorm2(u) + 1e-10)
    if (rel_change < tol && count > 5) break


    ## Output Control
    if (verbose && count %% 50 == 0) {
      cat(sprintf("Iter %d | Rel Change: %.2e\n", count, rel_change))
    }
  }

  if (verbose) cat("Total iterations:", count, "\n")
  return(u)
}
