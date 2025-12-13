#' Proximal Operator for the Fidelity Term
#'
#' @description
#' Computes the proximal operator for the data fidelity term.
#'
#' @param ubar Numeric vector. The variable to update.
#' @param g Numeric vector. The observed input image.
#' @param tau Numeric scalar. Step size parameter. Required for denoising/ROF tasks.
#' @param mask Numeric vector. The observation mask (0/1).
#' @param task_type Character. One of \code{"inpainting"}, \code{"denoising"}, or \code{"ROF_inpainting"}.
#'
#' @return A numeric vector of the same length as \code{ubar}.
#' @keywords internal
#' @noRd
prox <- function(ubar, g, tau = NULL, mask, task_type = c("inpainting", "denoising", "ROF_inpainting")) {
  task_type <- match.arg(task_type)
  if (task_type == "inpainting") {
    return(g * mask + ubar * (1 - mask))
  } else if (task_type == "denoising") {
    if (is.null(tau)) stop("tau must be provided for denoising")
    return((ubar + tau * g) / (1 + tau))
  } else if (task_type == "ROF_inpainting") {
    if (is.null(tau)) stop("tau must be provided for ROF_inpainting")
    uhat <- (ubar + tau * g) / (1 + tau)
    return(uhat * mask + ubar * (1 - mask))
  } else {
    stop("Need Valid Imaging Task Type: Either inpainting or denoising.")
  }
}


#' Primal-Dual Saddle Point Solver
#'
#' @description
#' Implements the Chambolle-Pock first-order primal-dual algorithm to solve the
#' Total Variation saddle point problem.
#'
#' @param im_vec Numeric vector of length $n^2$. Vectorized input image.
#' @param mask_vec Numeric vector of length $n^2$. Vectorized mask (0/1).
#' @param task Character. Task type: \code{"denoising"}, \code{"inpainting"}, or \code{"ROF_inpainting"}.
#' @param lmda Numeric. Regularization parameter $\lambda$.
#' @param u0 Numeric vector. Initial guess for the solution.
#' @param tol Numeric. Convergence tolerance for the relative change in primal and dual variables.
#' @param max_iter Integer. Maximum number of iterations.
#' @param verbose Logical. If \code{TRUE}, prints convergence metrics every 1000 iterations.
#' @param D Sparse matrix. The discrete gradient operator (stacked $D_x, D_y$).
#' @param Fmat Sparse matrix. Operator linking the auxiliary variable (currently identity).
#'
#' @return A numeric vector of length $n^2$ representing the restored image.
#'
#' @import Matrix
#' @keywords internal
#' @noRd
find_primaldual_saddle_point <- function(im_vec, mask_vec, task,
                              lmda, u0, tol,
                              D, Fmat,
                              max_iter, verbose = FALSE) {

  # --- Helper for Norm ---
  vecnorm2 <- function(x) sqrt(sum(x^2))

  # --- Input Checks ---
  n2 <- length(im_vec)
  n <- as.integer(sqrt(n2))

  if (n * n != n2) stop("im_vec length must be a perfect square")
  if (length(mask_vec) != n2) stop("mask_vec length mismatch")
  if (length(u0) != n2) stop("u0 length mismatch")

  # --- Setup Matrices ---
  im_copy <- as.numeric(im_vec)

  # --- Initialization ---
  u <- as.numeric(u0)
  q <- rep(1.0, 2 * n2)
  p <- rep(1.0, 2 * n2)

  # --- Stepsize Calculation ---
  row_sum_abs_D <- rowSums(abs(D))
  sigma_p <- 1 / max(row_sum_abs_D + 1)

  col_sum_abs_D <- colSums(abs(D))
  tau_u <- 1 / max(col_sum_abs_D)

  tau_q <- 1.0

  #### --- Main Loop --- ####
  count <- 0L
  rel_delta_u <- Inf
  rel_delta_p <- Inf
  eps <- 1e-10

  while (TRUE) {
    # Primal Update (u)
    Dt_p <- as.numeric(crossprod(D, p))

    unew <- prox(u - tau_u * Dt_p,
                 im_copy,
                 tau = tau_u,
                 mask = mask_vec,
                 task_type = task)

    ubar <- 2 * unew - u

    # q Update
    # Fp <- as.numeric(Fmat %*% p)
    Fp <- p

    q_in <- q + tau_q * Fp
    qnew <- shrink(q_in, tau = tau_q * lmda)
    qbar <- 2 * qnew - q

    # p Update
    D_ubar <- as.numeric(D %*% ubar)
    # t(F) %*% qbar. Since F is identity, this is just qbar.
    # Ft_qbar <- as.numeric(crossprod(Fmat, qbar))
    Ft_qbar <- qbar

    pnew <- p + sigma_p * (D_ubar - Ft_qbar)

    ## Convergence Checks
    if (count > max_iter) break

    rel_delta_u <- vecnorm2(unew - u) / (vecnorm2(u) + eps)
    rel_delta_p <- vecnorm2(pnew - p) / (vecnorm2(p) + eps)
    if (rel_delta_u < tol && rel_delta_p < tol) break

    # Commit Updates
    u <- unew
    q <- qnew
    p <- pnew
    count <- count + 1L

    ## Output Control
    if (verbose) {
      if ((count %% 1000) == 0L) {
        cat(sprintf("Iteration %d | dU: %.2e | dP: %.2e\n",
                    count, rel_delta_u, rel_delta_p))
      }
    }
  }

  if (verbose) cat("Total iterations:", count, "\n")
  return(u)
}
