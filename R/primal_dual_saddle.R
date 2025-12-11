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

shrink <- function(q, tau) {
  if (!is.numeric(q)) stop("q must be numeric")
  len_q <- length(q)
  if (len_q %% 2 != 0) stop("length(q) must be even (2 * n^2)")

  # Define split point for Stacked Layout
  n_pixels <- len_q / 2

  # Slice first half (x) and second half (y)
  qx <- q[1:n_pixels]
  qy <- q[(n_pixels + 1):len_q]

  # Compute two-norm per pixel
  two_norms <- sqrt(qx^2 + qy^2)

  # Tau handling (scalar vs vector)
  if (length(tau) == 1) {
    tau_eff <- tau
  } else if (length(tau) == n_pixels) {
    tau_eff <- tau
  } else {
    stop("tau must be scalar or length n^2")
  }

  # Scaling factor
  # Equivalent to: max(two_norms/tau, 1)
  denom <- pmax(two_norms / tau_eff, 1)
  scaling <- 1 - 1 / denom

  # Apply scaling
  qx_shr <- qx * scaling
  qy_shr <- qy * scaling

  # Concatenate back to Stacked Layout [qx_shr; qy_shr]
  return(c(qx_shr, qy_shr))
}

compute_D <- function(n) {
  if (n <= 1) stop("n must be > 1")
  # small n x n forward-difference matrix: diag0=-1, diag+1=1, last +1 entry zero
  diag0 <- rep(-1, n)
  diag1 <- c(rep(1, n - 1), 0)
  small <- bandSparse(n, k = c(0, 1), diagonals = list(diag0, diag1))
  Dx <- kronecker(Diagonal(n), small, make.dimnames = FALSE)
  Dy <- kronecker(small, Diagonal(n), make.dimnames = FALSE)
  D <- rbind(Dx, Dy)
  return(D)
}


compute_F <- function(n) {
  size <- 2L * n * n
  return(Diagonal(size))
}


#' find_primaldual_saddle_point
#'
#' @param im_vec
#' @param mask_vec
#' @param task
#' @param lmda
#' @param u0
#' @param tol
#' @param max_iter
#' @param verbose
#'
#' @returns
#' @export
#' @import Matrix
#'
#' @examples
find_primaldual_saddle_point <- function(im_vec, mask_vec, task,
                              lmda, u0, tol, max_iter, verbose = FALSE) {

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
  D <- compute_D(n) # returns stacked layout [Dx; Dy]
  Fmat <- compute_F(n)

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
