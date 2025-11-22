library(Matrix)

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
    return(g * mask + uhat * (1 - mask))
  } else {
    stop("Need Valid Imaging Task Type: Either inpainting or denoising.")
  }
}

shrink <- function(q, tau) {
  if (!is.numeric(q)) stop("q must be numeric")
  if (length(q) %% 2 != 0) stop("length(q) must be even (2 * n^2)")
  n2 <- length(q) / 2
  # split into two components
  qx <- q[seq(1, length(q), by = 2)]
  qy <- q[seq(2, length(q), by = 2)]
  # compute two-norm per pixel
  two_norms <- sqrt(qx^2 + qy^2)
  # tau handling
  if (length(tau) == 1) {
    tau_vec <- rep(tau, n2)
  } else if (length(tau) == n2) {
    tau_vec <- tau
  } else {
    stop("tau must be scalar or length n^2")
  }
  denom <- pmax(two_norms / tau_vec, 1)
  scaling <- 1 - 1 / denom
  qx_shr <- qx * scaling
  qy_shr <- qy * scaling
  # interleave back to [qx1, qy1, qx2, qy2, ...] to match input ordering
  out <- numeric(2 * n2)
  out[seq(1, length(out), by = 2)] <- qx_shr
  out[seq(2, length(out), by = 2)] <- qy_shr
  return(out)
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


# find_saddle_point: R translation of Python algorithm (L=1, FD only)
# im_vec, mask_vec: numeric vectors of length n^2
# dctzn_scheme: should be "FD"
# dctzn_values: ignored (kept for API compatibility)
# task: "inpainting" or "denoising" or "ROF_inpainting"
# lmda: regularization parameter (lambda)
# u0: initial u (numeric vector length n^2)
# stopping_threshold, maximum_iterations: numeric / integer
# update_interval: integer for console prints
find_saddle_point <- function(im_vec, mask_vec,
                              task = c("inpainting", "denoising", "ROF_inpainting"),
                              lmda = 1.0,
                              u0,
                              stopping_threshold = 1e-6,
                              maximum_iterations = 2000,
                              update_interval = 50) {
  vecnorm2 <- function(x) sqrt(sum((x)^2))

  task <- match.arg(task)
  n2 <- length(im_vec)
  n <- as.integer(sqrt(n2))
  if (n * n != n2) stop("im_vec length must be a perfect square")
  if (length(mask_vec) != n2) stop("mask_vec length mismatch")
  if (length(u0) != n2) stop("u0 length mismatch")

  # copies
  im_copy <- as.numeric(im_vec)
  D <- compute_D(n)           # (2*n^2) x (n^2) sparse
  F <- compute_F(n)           # (2*n^2) x (2*n^2) identity for FD

  # L = 1 (assumption)
  L <- 1L

  # initialize variables
  u <- as.numeric(u0)
  q <- rep(1.0, 2 * n2)       # q length = 2 * n^2, interleaved [qx1,qy1,qx2,qy2,...]
  p <- rep(1.0, 2 * n2)       # p length = 2 * n^2 (dual for data term)

  # constant stepsizes (using matrix absolute sums)
  # Build C = [D, F^T] horizontally
  # Note: F is identity so F^T is identity; but keep general form
  C <- cbind(D, t(F))
  # row sums of abs(C)
  row_sum_abs_C <- rowSums(abs(C))
  sigma_p <- 1 / max(row_sum_abs_C)
  # tau for u: 1 / max column sums of abs(D)
  col_sum_abs_D <- colSums(abs(D))
  tau_u <- 1 / max(col_sum_abs_D)
  # tau for q: 1 / max column sums of abs(F^T)  (F^T columns)
  col_sum_abs_Ft <- colSums(abs(t(F)))
  tau_q <- 1 / max(col_sum_abs_Ft)

  # main loop
  count <- 1L
  rel_delta_u <- Inf
  rel_delta_p <- Inf
  eps <- 1e-10

  while (TRUE) {
    # primal update (u)
    # D^T %*% p : yields length n^2 vector
    Dt_p <- as.numeric(t(D) %*% p)
    unew <- prox(u - tau_u * Dt_p, im_copy, tau = tau_u, mask = mask_vec, task_type = task)
    ubar <- 2 * unew - u

    # dual q update (shrink) : qnew = shrink(q + tau_q * F %*% p, tau_q * lmda)
    Fp <- as.numeric(F %*% p)
    q_in <- q + tau_q * Fp
    qnew <- shrink(q_in, tau = tau_q * lmda)
    qbar <- 2 * qnew - q

    # dual p update
    D_ubar <- as.numeric(D %*% ubar)
    Ft_qbar <- as.numeric(t(F) %*% qbar)
    pnew <- p + sigma_p * (D_ubar - Ft_qbar)

    # relative changes (2-norm)
    rel_delta_u <- vecnorm2(unew - u) / (vecnorm2(u) + eps)
    rel_delta_p <- vecnorm2(pnew - p) / (vecnorm2(p) + eps)

    # commit updates
    u <- unew
    q <- qnew
    p <- pnew

    # output control
    if ((count %% update_interval) == 0L) {
      cat("Iteration", count, "\n")
      cat(sprintf("Relative Norm change on u: %.6e\n", rel_delta_u))
      cat(sprintf("Relative Norm change on p: %.6e\n", rel_delta_p))
      cat(rep("-", 30), "\n")
    }
    count <- count + 1L

    # stopping criteria
    if ((count > maximum_iterations) ||
        (rel_delta_u < stopping_threshold && rel_delta_p < stopping_threshold)) {
      break
    }
  }

  cat("total iteration:", count, "\n")
  return(u)
}
