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


