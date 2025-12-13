#' Compute Finite Difference Operator (D)
#'
#' @description
#' Constructs the sparse discrete gradient operator matrix $D$,
#' to be acted on a vectorized image $u$ (length $n^2$).
#'
#' @param n Integer. The side length of the square image.
#'
#' @return A sparse matrix of class \code{dgCMatrix} with dimensions $2n^2 \times n^2$.
#'   The top half contains horizontal differences ($D_x$), and the bottom half
#'   contains vertical differences ($D_y$).
#'
#' @importFrom Matrix bandSparse kronecker Diagonal
#'
#' @keywords internal
#' @noRd
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

#' Compute Identity Operator (F)
#'
#' @description
#' Creates a sparse diagonal identity matrix.
#'
#' @param n Integer. The side length of the square image.
#'
#' @return A sparse diagonal matrix of dimensions $2n^2 \times 2n^2$.
#'
#' @keywords internal
#' @noRd
compute_F <- function(n) {
  size <- 2L * n * n
  return(Diagonal(size))
}

#' Isotropic Soft Thresholding (Shrinkage)
#'
#' @description
#' Applies generalized soft thresholding to the dual variables.
#'
#' @param q Numeric vector of length $2n^2$. Represents the stacked dual variables $[q_x; q_y]$.
#' @param tau Numeric scalar or vector. The threshold parameter.
#'
#' @return A numeric vector of the same length as \code{q}, containing the shrunk values.
#'
#' @keywords internal
#' @noRd
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
