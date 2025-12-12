#' compute_D
#'
#' @param n
#'
#' @returns
#'
#' @importFrom Matrix bandSparse kronecker Diagonal
#' @noRd
#'
#'
#' @examples
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


#' Title
#'
#' @param q
#' @param tau
#'
#' @returns
#' @noRd
#'
#' @examples
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
