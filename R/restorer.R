tv_restore <- function(img, mask = NULL, lambda,
                       method = c("primal_dual"),
                       # primal-dual params
                       tau, sigma, theta){
  # Match method string
  method <- match.arg(method, choices = c("primal_dual", "split_bregman"))

  # Call corresponding restore function
  if (method == "primal_dual") {
    return(.tv_restore_primal_dual(img = img, mask = mask, lambda = lambda,
                                   tau = tau, sigma = sigma, theta = theta))
  }
}

tv_restore_primal_dual <- function(img = img, mask = mask, lambda = lambda,
                       tau = tau, sigma = sigma, theta = theta){
  # Input Checks

  #### Initialize
  # Replace NA in img with 0 (doesn't matter at masked locations)
  f <- img
  f[!mask] <- 0

  # Initialize primal (u) and dual (p = (p_x, p_y))
  u <- f  # initial guess
  u_bar <- u
  p_x <- matrix(0, nr, nc)
  p_y <- matrix(0, nr, nc)

  # Helper functions: gradient and divergence
  grad <- function(u) {
    # forward finite differences
    dx <- rbind(diff(u, 1, 1), rep(0, nc))
    dy <- cbind(diff(u, 1, 2), rep(0, nr))
    list(dx = dx, dy = dy)
  }
  div <- function(px, py) {
    # divergence is negative adjoint of gradient
    # backward differences
    ddx <- rbind(px[1, ], px[2:nr, ] - px[1:(nr-1), ])  # px_i,j - px_{i-1,j}
    ddy <- cbind(py[,1], py[,2:nc] - py[,1:(nc-1)])     # py_i,j - py_{i,j-1}
    ddx + ddy
  }

  # Proximal operator for fidelity term: G(u) = (λ/2) * mask * (u - f)^2
  prox_G <- function(u_tilde) {
    # closed form: (u_tilde + τ λ mask f) / (1 + τ λ mask)
    numerator <- u_tilde + tau * lambda * mask * f
    denom     <- 1 + tau * lambda * mask
    numerator / denom
  }

  # Dual step prox for F^*: F(p) = ||p||_2, so F^*(q) is indicator of ||q|| ≤ 1
  # Prox of sigma * F^* is projection of (p + sigma * grad(u_bar)) onto unit ball
  prox_Fs <- function(px_tilde, py_tilde) {
    norm_p <- sqrt(px_tilde^2 + py_tilde^2)
    # For each pixel, if norm > 1, scale down
    scale <- pmin(1, 1 / norm_p)
    list(
      px = px_tilde * scale,
      py = py_tilde * scale
    )
  }

  #### Loop
}
