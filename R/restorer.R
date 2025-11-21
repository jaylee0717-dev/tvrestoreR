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
#' @examples
tv_restore <- function(img, mask = NULL, lambda,
                       method = "primal_dual",
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
  ###################MORE TO ADD LATER ###################

  # Match method string
  method <- match.arg(method, choices = c("primal_dual", "split_bregman"))

  # Call corresponding restore function
  if (method == "primal_dual") {
    return(tv_restore_primal_dual(img = img, mask = mask, lambda = lambda,
                                   tau = tau, sigma = sigma, theta = theta,
                                  max_iter = max_iter, tol = tol, verbose = verbose))
  }
}



#' Title
#'
#' @param img
#' @param mask
#' @param lambda
#' @param tau
#' @param sigma
#' @param theta
#'
#' @returns
#' @export
#'
#' @examples
tv_restore_primal_dual <- function(img, mask, lambda, tau, sigma, theta, max_iter, tol, verbose){
  # Assume all inputs are fine at this point

  #### Initialize
  nr <- nrow(img); nc <- ncol(img)
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
    nr <- nrow(u); nc <- ncol(u)
    # forward finite differences
    dx <- matrix(0, nr, nc)
    dy <- matrix(0, nr, nc)

    # difference in row direction (vertical): u[i+1, j] - u[i, j]
    dx[1:(nr - 1), ] <- u[2:nr, ] - u[1:(nr - 1), ]

    # difference in column direction (horizontal): u[i, j+1] - u[i, j]
    dy[, 1:(nc - 1)] <- u[, 2:nc] - u[, 1:(nc - 1)]

    list(dx = dx, dy = dy)
  }

  div <- function(px, py) {
    nr <- nrow(px); nc <- ncol(px)
    # divergence: backward difference
    ddx <- matrix(0, nr, nc)
    ddy <- matrix(0, nr, nc)

    # px difference: p_x[i,j] - p_x[i-1,j] (for i > 1)
    ddx[2:nr, ] <- px[2:nr, ] - px[1:(nr - 1), ]
    # for first row, ddx stays zero (or boundary condition)

    # py difference: p_y[i,j] - p_y[i,j-1] (for j > 1)
    ddy[, 2:nc] <- py[, 2:nc] - py[, 1:(nc - 1)]

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
  for (k in seq_len(max_iter)) {
    u_old <- u

    # dual
    g <- grad(u_bar)
    px_t <- p_x + sigma * g$dx
    py_t <- p_y + sigma * g$dy
    p_proj <- prox_Fs(px_t, py_t)
    p_x <- p_proj$px; p_y <- p_proj$py

    # primal
    div_p <- div(p_x, p_y)
    u_t <- u - tau * div_p
    u <- prox_G(u_t)

    # extrapolate
    u_bar <- u + theta * (u - u_old)

    # stopping
    rel <- sqrt(sum((u - u_old)^2)) / (sqrt(sum(u_old^2)) + 1e-12)
    if (verbose && (k %% 50 == 0)) {
      message("iter ", k, ": rel = ", signif(rel, 3))
    }
    if (rel < tol) {
      if (verbose) message("Converged at iteration ", k)
      break
    }
  }

  g_f <- grad(u)
  tv_val <- sum(sqrt(g_f$dx^2 + g_f$dy^2))
  fidelity <- 0.5 * lambda * sum((u[mask] - f[mask])^2)
  obj <- fidelity + tv_val

  # Return List of restored, iterations, and objective
  list(restored = u, iterations = k, objective = obj)
}
