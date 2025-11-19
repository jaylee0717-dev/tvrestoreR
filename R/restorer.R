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

  # Initialize

  # Loop
}
