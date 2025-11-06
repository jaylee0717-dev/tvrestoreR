generate_example_image <- function(size = c(64, 64), pattern = "gradient",
                                   noise_sigma = 0.0, missing_fraction = 0.0,
                                   mask_type = "disk", seed = NULL){
  # Set seed
  if (!is.null(seed)) set.seed(seed)
  nr <- size[1]
  nc <- size[2]
  # Generate base image
  if (pattern == "disk"){
    x <- seq(-1,1,length.out=nc)
    y <- seq(-1,1,length.out=nr)
    grid <- outer(y, x, FUN = function(yy, xx) sqrt(xx^2 + yy^2))
    orig <- ifelse(grid <= 0.8, 1, 0)
  }
  else if (pattern == "gradient"){
    orig <- outer(seq(0,1,length.out=nr), seq(0,1,length.out=nc), FUN = function(x, y) x)
  }
  else if (pattern == "checkerboard") {
    orig <- matrix(rep(rep(c(0,1), length.out = nc), length.out = nr),
                   nrow = nr, byrow = TRUE)
  # Add noise
  corrupted <- orig + noise_sigma * matrix(rnorm(nr * nc), nrow = nr, ncol = nc)
  # Generate Mask

  # Return Original image, Corrupted image, and mask pixels
}
