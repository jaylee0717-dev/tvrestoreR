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
  }
  else {
    stop("Unknown pattern type: ", pattern)
  }

  # Add noise
  corrupted <- orig + noise_sigma * matrix(rnorm(nr * nc), nrow = nr, ncol = nc)
  # Generate Mask
  mask <- matrix(TRUE, nrow = nr, ncol = nc)
  if (missing_fraction > 0) {
    n_missing <- round(nr * nc * missing_fraction)
    if (mask_type == "random_pixels") {
      idx <- sample(nr * nc, size = n_missing, replace = FALSE)
      mask[idx] <- FALSE
    }
    else if (mask_type == "random_blocks") {
      # create a random rectangular block
      block_rows <- sample(1:nr, size = round(sqrt(n_missing)), replace = FALSE)
      block_cols <- sample(1:nc, size = round(sqrt(n_missing)), replace = FALSE)
      mask[block_rows, block_cols] <- FALSE
    }
    else if (mask_type == "disk"){
      x <- seq(-1,1,length.out=nc)
      y <- seq(-1,1,length.out=nr)
      grid <- outer(y, x, FUN = function(yy, xx) sqrt(xx^2 + yy^2))
      r <- sqrt(missing_fraction)
      mask[grid <= r] <- FALSE
    }
    else {
      stop("Unknown mask_type: ", mask_type)
    }
  }
  corrupted[!mask] <- NA

  # Return Original image, Corrupted image, and mask pixel
  list(orig = orig, corrupted = corrupted, mask = mask)
}



plot_images <- function(imgs, titles = NULL, cols = NULL,
                        palette = grey(seq(0,1,length.out = 256)),
                        useRaster = TRUE, axes = FALSE, main = NULL) {
  n <- length(imgs)
  if (n < 1) stop("No images to plot.")
  if (!is.null(titles) && length(titles) != n) {
    stop("Length of titles must match number of images.")
  }
  # decide layout
  if (is.null(cols)) {
    cols <- floor(sqrt(n))
    if (cols < 1) cols <- 1
  }
  rows <- ceiling(n / cols)
  oldpar <- par(mfrow = c(rows, cols), mar = c(0,0,2,0))
  on.exit(par(oldpar))

  if (!is.null(main)) {
    # set outer title
    par(oma = c(0,0,2,0))
    title(main = main, outer = TRUE, line = 0)
  }

  for (i in seq_len(n)) {
    img <- imgs[[i]]
    if (!is.matrix(img)) {
      stop("Each image must be a matrix for grayscale plotting.")
    }
    # Flip & transpose for correct orientation
    image(t(img[nrow(img):1, ]),
          col = palette,
          useRaster = useRaster,
          axes = axes)
    if (!is.null(titles)) {
      mtext(titles[i], side = 3, line = 0.5, cex = 1)
    }
  }
  invisible(NULL)
}
