#' generate_example_images
#'
#' Generates an image with optional noise and missing‐pixel mask
#'
#' @param size Integer vector of length 2 giving number of rows and columns (default: c(64, 64)).
#' @param pattern Character. Type of base image to generate. Options include \code{"gradient"}, \code{"checkerboard"}, \code{"disk"} (default: \code{"gradient"}).
#' @param noise_sigma Numeric. Standard deviation of Gaussian noise to add to the image (default: 0.0, for no noise).
#' @param missing_fraction Numeric between 0 and 1. Fraction of pixels to mask (set to NA) randomly (default: 0.0, for no missing pixels).
#' @param mask_type Character. If \code{missing_fraction > 0}, type of mask to generate: \code{"random_pixels"} (default), \code{"random_blocks"}.
#' @param seed seed Integer (optional). Random seed for reproducibility.
#'
#' @returns A list with components:
#’ \item{original}{The original (clean) image matrix of size \code{size}.}
#’ \item{corrupted}{The image after adding noise and applying the mask (with NAs for missing pixels).}
#’ \item{mask}{Logical matrix of same size: \code{TRUE} means pixel is observed, \code{FALSE} means masked (missing).}
#'
#' @export
#'
#' @examples
#' ex <- generate_example_images(missing_fraction = 0.4, noise_sigma = 0.2)
#' image(ex$corrupted, main = "Example corrupted image", col=grey(seq(0,1,length.out = 256)))
generate_example_images <- function(size = c(64, 64), pattern = "gradient",
                                   noise_sigma = 0.0, missing_fraction = 0.0,
                                   mask_type = "disk", seed = NULL){
  # Set seed
  if (!is.null(seed)) set.seed(seed)
  nr <- size[1]
  nc <- size[2]
  # Generate base image
  if (pattern == "disk"){
    x <- seq(-1, 1, length.out = nc)
    y <- seq(-1, 1, length.out = nr)
    grid <- outer(y, x, FUN = function(yy, xx) sqrt(xx^2 + yy^2))
    original <- ifelse(grid <= 0.8, 1, 0)
  }
  else if (pattern == "gradient"){
    original <- outer(
      seq(0, 1, length.out = nr),
      seq(0, 1, length.out = nc),
      FUN = function(x, y)
        x
    )
  }
  else if (pattern == "checkerboard") {
    original <- outer(seq_len(nr), seq_len(nc),
                  FUN = function(i, j) (i + j) %% 2)
  }
  else {
    stop("Unknown pattern type: ", pattern)
  }

  # Add noise
  corrupted <- original + matrix(rnorm(nr * nc, sd = noise_sigma),
                             nrow = nr, ncol = nc)
  corrupted <- pmin(pmax(corrupted, 0), 1)

  # Generate Mask
  mask <- matrix(TRUE, nrow = nr, ncol = nc)
  if (missing_fraction > 0) {
    n_missing <- round(nr * nc * missing_fraction)
    if (mask_type == "random_pixels") {
      idx <- sample(nr * nc, size = n_missing, replace = FALSE)
      mask[idx] <- FALSE
    }
    else if (mask_type == "random_blocks") {
      # create random rectangular blocks
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
  list(original = original, corrupted = corrupted, mask = mask)
}




#' plot_images
#' Plot one or more grayscale images in a grid with titles
#'
#' @param imgs A list of image matrices (rows × cols) for grayscale images.
#' @param titles Character vector of same length as \code{imgs}. Titles to display above each image. If \code{NULL}, no individual titles.
#' @param cols Integer or \code{NULL}. Number of columns in the grid layout (by row-major order). If \code{NULL}, it uses \code{floor(sqrt(length(imgs)))} columns.
#' @param palette Character vector. Colour ramp for grayscale (default: \code{grey(seq(0,1,length.out=256))}).
#' @param useRaster Logical. Whether to use raster graphics (default: \code{TRUE}).
#' @param axes Logical. Whether to draw axes around each image (default: \code{FALSE}).
#' @param main Character scalar. Overall main title for the combined grid (default: \code{NULL}).
#’
#'
#' @export
#'
#' @examples
#' a <- matrix(rnorm(100), nrow=10)
#' b <- matrix(runif(160), nrow=10)
#' plot_images(list(a, b), titles=c("Norm", "Unif"))
plot_images <- function(imgs, titles = NULL, cols = NULL,
                        palette = grey(seq(0,1,length.out = 256)),
                        na_color = "red", useRaster = TRUE, axes = FALSE, main = NULL) {
  ### Check Input
  n <- length(imgs)
  if (n < 1) stop("No images to plot.")
  if (!is.null(titles) && length(titles) != n) {
    stop("Length of titles must match number of images.")
  }
  if (is.null(cols)) {
    cols <- floor(sqrt(n)); if (cols < 1) cols <- 1
  }
  rows <- ceiling(n / cols)
  oldpar <- par(mfrow = c(rows, cols), mar = c(0,0,2,0))
  on.exit(par(oldpar))

  if (!is.null(main)) {
    par(oma = c(0,0,2,0))
    title(main = main, outer = TRUE, line = 0)
  }

  ### Loop for every images in imgs
  for (i in seq_len(n)) {
    img <- imgs[[i]]
    if (!is.matrix(img)) stop("Each image must be a matrix for grayscale plotting.")

    # Create “special” value for NA
    special_val <- 1.1

    # Prepare modified image
    img2 <- img
    img2[is.na(img2)] <- 1.1

    # Extended palette
    new_palette <- c(palette, na_color)

    # Create breaks of length = length(new_palette) + 1
    breaks <- c(
      seq(0,1, length.out = length(new_palette)),
      1.1
    )

    image(
      t(img2[nrow(img2):1, ]),
      col = new_palette,
      breaks = breaks,
      useRaster = useRaster,
      axes = axes,
    )

    if (!is.null(titles)) {
      mtext(titles[i], side = 3, line = 0.5, cex = 1)
    }
  }

  invisible(NULL)
}


