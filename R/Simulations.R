#' Generate Synthetic Images for Restoration
#'
#' @description
#' Creates synthetic test images with optional Gaussian noise and missing pixel masks.
#'
#' @param size Integer vector of length 2. Dimensions (rows, cols).
#' @param pattern Character String detailing what the base image should look.
#' Options: \code{"gradient"}, \code{"checkerboard"}, \code{"disk"}, \code{"wedge"}, or \code{"wedge60"}.
#' @param noise_sigma Numeric. Standard deviation of Gaussian noise to add.
#' @param missing_fraction Numeric (0 to 1). Fraction of pixels to mask as \code{NA}.
#' @param mask_type Character String, for strategy for masking. Options: \code{"random_pixels"},
#'   \code{"random_blocks"}, or \code{"disk"}.
#' @param seed Integer. Optional random seed for reproducibility.
#'
#' @export
#'
#' @returns A list with components:
#’ \item{original}{The original (clean) ground-truth image matrix of size \code{size}.}
#’ \item{corrupted}{The image after adding noise and applying the mask (with NAs for missing pixels).}
#’ \item{mask}{Logical matrix of same size: \code{TRUE} means pixel is observed, \code{FALSE} means masked (missing).}
#'
#' @examples
#' ex <- generate_example_images(missing_fraction = 0.4, noise_sigma = 0.2)
#' image(ex$corrupted, main = "Example corrupted image", col=grey(seq(0,1,length.out = 256)))
generate_example_images <- function(size = c(64, 64), pattern = "wedge",
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
  else if (pattern == "wedge") {
    ii <- matrix(rep(1:nr, each = nc), nrow = nr)
    jj <- matrix(rep(1:nc, nr), nrow = nr)

    original <- ifelse(ii > jj, 1, 0)
  }
  else if (pattern == "wedge60") {
    angle_deg <- 60
    theta <- angle_deg * pi/180

    x <- matrix(rep(1:nc, each = nr), nrow = nr) - (nc + 1)/2
    y <- matrix(rep(1:nr, nc), nrow = nr) - (nr + 1)/2

    line_val <- x * cos(theta) + y * sin(theta)
    original <- ifelse(line_val >= 0, 1, 0)
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

  # Return Original image, Corrupted image, and True/False mask
  list(original = original, corrupted = corrupted, mask = mask)
}




#' Grid Visualization of Image Lists
#'
#' @description
#' Plots a list of image matrices in a grid layout. Automatically handles \code{NA} values
#' by plotting them in a distinct color.
#'
#' @param imgs A list of image matrices (rows × cols) for grayscale images.
#' @param titles Character vector of same length as \code{imgs}. Titles to display above each image.
#' If \code{NULL}, no individual titles.
#' @param cols Integer or \code{NULL}. Number of columns in the grid layout (by row-major order).
#' If \code{NULL}, it uses \code{floor(sqrt(length(imgs)))} columns.
#' @param palette Character vector. Color ramp for grayscale.
#' @param useRaster Logical. Whether to use raster graphics
#' @param axes Logical. Whether to draw axes around each image
#' @param main Character scalar. Overall main title for the combined grid
#' @param na_color Character string. Sets which color "N/A" will be represented in.
#'
#' @importFrom grDevices grey
#' @importFrom graphics par title mtext image
#' @export
#'
#' @examples
#' a <- matrix(runif(100), 10, 10)
#' b <- matrix(runif(100), 10, 10)
#' b[1:5,1:5] <- NA
#' plot_images(list(a, b), titles = c("Clean", "Missing Data"), na_color = "red")
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
    img2[is.na(img2)] <- special_val

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


