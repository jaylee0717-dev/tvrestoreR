im1 <- generate_example_image()

# Set up a greyscale palette
pal <- grey(seq(0, 1, length.out = 256))

# Plot the corrupted image
par(mar = c(2,2,2,2))  # minimal margins
image(t(im1$corrupted[nrow(im1$corrupted):1,]),
      col = pal,
      useRaster = TRUE,
      main = "Corrupted Toy Image")

im2 <- generate_example_image(missing_fraction = 0.4, noise_sigma = 0.2)
image(t(im2$corrupted[nrow(im2$corrupted):1,]),
      col = pal,
      useRaster = TRUE,
      main = "Corrupted Toy Image")

plot_images(list(im1$corrupted, im2$corrupted), titles=c("Original", "Corrupt"))
