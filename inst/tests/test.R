# Set up a greyscale palette
pal <- grey(seq(0, 1, length.out = 256))

# Im1: default grayscale
im1 <- generate_example_images()
par(mar = c(2,2,2,2))  # minimal margins
image(t(im1$corrupted[nrow(im1$corrupted):1,]),
      col = pal,
      useRaster = TRUE,
      main = "Corrupted Default Image")

# Im2: add missing and noise
im2 <- generate_example_images(missing_fraction = 0.1, noise_sigma = 0.2)
plot_images(list(im2$original, im2$corrupted), titles=c("Original", "Corrupt"))

# Im3: disk with blocks
im3 <- generate_example_images(pattern= "disk", mask_type = "random_blocks", missing_fraction = 0.2, noise_sigma = 0.2)
plot_images(list(im3$original, im3$corrupted), titles=c("Original", "Corrupt"))

