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
im2 <- generate_example_images(size= c(128, 128), missing_fraction = 0., noise_sigma = 0.05)
plot_images(list(im2$original, im2$corrupted), titles=c("Original", "Corrupt"))

# Im3: disk with blocks
im3 <- generate_example_images(pattern= "disk", mask_type = "random_blocks", missing_fraction = 0.1, noise_sigma = 0.1)
plot_images(list(im3$original, im3$corrupted), titles=c("Original", "Corrupt"))

# Im4: checkerboard with random pixels
im4 <- generate_example_images(pattern= "checkerboard", mask_type = "random_pixels", missing_fraction = 0.2, noise_sigma = 0.2)

plot_images(list(im4$original, im4$corrupted), titles=c("Original", "Corrupt"))


#
im5 <- generate_example_images(pattern= "disk", mask_type = "disk", missing_fraction = 0.1, noise_sigma = 0.1)
plot_images(list(im5$original, im5$corrupted), titles=c("Original", "Corrupt"))


# Re2: Im2 restored

re2 <- tv_restore(im2$corrupted, mask=im2$mask, lambda = 1e-12, tau = 1e6, sigma = 1e-1, theta = 1e0)
plot_images(list(re2$restored, im2$corrupted, im2$corrupted-re2$restored), titles = c("Restored", "Corrupted", "Diff"))
  diff <- im2$corrupted-re2$restored
max(diff)
max(re2$restored)
a <- re2$restored
image(t(re2$restored[nrow(re2$restored):1,]),
      col = pal,
      useRaster = TRUE,
      main = "Restored Image")
b <- im2$corrupted


re3 <- tv_restore(im3$corrupted, mask=im3$mask, lambda = 1e-9, tau = 0.125, sigma = 0.125, theta = 1)
plot_images(list(re3$restored, im3$original, im3$corrupted), titles = c("Restored", "Original", "Corrupted"))
c <- im3$corrupted
