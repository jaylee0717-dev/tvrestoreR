devtools::load_all(".")

################################################################

# Im3: disk with blocks
im3 <- generate_example_images(size= c(128,128), pattern= "wedge", mask_type = "disk", missing_fraction = 0., noise_sigma = 0.1)
plot_images(list(im3$original, im3$corrupted), titles=c("Original", "Corrupt"))

re3_2 <- tv_restore(im3$corrupted, mask=im3$mask, task="inpainting",
                    lambda = 1, max_iter = 20000, tol = 5e-04, verbose = TRUE)
plot_images(list(re3_2, im3$original, im3$corrupted), titles = c("Restored", "Original", "Corrupted"))


re3 <- tv_restore(im3$corrupted, mask=im3$mask, task="ROF_inpainting",
                  lambda = 1, max_iter = 20000, tol = 5e-04, verbose = TRUE)
plot_images(list(re3, im3$original, im3$corrupted), titles = c("Restored", "Original", "Corrupted"))


re3_3 <- tv_restore(im3$corrupted, mask=im3$mask, task="denoising",
                    lambda = 0.001, max_iter = 20000, tol = 5e-04, verbose = TRUE)
plot_images(list(re3_3, im3$original, im3$corrupted), titles = c("Restored", "Original", "Corrupted"))


# Im4: find lambda tests
im4 <- generate_example_images(size= c(128,128), pattern= "wedge", mask_type = "disk", missing_fraction = 0., noise_sigma = 0.1)
plot_images(list(im4$original, im4$corrupted), titles=c("Original", "Corrupt"))

out <- find_lambda(im4$corrupted, im4$mask, task="denoising", lambda_grid = 10^seq(-2.5, 1, length.out = 3),verbose=TRUE, ground_truth = 0.1)
plot_images(list(out$restored, im3$original, im3$corrupted), titles = c("Restored", "Original", "Corrupted"))

##############################################################################
# Test Split Bregman with im3
im3 <- generate_example_images(size= c(128,128), pattern= "wedge", mask_type = "disk", missing_fraction = 0.15, noise_sigma = 0.1)
plot_images(list(im3$original, im3$corrupted), titles=c("Original", "Corrupt"))

re3_3 <- tv_restore(im3$corrupted, mask=im3$mask, task="denoising", method="split_bregman",
                    lambda = 1, max_iter = 20000, tol = 1e-04, verbose = TRUE)
plot_images(list(re3_3, im3$original, im3$corrupted), titles = c("Restored", "Original", "Corrupted"))


re3_2 <- tv_restore(im3$corrupted, mask=im3$mask, task="inpainting", method = "split_bregman",
                    lambda = 1, max_iter = 20000, tol = 5e-04, verbose = TRUE)
plot_images(list(re3_2, im3$original, im3$corrupted), titles = c("Restored", "Original", "Corrupted"))


re3 <- tv_restore(im3$corrupted, mask=im3$mask, task="ROF_inpainting",method = "split_bregman",
                  lambda = 10, max_iter = 20000, tol = 5e-04, verbose = TRUE)
plot_images(list(re3, im3$original, im3$corrupted), titles = c("Restored", "Original", "Corrupted"))



#############################
# Test Image List input
im3 <- generate_example_images(size= c(128,128), pattern= "wedge", mask_type = "disk", missing_fraction = 0.1, noise_sigma = 0.1)
im4 <- generate_example_images(size= c(128,128), pattern= "gradient", mask_type = "disk", missing_fraction = 0.15, noise_sigma = 0.2)

a <- list(im3$corrupted,im4$corrupted)
re_list <- tv_restore(a, mask=im3$mask, task="ROF_inpainting", method = "primal_dual",
                               lambda = 1, max_iter = 20000, tol = 5e-04, verbose = TRUE)

plot_images(list(re_list[[1]], im3$original, im3$corrupted), titles = c("Restored", "Original", "Corrupted"))
plot_images(list(re_list[[2]], im4$original, im4$corrupted), titles = c("Restored", "Original", "Corrupted"))



###################################
# Test Lambda finding
best <- find_lambda(im3$corrupted, mask = im3$mask, lambda_grid = 10^seq(0, 1, length.out=3),
            task="inpainting", method = "split_bregman", ground_truth = 0.1)
plot_images(list(best$restored, im3$original, im3$corrupted), titles = c("Restored", "Original", "Corrupted"))
