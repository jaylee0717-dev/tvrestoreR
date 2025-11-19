
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tvrestoreR

<!-- badges: start -->

<!-- badges: end -->

The goal of tvrestoreR is to provide for a unified, flexible package
that handles both denoising and inpainting of grayscale images under TV
regularization, with good performance, diagnostics, and ease of use.

## Work In Progress – What is Left

So Far: functions for creating and plotting example images What is Left:
Implementation of restore function, and gridsearch on lamba function.

## Installation

You can install the development version of tvrestoreR from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("jaylee0717-dev/tvrestoreR")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(tvrestoreR)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
