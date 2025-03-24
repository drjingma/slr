# About

`slr` is an R package for selection of balances in compositional data analysis. It implements a supervised log ratio method for the identification of two groups of variables whose log ratio is associated with a continuous or binary response variable. 

# Getting Started

## How to install

`devtools::install_github(repo = "drjingma/slr",build_vignettes=T)`

## How to use

Once installed, you can call the library and check out the help functions

`library(slr)`

`?slr`

You can read the associated vignette in R directly if `build_vignettes=T`:  

`browseVignettes(package='slr')`

You can also find the vignette [here](https://htmlpreview.github.io/?https://github.com/drjingma/slr/blob/master/vignettes/slr.html). In this vignette, you will see that `slr` performs comparably to `selbal` on an HIV data set and runs much faster than `selbal` on a UC data set where the number of features is over 400. 

## Note

- In the `cv.slr` function, if `plot=TRUE`, the function will display a plot of the cross-validation deviance measure as a function of the tuning parameter. The error bars of the CV deviance are also displayed. However, at the boundary of tuning parameters, the length of the error bar could be zero, which means that the standard error of the CV deviance across folds is 0. This may render the warning message: 

        1: In doTryCatch(return(expr), name, parentenv, handler) :
          zero-length arrow is of indeterminate angle and so skipped

- The current implementation of `slr` works for continuous, binary, or survival response variables. 
