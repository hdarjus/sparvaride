
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Package `sparvaride`

<!-- badges: start -->

[![R-CMD-check](https://github.com/hdarjus/econometric.factor.identification/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/hdarjus/econometric.factor.identification/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The package implements the variance identification algorithm for sparse
factor analysis described in the paper “Cover It Up\! Bipartite Graphs
Uncover Identifiability in Sparse Factor Analysis” by Darjus Hosszejni
and Sylvia Frühwirth-Schnatter. The paper is available at
[arXiv](https://arxiv.org/abs/2211.00671).

The package is still under development and the API is subject to change.

## Installation

You can install the development version of `sparvaride` from
[GitHub](https://github.com/hdarjus/sparvaride) with:

``` r
# install.packages("devtools")
devtools::install_github("hdarjus/sparvaride")
```

## The `counting_rule_holds` Function

We can check whether the 3579 counting rule holds for a given binary
matrix `delta` using the `counting_rule_holds` function in the
`sparvaride` package.

``` r
library(sparvaride)
```

We define two matrices as above in R:

``` r
delta1 <-
  matrix(c(1, 0, 0,
           0, 1, 0,
           0, 0, 1,
           1, 1, 1,
           1, 0, 1,
           1, 0, 1,
           1, 0, 1),
         nrow = 6, ncol = 3,
         byrow = TRUE)
#> Warning in matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, :
#> data length [21] is not a sub-multiple or multiple of the number of rows [6]
delta2 <-
  matrix(c(1, 0, 0,
           0, 1, 0,
           0, 0, 1,
           1, 1, 1,
           1, 0, 1,
           1, 1, 1,
           1, 0, 1),
         nrow = 6, ncol = 3,
         byrow = TRUE)
#> Warning in matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, :
#> data length [21] is not a sub-multiple or multiple of the number of rows [6]
```

Then, we call the `counting_rule_holds` function on these matrices:

``` r
counting_rule_holds(delta1)
#> [1] FALSE
counting_rule_holds(delta2)
#> [1] FALSE
```

## Citation

For citing our work, please check the `citation` function in R:

``` r
citation("sparvaride")
#> 
#> To cite sparvaride in publications use:
#> 
#>   Hosszejni D, Frühwirth-Schnatter S (2022). "Cover It Up! Bipartite
#>   Graphs Uncover Identifiability in Sparse Factor Analysis."
#>   doi:10.48550/arXiv.2211.00671
#>   <https://doi.org/10.48550/arXiv.2211.00671>, arXiv: 2211.00671.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Unpublished{,
#>     title = {Cover It Up! Bipartite Graphs Uncover Identifiability in Sparse Factor Analysis},
#>     author = {Darjus Hosszejni and Sylvia Frühwirth-Schnatter},
#>     year = {2022},
#>     note = {arXiv: 2211.00671},
#>     doi = {10.48550/arXiv.2211.00671},
#>   }
```
