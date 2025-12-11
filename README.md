
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Package `sparvaride`

The package implements the variance identification algorithm for sparse
factor analysis described in the paper “Cover It Up\! Bipartite Graphs
Uncover Identifiability in Sparse Factor Analysis” by Darjus Hosszejni
and Sylvia Frühwirth-Schnatter. The paper is published in the [Journal
of Multivariate
Analysis](https://doi.org/10.1016/j.jmva.2025.105536).

The package is still under development and the API is subject to change.
For a Matlab implementation, see [`sparvaride-matlab`](https://github.com/hdarjus/sparvaride-matlab).

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
         nrow = 7, ncol = 3,
         byrow = TRUE)
delta2 <-
  matrix(c(1, 0, 0,
           0, 1, 0,
           0, 0, 1,
           1, 1, 1,
           1, 0, 1,
           1, 1, 1,
           1, 0, 1),
         nrow = 7, ncol = 3,
         byrow = TRUE)
```

Then, we call the `counting_rule_holds` function on these matrices:

``` r
counting_rule_holds(delta1)
#> [1] FALSE
counting_rule_holds(delta2)
#> [1] TRUE
```

## Citation

For citing our work, please check the `citation` function in R:

``` r
citation("sparvaride")
#> 
#> To cite sparvaride in publications use:
#> 
#>   Hosszejni D, Frühwirth-Schnatter S (2026). "Cover It Up! Bipartite
#>   Graphs Uncover Identifiability in Sparse Factor Analysis." _Journal
#>   of Multivariate Analysis_, *211*, 105536. ISSN 0047-259X,
#>   doi:10.1016/j.jmva.2025.105536
#>   <https://doi.org/10.1016/j.jmva.2025.105536>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {Cover It Up! Bipartite Graphs Uncover Identifiability in Sparse Factor Analysis},
#>     author = {Darjus Hosszejni and Sylvia Frühwirth-Schnatter},
#>     journal = {Journal of Multivariate Analysis},
#>     year = {2026},
#>     volume = {211},
#>     pages = {105536},
#>     issn = {0047-259X},
#>     doi = {10.1016/j.jmva.2025.105536},
#>   }
```
