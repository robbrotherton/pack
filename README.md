
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pack

<!-- badges: start -->
<!-- badges: end -->

This package provides some shape-packing functions. It’s a work in
progress, full of bugs and inefficiencies.

## Installation

You can install the development version of pack like so:

``` r
pak::pak("robbrotherton/pack")
```

## Example

``` r
library(pack)
library(ggplot2)
#> Warning: package 'ggplot2' was built under R version 4.2.2

set.seed(123)

container <- draw::polygon(sides = 6, radius = 100)

radii <- c(rep(30, 2),
           rep(20, 4),
           rep(10, 20),
           rep(5, 200), 
           rep(3, 300),
           rep(2, 500),
           rep(1, 500),
           rep(0.5, 1000))

circles <- pack_circles_all(container, 
                        radii = radii, 
                        existing_circles = data.frame(x = 0, y = 0, r = 0)) |> 
  dplyr::mutate(fill = factor(round(runif(dplyr::n(), min = 1, max = 4), 0)))
#> made 1399

ggplot(circles) +
  geom_polygon(data = container, aes(x, y)) +
  ggforce::geom_circle(aes(x0 = x, y0 = y, r = r, fill = fill)) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "none")
#> Warning: Using the `size` aesthetic in this geom was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` in the `default_aes` field and elsewhere instead.
```

<img src="man/figures/README-example-1.png" width="100%" />
