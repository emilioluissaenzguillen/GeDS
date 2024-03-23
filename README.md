# GeDS
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/GeDS)](https://cran.r-project.org/package=GeDS) ![downloads](https://cranlogs.r-pkg.org/badges/grand-total/GeDS)

Geometrically Designed Spline ('GeDS') Regression is a non-parametric geometrically motivated method for fitting variable knots spline predictor models in one or two independent variables, in the context of generalized (non-)linear models. 'GeDS' estimates the number and position of the knots and the order of the spline, assuming the response variable has a distribution from the exponential family. A description of the method can be found in [Kaishev et al. (2016)](https://link.springer.com/article/10.1007/s00180-015-0621-7) and [Dimitrova et al. (2017)](https://www.sciencedirect.com/science/article/pii/S0096300322005677?via%3Dihub).

Installation:
=============

To install the stable version on R CRAN:

``` r
    install.packages("GeDS")
```

To install the latest development version:

``` r
    install.packages("devtools")
    devtools::install_github("emilioluissaenzguillen/GeDS")
```

License:
========

This package is free and open source software, licensed under GPL-3

