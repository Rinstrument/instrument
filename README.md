<!-- Improved compatibility of back to top link: See: https://github.com/othneildrew/Best-README-Template/pull/73 -->
<a name="readme-top"></a>
<!-- PROJECT LOGO -->

---

<!-- <p align="right">
  <a href="https://github.com/inirt/theta2">
    <img src="https://github.com/inirt/.github/blob/master/images/hex-inirt.png" alt="Logo" width="80" height="80">
  </a>
</p> -->

<img align="right" src="https://github.com/Rinstrument/instrument/blob/master/www/hexsticker.png" width="200px">

# R instrument quick guide

## Overview

The R package `instrument` is an item response theory modeling software whose purpose is: 

 - Fit a variety of IRT models including the univdimensional, multidimensional, and higher-order models
 - Specify regression models in the context of IRT with both fixed and random effects (i.e., mixed modeling)
 - Simple model syntax to describe IRT models with and without regression

## Documentation & Source

 - The `instrument` R package source code can be found at [github.com/inirt/theta2](https://github.com/inirt/theta2)

 - Full documentation and tutorials can be found at [inirt.github.io/doc/](https://inirt.github.io/doc/)

---

## Installation

1. Since this is an R pacakge, the user first needs to install R from <a href="https://www.r-project.org/">the R project website</a>.

2. In addition to R, a good editor is recommended such as RStudio.

3. The key dependency of this package is RStan. Use the development version of the software, which can be installed with:

``` r
install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
```

4. Once R and RStan are installed, open an R console and install the package from Github using remotes (recommended):

``` r
# remotes
install.packages("remotes")
# install from Github
remotes::install_github("Rinstrument/instrument")
```

Alternatively, install from CRAN with (coming soon):

``` r
install.packages("instrument")
```

---
