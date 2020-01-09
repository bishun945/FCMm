
# FCMm <a href='https://github.com/bishun945/FCM-m'><img src='man/figures/logo.png' align='right' height="139" /></a>

**Author**: Shun Bi  
**Date**: 2020/01/09  
**E-mail**: <bishun1994@foxmail.com>

## Overview

`FCMm` is a package for fuzzy clustering water spectra (or called water
color). Given that the most of water color spectra data sets are
considered as the high dimensional set, the advantage of this method is
making FCM assign the membership (sum as 1) harder, ensuring the desired
water type are restricted to its belongings (not too soft). It is
possible to cluster the harm algal bloom water type which can not be
produced by FCM with `m=2`.

  - If you want to cluster your own data sets, it provides an improved
    Fuzzy Cluster Method (FCM) by optimizing the fuzzifier value
    (default but not good being 2).
  - You can also use the built-in cluster of inland waters produced by
    [Bi *et al.*
    (2019)](https://www.osapublishing.org/oe/abstract.cfm?uri=oe-27-24-34838)
    and can simply obtain the Chlorophyll-a concentration by blending
    three algorithms with relatively low bias.  
  - It supports raster (or called imagery) processing (see more details
    in help documents or vignettes).
  - It includes several data sets about water color spectra and
    corresponding water quality parameters and a testing image raster
    (see help documents for details).

## Installation

``` r
# Install the FCMm by typing (Dont RUN it right now):
install.packages("FCMm")
# To get bugs fixed version, you can install by (Use this one):
devtools::install_github("bishun945/FCMm")
```

## Usage

``` r
# Load testing data
library(FCMm)
library(ggplot2)
data("WaterSpec35")
data("Bi_clusters")
Rrs <- WaterSpec35[,3:17]
# Plot the spectra
plot_spec_from_df(Rrs) + 
  labs(x='Wavelength (nm)',y=expression(Rrs~(sr^-1))) + 
  theme_bw() + 
  theme(legend.position='none', text=element_text(size=18))
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="60%" style="display: block; margin: auto;" />

``` r
# Applying FCMm
result <- apply_FCM_m(Rrs=Rrs, option.plot=TRUE)
plot(result$p.group)
```

<img src="man/figures/README-unnamed-chunk-2-2.png" width="60%" style="display: block; margin: auto;" />

## Getting help

  - About this package, I have written four vignettes to present the
    usage of `FCMm`. Please read them carefully if you want to use this
    package for your research. Also, e-mail me via
    [bishun1994@foxmail.com](bishun1994@foxmail.com) without hesitation
    if you have any questions or find any bug about it. They are:
      - [Cluster\_number\_determination](./doc/Cluster_number_determination.html)
      - [Fuzzy\_cluster\_training](./doc/Fuzzy_cluster_training.html)
      - [New\_Data\_Running\_FCMm](./doc/New_Data_Running_FCMm.html)
      - [Raster\_Running\_FCMm](./doc/Raster_Running_FCMm.html)
  - If you are more interested in the application of FCM-m about inland
    water spectra, I recommend you to read [Bi *et al.*
    (2019)](https://www.osapublishing.org/oe/abstract.cfm?uri=oe-27-24-34838)
    for more details.
  - If you want to know some theoretical knowledge about FCM in
    mathematics, you could read some researches like [Dembele *et al.*
    (2018)](https://link.springer.com/article/10.1007/s11634-008-0032-5).
  - More about FCM in remote sensing applications, you can read [Moore
    *et al.*
    (2014)](https://www.sciencedirect.com/science/article/pii/S0034425713004434)
    and [Jackson *et al.*
    (2017)](https://www.sciencedirect.com/science/article/pii/S0034425717301396)
    which focus on Case-II and Case-I waters, respectively.
  - See more details about optical water types of inland waters in
    [Spyrakos *et al.*
    (2018)](https://aslopubs.onlinelibrary.wiley.com/doi/abs/10.1002/lno.10674)
  - Hope you will enjoy using this package and have a nice day.
