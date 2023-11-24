
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TrendLSW

<!-- badges: start -->

[![R-CMD-check](https://github.com/EuanMcGonigle/TrendLSW/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/EuanMcGonigle/TrendLSW/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/EuanMcGonigle/TrendLSW/branch/main/graph/badge.svg)](https://app.codecov.io/gh/EuanMcGonigle/TrendLSW?branch=main)
<!-- badges: end -->

Implements wavelet methods for analysis of nonstationary time series.
See

> McGonigle, E. T., Killick, R., and Nunes, M. (2022). Trend locally
> stationary wavelet processes. *Journal of Time Series Analysis*,
> 43(6), 895-917.

> McGonigle, E. T., Killick, R., and Nunes, M. (2022). Modelling
> time-varying first and second-order structure of time series via
> wavelets and differencing. *Electronic Journal of Statistics*, 6(2),
> 4398-4448.

for full details.

## Installation

To install `TrendLSW` from GitHub:

    devtools::install_github("https://github.com/EuanMcGonigle/TrendLSW")

## Usage

For detailed examples, see the help files within the package. We can
generate a small example for performing trend estimation as follows:

    set.seed(1)

    noise <- rnorm(512)
    trend <- seq(from = 0, to = 5,length = 512)
    x <- trend + noise

Apply the `TLSW` function:

    x.TLSW <- TLSW.est(x)

Visualise the estimated trend and spectrum:

    plot(x.TLSW)
