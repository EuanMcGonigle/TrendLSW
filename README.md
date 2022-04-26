# TrendLSW
TrendLSW R Package
# fvarseg
Implements wavelet methods for analysis of nonstionary time series. See 

> Trend locally stationary wavelet processes
> Modelling time-varying first and second-order structure of time Series via wavelets and differencing

by Euan T. McGonigle, Rebecca Killick and Matthew A. Nunes. See [here](https://onlinelibrary.wiley.com/doi/10.1111/jtsa.12643) and [arXiv:2108.07750](https://arxiv.org/abs/2108.07550) for full details.

## Installation

To install `TrendLSW` from GitHub:

```
devtools::install_github("https://github.com/EuanMcGonigle/TrendLSW")
```

## Usage

For detailed examples, see the help files within the package. We can generate a small example for performing trend estimation as follows:

```
set.seed(1)

noise = rnorm(512)
trend = seq(from = 0, to = 5,length = 512)
x = trend+noise

plot.ts(x, lty = 1, col = 8)
lines(trend, col = 2, lwd = 2)
lines(trend.est,col = 4, lwd = 2, lty = 2)
````

Apply `wav.trend.est`:
```
trend.est = wav.trend.est(x, filter.number = 4, family = "DaubLeAsymm", boundary.handle =TRUE)
```

Visulaise the estimated trens function and underlying truth:
```
plot.ts(x, lty = 1, col = 8)
lines(trend, col = 2, lwd = 2)
lines(trend.est,col = 4, lwd = 2, lty = 2)
```


