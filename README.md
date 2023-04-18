# TrendLSW

Implements wavelet methods for analysis of nonstationary time series. See 

> McGonigle, E. T., Killick, R., and Nunes, M. (2022). Trend locally stationary wavelet processes. *Journal of Time Series Analysis*.
> 
> McGonigle, E. T., Killick, R., and Nunes, M. (2022). Modelling time-varying first and second-order structure of time series via wavelets and differencing. *Electronic Journal of Statistics* [link to paper](https://projecteuclid.org/journals/electronic-journal-of-statistics/volume-16/issue-2/Modelling-time-varying-first-and-second-order-structure-of-time/10.1214/22-EJS2044.full).

for full details.

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

````

Apply `wav.trend.est`:
```
trend.est = wav.trend.est(x, filter.number = 4, family = "DaubLeAsymm", boundary.handle =TRUE)
```

Visualise the estimated trend function and underlying truth:
```
plot.ts(x, lty = 1, col = 8)
lines(trend, col = 2, lwd = 2)
lines(trend.est$trend.est,col = 4, lwd = 2, lty = 2)
```


