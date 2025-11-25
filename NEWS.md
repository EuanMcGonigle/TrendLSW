# TrendLSW 1.0.4

* Updated the example in the TLSW() function to a smaller sample size to comply with CRAN checks.

# TrendLSW 1.0.3

* The default values for S.filter.number and S.family in the `TLSW()` function 
have been changed to T.filer.number and T.family respectively.

* The default value for the `plot.CI` argument in `plot.TLSW()` has now been removed. 
If a CI was calculated for the trend, it will be plotted by default, unless the
user sets `plot.CI = FALSE`. In addition, `plot.TLSW()` now warns the user if 
argument `plot.CI = TRUE` but no confidence interval was calculated. 

* Fixed issue in `plot.TLSW()` where data was not displayed in trend plot when 
confidence interval not requested.

# TrendLSW 1.0.2

* Added new data: z.acc and z.labels. See the associated help files for more 
details.

# TrendLSW 1.0.1

* Updated description field in Description file.
* Updated documentation for plot.TLSW.

# TrendLSW 1.0.0

* Initial CRAN submission.
