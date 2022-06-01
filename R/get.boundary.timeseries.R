get.boundary.timeseries <- function(data, type = "TLSW") {

  # this function takes a time series and produces a 4*data.len length series that adds
  # boundary handling to both sides.

  data.len <- length(data)
  J <- wavethresh::IsPowerOfTwo(data.len)

  s <- seq(from = 0, to = (data.len - 1) / data.len, length = data.len)

  L <- lm(data ~ poly(s, 3, raw = TRUE))

  bh.right <- predict(L, newdata = data.frame(s = 1))

  bh.left <- predict(L, newdata = data.frame(s = -1 / data.len))

  if (type == "LSW.diff") {
    bh.series1 <- c(data - bh.right + bh.left, data, data + bh.right - bh.left)
  } else {
    bh.series1 <- c(-rev(data) + 2 * bh.left, data, -rev(data) + 2 * bh.right)
  }

  if (is.na(J) == TRUE) {
    l <- 2^floor(log2(length(bh.series1)))
    k <- floor((3 * data.len - l) / 2)
    if (data.len %% 2 == 0) {
      bh.series2 <- bh.series1[(k + 1):(3 * data.len - k)]
    } else {
      bh.series2 <- bh.series1[(k + 1):(3 * data.len - k - 1)]
    }

    return(bh.series2)
  } else {
    l <- length(bh.series1)

    bh.series2 <- c((data)[(data.len / 2 + 1):data.len] - abs(2 * bh.right - 2 * bh.left), bh.series1, (data)[1:(data.len / 2)] + abs(2 * bh.right - 2 * bh.left))

    return(bh.series2)
  }
}
