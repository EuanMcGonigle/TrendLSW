spec.plot <- function (x, xlabvals, xlabchars, ylabchars, first.level = 0, n,
          main = "Spectral Estimate", scaling = "global",
          rhlab = FALSE, sub = "", NotPlotVal = 0.005, xlab = "Time",
          ylab = "Scale", aspect = "Identity", ...)
{
  ctmp <- class(x)
  if (is.null(ctmp))
    stop("wd has no class")
  else if (ctmp != "wd")
    stop("wd is not of class wd")
  levels <- nlevelsWT(x)
  nlevels <- levels - first.level
  type <- x$type

  if (aspect != "Identity"){
    sub <- paste(sub, "(", aspect, ")")
  }

  plot(c(0, 0, n, n), c(0, nlevels + 1, nlevels + 1, 0), type = "n",
       xlab = xlab, ylab = ylab, main = main, yaxt = "n", xaxt = "n",
       sub = sub, ...)

  if(missing(ylabchars)){
    ylabchars <- 1:nlevels
  }

  axis(2, at = 1:(nlevels), labels = ylabchars)
  if (missing(xlabchars)) {
    if (missing(xlabvals)) {
       axx <- c(0, 2^(levels - 2), 2^(levels - 1),
                    2^(levels - 1) + 2^(levels - 2), 2^levels)
      axis(1, at = axx)
    }
    else {
      lx <- pretty(xlabvals, n = 4)
      cat("lx is ", lx, "\n")
      if (lx[1] < min(xlabvals))
        lx[1] <- min(xlabvals)
      if (lx[length(lx)] > max(xlabvals))
        lx[length(lx)] <- max(xlabvals)
      cat("lx is ", lx, "\n")
      xix <- NULL
      for (i in 1:length(lx)) {
        u <- (xlabvals - lx[i])^2
        xix <- c(xix, (1:length(u))[u == min(u)])
      }
      axx <- xix
      axl <- signif(lx, digits = 2)
      axis(1, at = axx, labels = axl)
    }
  }
  else axis(1, at = xlabvals, labels = xlabchars)
  myxx <- 1:n
  height <- 1
  first.last.d <- x$fl.dbase$first.last.d
  axr <- NULL
  if (scaling == "global") {
    my <- 0
    for (i in ((levels - 1):first.level)) {
      y <- accessD(x, i, aspect = aspect)
      my <- max(c(my, abs(y)))
    }
  }
  if (scaling == "compensated") {
    my <- 0
    for (i in ((levels - 1):first.level)) {
      y <- accessD(x, i, aspect = aspect) * 2^(i/2)
      my <- max(c(my, abs(y)))
    }
  }
  if (scaling == "super") {
    my <- 0
    for (i in ((levels - 1):first.level)) {
      y <- accessD(x, i, aspect = aspect) * 2^i
      my <- max(c(my, abs(y)))
    }
  }
  shift <- 1
  for (i in ((levels - 1):first.level)) {
    y <- accessD(x, i, aspect = aspect)
    if (type == "wavelet")
      n <- 2^i
    else {
      y <- y[c((n - shift + 1):n, 1:(n - shift))]
      shift <- shift * 2
    }
    xplot <- myxx
    ly <- length(y)
    if (scaling == "by.level")
      my <- max(abs(y))
    if (scaling == "compensated")
      y <- y * 2^(i/2)
    if (scaling == "super")
      y <- y * 2^i
    if (my == 0) {
      y <- rep(0, length(y))
    }
    else y <- (0.5 * y)/my
    axr <- c(axr, my)
    if (max(abs(y)) > NotPlotVal)
      segments(xplot, height, xplot, height + y)
    if (i != first.level) {
      if (type == "wavelet") {
        x1 <- myxx[seq(1, n - 1, 2)]
        x2 <- myxx[seq(2, n, 2)]
        myxx <- (x1 + x2)/2
      }
      height <- height + 1
    }
  }
  if (rhlab == TRUE)
    axis(4, at = 1:length(axr), labels = signif(axr, digits = 3))
  axr
}
