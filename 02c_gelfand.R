library(LaplacesDemon)

gelfandplot<- function (x, k = 3, pen = FALSE, xlab = xlab) 
{
  if (missing(x)) 
    stop("The x argument is required.")
  if (!is.vector(x)) 
    x <- as.vector(x)
  if (k < 2) 
    k <- 2
  if (k > length(x)/2) 
    k <- round(length(x)/2)
  if ({
    length(x)/k
  } < 2) 
    stop("k is too large relative to length(x).")
  quantiles <- seq(from = 0, to = 1, by = 1/k)
  breaks <- round(as.vector(quantiles) * length(x))
  breaks <- breaks[-1]
  d.temp <- stats::density(x)
  d <- array(c(d.temp$x, d.temp$y), dim = c(length(d.temp$x), 
                                            2, length(breaks)))
  d.temp <- stats::density(x[1:breaks[1]])
  d[, , 1] <- c(d.temp$x, d.temp$y)
  for (i in 2:length(breaks)) {
    d.temp <- stats::density(x[1:breaks[i]])
    d[, , i] <- c(d.temp$x, d.temp$y)
  }
  ymax <- max(d[, 2, ])
  col.list <- c("red", "green", "blue", "yellow", "purple", 
                "orange", "brown", "gray", "burlywood", "aquamarine")
  col.list <- rep(col.list, len = length(breaks))
  rgb.temp <- as.vector(col2rgb(col.list[1]))
  mycol <- rgb(red = rgb.temp[1], green = rgb.temp[2], blue = rgb.temp[3], 
               alpha = 50, maxColorValue = 255)
  plot(d[, 1, 1], d[, 2, 1], type = "l", col = mycol, xlim = c(range(d[, 
                                                                       1, ])), ylim = c(0, ymax), main = "Gelfand Diagnostic", 
       xlab = xlab, ylab = "Density")
  polygon(x = d[, 1, 1], y = d[, 2, 1], col = mycol, border = NULL)
  for (i in 2:length(breaks)) {
    rgb.temp <- as.vector(col2rgb(col.list[i]))
    mycol <- rgb(red = rgb.temp[1], green = rgb.temp[2], 
                 blue = rgb.temp[3], alpha = 50, maxColorValue = 255)
    lines(d[, 1, i], d[, 2, i], col = mycol)
    polygon(x = d[, 1, i], y = d[, 2, i], col = mycol, border = mycol)
    lines(d[, 1, i], d[, 2, i], lty = i)
  }
  if (pen == TRUE) 
    abline(v = mean(range(d[, 1, ])), col = "black", lwd = 10)
  return(invisible(x))
}

nums = c(1:4,Region +4, 2*Region+4)

par(mfrow = c(2,3))
for(i in nums) 
  gelfandplot(Posterior[,i], xlab = colnames(Posterior)[i])
