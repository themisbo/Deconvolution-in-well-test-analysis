
par(mfrow = c(2,3))

hist(Posterior[,1], xlim = c(1, 4), freq = F, col = "red", main = "TM", xlab = "value")
lines(seq(1, 4, length.out = 10000), dnorm(seq(1, 4, length.out = 10000), 2, .2), col = "blue" )

hist(Posterior[,2], xlim = c(1.5, 3.5), freq = F, col = "red", main = "PM", xlab = "value")
lines(seq(1.5, 3.5, length.out = 10000), dnorm(seq(1.5,3.5, length.out = 10000), 1.5, .2), col = "blue" )

hist(Posterior[,3], xlim = c(0, 50), freq = F, col = "red", main = "CDe2Skin", xlab = "value")
lines(seq(0, 50, length.out = 10000), dgamma(seq(0, 50, length.out = 10000), 1, .2), col = "blue" )

hist(Posterior[,4], xlim = c(0, 10), freq = F, col = "red", main = "RD1", xlab = "value")
lines(seq(0, 10, length.out = 10000), dnorm(seq(0, 10, length.out = 10000), 2, 1), col = "blue" )

hist(Posterior[,Region+4], xlim = c(-7, 7), freq = F, col = "red", main = "eta1", xlab = "value")
lines(seq(-7, 7, length.out = 10000), dnorm(seq(-7, 7, length.out = 10000), 0, 1), col = "blue" )

hist(Posterior[,2*Region+4], xlim = c(-7, 7), freq = F, col = "red", main = "M1", xlab = "value")
lines(seq(-7, 7, length.out = 10000), dnorm(seq(-7, 7, length.out = 10000), 0, 1), col = "blue" )
