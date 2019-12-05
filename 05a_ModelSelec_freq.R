library(mdeconv)
library(tidyverse)
library(BayesianTools)

load("oil1r_DEzs.RData")
Posterior_1 <- getSample(out, start = 30001, coda = FALSE, parametersOnly = FALSE)

load("oil2r_DEzs.RData")
Posterior_2 <- getSample(out, start = 30001, coda = FALSE, parametersOnly = FALSE)

load("oil3r_DEzs.RData")
Posterior_3 <- getSample(out, start = 30001, coda = FALSE, parametersOnly = FALSE)

load("oil4r_DEzs.RData")
Posterior_4 <- getSample(out, start = 30001, coda = FALSE, parametersOnly = FALSE)

# AIC
AIC <- function(x) min(2*(ncol(x) - 3) -2*x[, ncol(x)-2])
AIC_1r <- AIC(Posterior_1)
AIC_2r <- AIC(Posterior_2)
AIC_3r <- AIC(Posterior_3)
AIC_4r <- AIC(Posterior_4)

AICs <- c(AIC_1r, AIC_2r, AIC_3r, AIC_4r)
plot(AICs, ylab = "AIC", xlab = "Region")
lines(AICs)
points(which.min(AICs), min(AICs), col = 2, pch=19)

# BIC
BIC <- function(x) min(log(nrow(x))*(ncol(x) - 3) -2*x[, ncol(x)-2])
BIC_1r <- BIC(Posterior_1)
BIC_2r <- BIC(Posterior_2)
BIC_3r <- BIC(Posterior_3)
BIC_4r <- BIC(Posterior_4)

BICs <- c(BIC_1r, BIC_2r, BIC_3r, BIC_4r)
plot(BICs, ylab = "BIC", xlab = "Region")
lines(BICs)
points(which.min(BICs), min(BICs), col = 2, pch=19)

# DIC
DIC <- function(x) {
  Dev = -2*x[, ncol(x)-2]
  DIC <- mean(Dev) + var(Dev)/2
  return(DIC)
}

DIC_1r <- DIC(Posterior_1)
DIC_2r <- DIC(Posterior_2)
DIC_3r <- DIC(Posterior_3)
DIC_4r <- DIC(Posterior_4)

DICs <- c(DIC_1r, DIC_2r, DIC_3r, DIC_4r)
plot(DICs, ylab = "DIC", xlab = "Region")
lines(DICs)
points(which.min(DICs), min(DICs), col = 2, pch=19)

xtable(tibble(AIC = AICs, BIC = BICs, DIC = DICs) %>%
         rownames_to_column %>% 
         gather(var, value, -rowname) %>% 
         spread(rowname, value) )

