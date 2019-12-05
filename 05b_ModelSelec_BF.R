# Erase everything in memory
rm(list=ls())

library(BayesianTools)
library(mdeconv)
library(Matrix)

load("real_TLS.RData")
obsPress <- real_TLS$swPress
obsRate    <- real_TLS$swRate
trueResp    <- real_TLS$TlsResp
which_dd <- which(obsRate$Rate != 0)

q = obsRate$Rate[which_dd]

p0 <- max(obsPress$Press)
sigmaq = .05*max(q)
sigmap0 = 10

N = length(which_dd)
m = length(obsPress$Press)
p = obsPress$Press
N_Nodes = 30
t0 = 1e-3
tN = ceiling( obsPress[nrow(obsPress)]$Time )
tau = seq(log(t0), log(tN), length.out = N_Nodes)
t = exp(tau)

load("oil1r_DEzs.RData")
start_value = 30001
Posterior <- getSample(out, start = start_value, coda = FALSE) 
Region = (ncol(Posterior) - 4)/3
source("functions.R")
out1 <- out
M1 = marginalLikelihood(out1, start = start_value)

load("oil2r_DEzs.RData")
start_value = 30001
Posterior <- getSample(out, start = start_value, coda = FALSE) 
Region = (ncol(Posterior) - 4)/3
source("functions.R")
out2 <- out
M2 = marginalLikelihood(out2, start = start_value)

load("oil3r_DEzs.RData")
start_value = 30001
Posterior <- getSample(out, start = start_value, coda = FALSE) 
Region = (ncol(Posterior) - 4)/3
source("functions.R")
out3 <- out
M3 = marginalLikelihood(out3, start = start_value)

load("oil4r_DEzs.RData")
start_value = 30001
Posterior <- getSample(out, start = start_value, coda = FALSE) 
Region = (ncol(Posterior) - 4)/3
source("functions.R")
out4 <- out
M4 = marginalLikelihood(out4, start = start_value)

BF_12r = exp(M1$ln.ML - M2$ln.ML )
BF_13r = exp(M1$ln.ML - M3$ln.ML )
BF_14r = exp(M1$ln.ML - M4$ln.ML )
BF_23r = exp(M2$ln.ML - M3$ln.ML )
BF_24r = exp(M2$ln.ML - M4$ln.ML )
BF_34r = exp(M3$ln.ML - M4$ln.ML )

BFs <- c(BF_12r, BF_13r, BF_14r, BF_23r, BF_24r, BF_34r)
BFs

df <- data.frame(B10 = c("21", "31", "41", "32", "42", "43"),
                 outcome = 2*c(M2$ln.ML - M1$ln.ML, M3$ln.ML - M1$ln.ML, M4$ln.ML - M1$ln.ML, M3$ln.ML - M2$ln.ML,
                               M4$ln.ML - M2$ln.ML, M4$ln.ML - M3$ln.ML))

df$category <- cut(df$outcome, breaks=c(-Inf, 0, 2, 6, 10, Inf), 
                   labels=c("H_1 is better","Not worth more than a mention","Positive", "Strong", "Very Strong"))

ggplot(df, aes(B10, outcome)) +
  geom_col(aes(fill = category)) +
  coord_cartesian(ylim = c(-1, 11))

M1ln <- M1$ln.ML
M2ln <- M2$ln.ML
M3ln <- M3$ln.ML
M4ln <- M4$ln.ML

Mln_vec <- c(M1ln, M2ln, M3ln, M4ln)
maxMln <- max(Mln_vec)

Mln_vec_stand <- 2*(maxMln - Mln_vec)

plot(Mln_vec_stand, ylim = c(0,150), xlab = "Region", ylab = "Posterior weight")
lines(Mln_vec_stand)
abline(h = 2, lty = 2, col = "red")
abline(h = 6, lty = 2, col = "blue")
abline(h = 10, lty = 2, col = "green")


BF_matrix <- matrix(NA, nrow = 4, ncol = 4)
diag(BF_matrix) <- 1
BF_matrix[2,1] <- BF_12r
BF_matrix[3,1] <- BF_13r
BF_matrix[4,1] <- BF_14r
BF_matrix[3,2] <- BF_23r
BF_matrix[4,2] <- BF_24r
BF_matrix[4,3] <- BF_34r

library(xtable)

xtable(BF_matrix, display = c("s","e","e","e","e"), digits = 4)
