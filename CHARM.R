# Erase everything in memory
rm(list=ls())

library(mdeconv)
library(Matrix)
library(LaplacesDemon)

load("real_TLS.RData")

obsPress <- real_TLS$swPress
obsRate    <- real_TLS$swRate
trueResp    <- real_TLS$TlsResp

which_dd <- which(obsRate$Rate != 0)

set.seed(123)

# Initial pressure p0
p0 <- max(obsPress$Press)

# Number of regions (transitions)
Region = 1

N = length(which_dd)
m = length(obsPress$Press)
q = obsRate$Rate[which_dd]
p = obsPress$Press
N_Nodes = 30
t0 = 1e-3
tN = ceiling( obsPress[nrow(obsPress)]$Time )
tau = seq(log(t0), log(tN), length.out = N_Nodes)
t = exp(tau)

source("functions.R")

#sigmap = .5
sigmaq = .05*max(q)
sigmap0 = 10

Data <- list(N = length(which_dd), m =length(obsPress$Press),
             q = q, p = p,
             p0 = p0, tau = tau, tN = tN,
             mon.names = "LP",
             parm.names = c("TM", "PM", "CDe2Skin",
                            paste0("RD", 1:Region, "_sqrt(CD)"), paste0("eta", 1:Region), 
                            paste0("M", 1:Region), "sigma_p" ) )

Model <- function(parm, data) {
  
  parm[3*Region + 4] <- interval(parm[3*Region + 4], 0, 5)
  sigmap  <- parm[3*Region + 4]
  
  TM     <- 10^parm[1]
  PM        <- 10^parm[2]
  CDe2Skin  <- 10^parm[3]
  
  RD        <- cumsum( 10^( parm[(1:Region) + 3] )  )
  eta       <- 10^( parm[Region + 1:Region + 3]  )
  M         <- 10^( parm[2*Region + 1:Region + 3] )
  
  # Dimensionless time
  tD = TM*t
  
  # Dimensionless pressure derivative in Laplace space
  dpwD = function(s){
    
    # Wellbore region solution parameters
    sqrt_z <- sqrt(s/CDe2Skin)
    
    a11 = s*besselI(sqrt_z, 0, expon.scaled = TRUE) - sqrt_z*besselI(sqrt_z, 1, expon.scaled = TRUE)
    a12 = s*besselK(sqrt_z, 0, expon.scaled = TRUE) + sqrt_z*besselK(sqrt_z, 1, expon.scaled = TRUE)
    
    # a11 element sub-determinant
    bb = array(dim = c(2*Region, 2*Region, dim(s)) )
    
    RD1_s       <- RD[1]*sqrt(s)
    RDN_etaN_s  <- ( RD[Region]*sqrt(eta[Region]) )*sqrt(s)
    
    # Solution parameters
    sqeta  <- sqrt(eta[1:(Region-1)])
    Msqeta <- M[2:Region]*sqeta
    
    bb[1,1,,] <-       besselK(RD1_s, 0, expon.scaled = TRUE) # A - 21
    bb[2,1,,] <- -M[1]*besselK(RD1_s, 1, expon.scaled = TRUE) # A - 23
    
    if(Region != 1){
      
      RD1i_etai_s <- outer(sqrt(s), RD[1:(Region-1)]*sqeta, "*")
      RD2i_etai_s <- outer(sqrt(s), RD[2:Region]*sqeta, "*")
      
      bb[ind_total[[1]]] <-       -besselI(RD1i_etai_s, 0, expon.scaled = TRUE) # A - 26
      bb[ind_total[[2]]] <-       -besselK(RD1i_etai_s, 0, expon.scaled = TRUE) # A - 27
      bb[ind_total[[3]]] <-  sweep(besselI(RD1i_etai_s, 1, expon.scaled = TRUE), MARGIN = 3, -sqeta, "*") # A - 30
      bb[ind_total[[4]]] <-  sweep(besselK(RD1i_etai_s, 1, expon.scaled = TRUE), MARGIN = 3,  sqeta, "*") # A - 31
      
      bb[ind_total[[5]]] <-         besselI(RD2i_etai_s, 0, expon.scaled = TRUE) # A - 24
      bb[ind_total[[6]]] <-         besselK(RD2i_etai_s, 0, expon.scaled = TRUE) # A - 25
      bb[ind_total[[7]]] <-   sweep(besselI(RD2i_etai_s, 1, expon.scaled = TRUE), MARGIN = 3, Msqeta, "*") # A - 28
      bb[ind_total[[8]]] <-   sweep(besselK(RD2i_etai_s, 1, expon.scaled = TRUE), MARGIN = 3,-Msqeta, "*") # A - 29
    }
    
    bb[2*Region-1, 2*Region,,] <-                -besselK(RDN_etaN_s, 0, expon.scaled = TRUE) # A - 27
    bb[2*Region, 2*Region,,] <- sqrt(eta[Region])*besselK(RDN_etaN_s, 1, expon.scaled = TRUE) # A - 31
    
    # a12 element sub-determinant
    cc = bb
    cc[1,1,,] <-      besselI(RD1_s, 0, expon.scaled = TRUE) # A - 20
    cc[2,1,,] <- M[1]*besselI(RD1_s, 1, expon.scaled = TRUE) # A - 22
    
    scale_exp = exp(2*sqrt_z - 2*RD1_s)*( detBTM(bb,s)/detBTM(cc,s) )
    Bscale = scale_exp*a11 - a12
    Ascale = Bscale/scale_exp
    
    # Solution of the Dimensionless pressure derivative in Laplace space
    return( besselI(sqrt_z, 0, expon.scaled = TRUE)/Ascale
            -besselK(sqrt_z, 0, expon.scaled = TRUE)/Bscale )
  }
  
  # Dimensionless pressure derivative in time space
  dpwD_time <- stehfest(dpwD, tD)
  
  # Proposed Response
  log_tg <- log(tD*dpwD_time/PM)
  
  if(anyNA(log_tg)) log_tg = rep(1, length(log_tg))
  
  # Response object
  resp            <- Response(log_tg, tau, 1, 0, tN)
  
  # Convolution
  C  <- Themis.calc.CMatrix(resp, I, c(m, N) )
  
  A1 = m/(sigmap^2) + 1/(sigmap0^2)
  A_1  = 1/A1 # sigmap^2*sigmap0^2/(sigmap^2 + m*sigmap0^2)
  B1 = - as.matrix(rep(1/sigmap^2, m)%*%C)
  D1 = as.matrix(diag(1/sigmaq^2, N) + t(C)%*%C/sigmap^2)
  
  A_matrix = rbind(cbind(A1, B1), cbind(t(B1), D1))
  
  b_sq1 <- sum(p)/sigmap^2 + p0/(sigmap0^2)
  b_sq2 <- t(q/sigmaq^2) - t(p/sigmap^2) %*% C
  b_sq = cbind(b_sq1, b_sq2)
  
  Sig_22 = chol2inv(chol(D1-t(B1)%*%A_1%*%B1))
  Sig_12 = -A_1%*%B1%*%Sig_22
  Sig_21 = t(Sig_12)
  Sig_11 = A_1 - Sig_12 %*% t(B1) %*% A_1
  
  Sigma_cond <- rbind(cbind(Sig_11, Sig_12), cbind(Sig_21, Sig_22))
  mu_cond = Sigma_cond %*% t(b_sq)
  
  c1 = - m*log(sigmap) - N*log(sigmaq) - log(sigmap0) -m/2*log(2*pi) -N/2*log(2*pi) -log(2*pi)/2 
  c2 = -( sum(p^2)/sigmap^2 + sum(q^2)/sigmaq^2 + p0^2/sigmap0^2  )/2
  c3 = drop(t(mu_cond) %*% A_matrix %*% mu_cond)/2
  c4 = log( det( 2*pi*Sigma_cond) )/2
  
  llz = c1 + c2 + c3 + c4
  if(llz == Inf) llz <- -1e12
  
  d1a = dnorm(parm[1], 2, 0.2, log =TRUE)
  d1b = dnorm(parm[2], 1.5, 0.2, log =TRUE)
  d1c = dgamma(parm[3], 1, .2, log =TRUE)
  d2 = dnorm(parm[4], 2, 1, log =TRUE) + 
    if(Region > 1) sum(dnorm(parm[(2:Region) + 3], 2, 1, log =TRUE)) else 0
    d3 = sum(dnorm(parm[(Region + 4):(3*Region + 3)], 0, 1, log =TRUE))
  d4 = dunif(parm[3*Region + 4], 0, 5, log =TRUE)
  
  lld <- d1a + d1b + d1c + d2 + d3 + d4
  if(lld == -Inf) lld <- -1e12
  
  llp <- llz + lld
  
  Modelout <- list(LP = llp, Dev = -2 * llz, Monitor = llp, yhat = NA, parm = parm)
  return(Modelout)
}

Initial.Values      <-  c(3, 2.5, 5, rep(1.3, Region), rep(0.2, 2*Region), 3)

Fit <- LaplacesDemon(Model, Data=Data, Initial.Values, Covar=NULL,
                     Iterations=2e4, Status=100, Thinning=1,
                     Algorithm="CHARM", Specs=list(alpha.star=0.44),
                     Debug=list(DB.chol=FALSE, DB.eigen=FALSE, DB.MCSE=FALSE, DB.Model=FALSE))

save(Fit, file = paste0("oil", Region, "r_CHARMlog10.RData"))
