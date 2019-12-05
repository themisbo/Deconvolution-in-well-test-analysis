
library(BayesianTools)

Model <- function(parm) {
  
  sigmap  <- parm[3*Region + 4]
  
  TM     <- 10^parm[1]
  PM        <- 10^parm[2]
  CDe2Skin  <- 10^parm[3]
  
  RD        <- cumsum( 10^( parm[(1:Region) + 3] )  )
  eta       <- 10^( parm[Region + 1:Region + 3]  )
  M         <- 10^( parm[2*Region + 1:Region + 3] )
  
  tD = TM*t
  
  dpwD = function(s){
    
    sqrt_z <- sqrt(s/CDe2Skin)
    
    a11 = s*besselI(sqrt_z, 0, expon.scaled = TRUE) - sqrt_z*besselI(sqrt_z, 1, expon.scaled = TRUE)
    a12 = s*besselK(sqrt_z, 0, expon.scaled = TRUE) + sqrt_z*besselK(sqrt_z, 1, expon.scaled = TRUE)
    
    bb = array(dim = c(2*Region, 2*Region, dim(s)) )
    
    RD1_s       <- RD[1]*sqrt(s)
    RDN_etaN_s  <- ( RD[Region]*sqrt(eta[Region]) )*sqrt(s)
    
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
    
    cc = bb
    cc[1,1,,] <-      besselI(RD1_s, 0, expon.scaled = TRUE) # A - 20
    cc[2,1,,] <- M[1]*besselI(RD1_s, 1, expon.scaled = TRUE) # A - 22
    
    scale_exp = exp(2*sqrt_z - 2*RD1_s)*( detBTM(bb,s)/detBTM(cc,s) )
    Bscale = scale_exp*a11 - a12
    Ascale = Bscale/scale_exp
    
    return( besselI(sqrt_z, 0, expon.scaled = TRUE)/Ascale
            -besselK(sqrt_z, 0, expon.scaled = TRUE)/Bscale )
  }
  
  dpwD_time <- stehfest(dpwD, tD)
  
  log_tg <- log(tD*dpwD_time/PM)
  
  if(anyNA(log_tg)) log_tg = rep(1, length(log_tg))
  
  resp            <- Response(log_tg, tau, 1, 0, tN)
  
  C  <- Themis.calc.CMatrix(resp, I, c(m, N) )
  
  A1 = m/(sigmap^2) + 1/(sigmap0^2)
  A_1  = 1/A1
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
  
  return(llz)
}

density = function(par){
  d1a = dnorm(par[1], 2, .2, log =TRUE)
  d1b = dnorm(par[2], 1.5, .2, log =TRUE)
  d1c = dgamma(par[3], 1, .2, log =TRUE)
  d2 = dnorm(par[4], 2, 1, log =TRUE) + 
    if(Region > 1) sum(dnorm(par[(2:Region) + 3], 2, 1, log =TRUE)) else 0
    d3 = sum(dnorm(par[(Region + 4):(3*Region + 3)], 0, 1, log =TRUE))
  d4 = dunif(par[3*Region + 4], 0, 5, log =TRUE)
  return(d1a + d1b + d1c + d2 + d3 + d4)
}

sampler = function(n=1){
  d1a = matrix(rnorm(n, 2, .2), ncol = 1, byrow = T)
  d1b = matrix(rnorm(n, 1.5, .2), ncol = 1, byrow = T)
  d1c = matrix(rgamma(n, 1, .2), ncol = 1, byrow = T)
  d2 = matrix(rnorm(n, 2, 1), ncol = 1, byrow = T)
  if(Region > 1) d2 = cbind(d2, matrix(rnorm((Region-1)*n, 2, 1), ncol = Region-1, byrow = T) )
    d3 = matrix(rnorm(2*Region*n, 0, 1), ncol = 2*Region, byrow = T)
  d4 = matrix(runif(n, 0, 5), ncol = 1, byrow = T)
  return(cbind(d1a, d1b, d1c, d2, d3, d4))
}

low = c(rep(-3, 3), rep(-100, Region), rep(-10, 2*Region), 1e-8 )
up = c(10, 10, 60, rep(100, Region), rep(10, 2*Region), 5 )

prior <- createPrior(density = density, sampler = sampler, lower = low, upper = up, best = NULL)

bayesianSetup <- createBayesianSetup(likelihood = Model, prior = prior, parallel = FALSE)

sv = rbind(sampler(2), Fit$Posterior1[which.max(Fit$Monitor),])

settings <- list(iterations = 5e6,  message = TRUE, thin = 50, startValue = sv)

out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DEzs", settings = settings)

save(out, file = paste0("oil", Region, "r_DEzspriorlog10.RData"))
