
library(MASS)

sigmaq = .05*max(q)
sigmap0 = 10

p0q_samples <- function(parm) {

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
  
  # Convolution
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
  
  p0q = as.numeric( mvrnorm(1, mu_cond, Sigma_cond) )
  
  return(p0q)
}

p0q_res <- t(apply(Posterior, 1, p0q_samples))

save(p0q_res, file = paste0("oil", Region, "r_DEzsRATES.RData"))
