
obsPress <- real_TLS$swPress
obsRate    <- real_TLS$swRate
trueResp    <- real_TLS$TlsResp
which_dd <- which(obsRate$Rate != 0)
which_bu <- which(obsRate$Rate == 0)


tau_q <- c(obsRate$Start, last(obsRate$End))

N_Nodes = 30
t0 = 1e-3
tN = ceiling( obsPress[nrow(obsPress)]$Time )
tau = seq(log(t0), log(tN), length.out = N_Nodes)
t = exp(tau)

source("functions.R")

 dat_p0q <- as_tibble(p0q_res)

 names(dat_p0q) <- c("p0", paste0("q", 1:N) )

dat_full <- bind_cols(dat, dat_p0q)

p0q_samples_plot <- function(parm) {
  
  sigmap  <- parm[3*Region + 4]
  
  # Early time parameters
  TM     <- 10^parm[1]
  PM        <- 10^parm[2]
  CDe2Skin  <- 10^parm[3]
  
  # Late time parameters
  RD        <- cumsum( 10^( parm[(1:Region) + 3]))
  eta       <- 10^( parm[Region + 1:Region + 3])
  M         <- 10^( parm[2*Region + 1:Region + 3])
  
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
  
  return(list("p0q" = p0q, "mu_cond" = as.numeric(mu_cond)))
}

MAP_point[(3*Region + 5):(3*Region + N + 5)] <- p0q_samples_plot(MAP_point)$mu_cond

tau_dataRate <- c(obsRate$Start, last(obsRate$End))

df_data <- tibble(tau_qdata = tau_dataRate, q_data = c(obsRate$Rate, 0))#, q_true = trueRate$Rate)
df_MAP <- df_data
df_MAP$q_data[which_dd] <- MAP_point[(3*Region + 6):( 3*Region + 5 + N )]

df_Rate <- dat_full %>%
  dplyr::select(starts_with("q"))

for (i in which_bu){
  df_Rate <- df_Rate %>%
    add_column(!!(paste0('qb', i)) := 0, .before = i)
}

df_Rate <- df_Rate %>%
  add_column(qlast = last(.))
names(df_Rate) <- tau_q

df_Rate_thin <- df_Rate %>%
  slice(seq(1, nrow(.), by = 20))

df_Rate2 <- df_Rate_thin %>%
  gather(tau_q, value, convert = TRUE)

df_Rate2$id <- 1:nrow(df_Rate_thin)

ggplot(df_Rate2, aes(x = tau_q, y = value)) +
  geom_step(aes(group = id, col = "Samples"), alpha = 0.1) +
  geom_step(data=df_MAP, aes(x = tau_q, y = q_data, col = "MAP")) +
  geom_step(data=df_data, aes(x = tau_qdata, y = q_data, col = "Data"))+
  scale_color_manual(values = c(Samples = uncty_col, Data = 'red', MAP = MAP_col)) + 
  labs(x = "Time", y= "q")

ggplot(df_Rate2, aes(x = tau_q, y = value)) +
  geom_step(aes(group = id), col = uncty_col, alpha = 0.1) +
  geom_step(data=df_MAP, aes(x = tau_q, y = q_data), col = MAP_col) +
  geom_step(data=df_data, aes(x = tau_qdata, y = q_data), col = 'red')+
  theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size=20) ) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size=20) ) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  labs(x = "Time", y= "q")
