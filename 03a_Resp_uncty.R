
if(Region == 1){
  uncty_col = "deepskyblue"
  MAP_col = "navyblue"
}else if(Region == 2){
  uncty_col = "darkolivegreen2"
  MAP_col = "darkolivegreen4"
}else if(Region == 3){
  uncty_col = "gold3"
  MAP_col = "gold4"
}else if(Region == 4){
  uncty_col = "violet"
  MAP_col = "violetred4"
}

 N_Nodes = 100
 t0 = 1e-3
 tN = ceiling( obsPress[nrow(obsPress)]$Time )
 tau = seq(log(t0), log(tN), length.out = N_Nodes)
 t = exp(tau)

 source("functions.R")

MAP_point <- MAP(out)$parametersMAP

df_Resp <- dat %>%
  dplyr::select(1:(3*(Region + 1))) %>%
  slice(seq(1, nrow(.), by = 20)) 

get_z <- function(parm){

  TM     <- 10^parm[1]
  PM        <- 10^parm[2]
  CDe2Skin  <- 10^parm[3]
  RD        <- cumsum(  10^( parm[(1:Region) + 3] ) )
  eta       <- 10^( parm[Region + 1:Region + 3] )
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
  dpwD_time <- stehfest(dpwD,tD)
  
  # Proposed Response
  return(log(tD*dpwD_time/PM))}

z <- apply(df_Resp, 1, get_z)

df_z <- as_tibble(t(z))
names(df_z) <- tau
df_z2 <- df_z %>%
  gather(tau, value, convert = TRUE)

df_data1 <- tibble(tau, MAP_z = get_z(MAP_point))
df_data2 <- tibble(tau = trueResp$Node, true_z = trueResp$Z)

df_z2$id <- 1:nrow(df_z)

ggplot(df_z2, aes(x = tau, y = value)) +
  geom_line(aes(group = id), col = uncty_col, alpha = 0.1) +
  geom_line(data=df_data1, aes(x = tau, y = MAP_z), col = MAP_col) +
  geom_line(data=df_data2, aes(x = tau, y = true_z), col = 'grey30')+
  theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size=20) ) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size=20) ) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  labs(y = "z", x = "log(t)")
