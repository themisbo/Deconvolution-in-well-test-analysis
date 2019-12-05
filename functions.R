##### Myfun all #######
# 1. Themis.create.Superposition
# 2. Themis.SuperposItem
# 3. Themis.calc.CMatrix

# Laplace Numerical inversion
# Stehfest parameters
NN=12

VV_fun <- function(NN){
  VV=0
  qr=0
  
  for(i in 1:NN){
    ss <- 0
    qq <- 0
    qq <- seq(floor((i+1)/2), min(i,NN/2), by=1)
    for(k in qq){
      ss[k] = sum(k^(NN/2)*factorial(2*k)/(factorial(NN/2-k)*factorial(k)*factorial(k-1)*factorial(i-k)*factorial(2*k-i)))
      qr[i] <- sum(ss, na.rm = TRUE)
    }
    VV[i] = ((-1)^(NN/2+i))*qr[i]
  }
  return(VV)
}

VV <- VV_fun(NN)

# Stehfest function
stehfest <- function(g,t){
  step1 <- g(outer(log(2)/t, 1:NN, "*"))
  step2 <- sweep(step1, MARGIN=2, VV, "*")
  step3 <- rowSums(step2)
  step4 <- log(2)/t*step3
  return(step4)
}


prod_ind <- function(tau, NN){
  ind_total <- list()
  grd <- expand.grid(1:length(tau), 1:NN)
  
  ind11 <- cbind(2*(1:(Region-1))-1, 2*(1:(Region-1)))   # (A - 26)
  ind11 <- cbind(ind11[rep(seq_len(nrow(ind11)), each=length(tau)*NN), ] , grd)
  ind_total[[1]] <- as.matrix(ind11)
  
  ind22 <- cbind(2*(1:(Region-1))-1, 2*(1:(Region-1))+1) # (A - 27)
  ind22 <- cbind(ind22[rep(seq_len(nrow(ind22)), each=length(tau)*NN), ] , grd)
  ind_total[[2]] <- as.matrix(ind22)
  
  ind33 <- cbind(2*(1:(Region-1)), 2*(1:(Region-1)))     # (A - 30)
  ind33 <- cbind(ind33[rep(seq_len(nrow(ind33)), each=length(tau)*NN), ] , grd)
  ind_total[[3]] <- as.matrix(ind33)
  
  ind44 <- cbind(2*(1:(Region-1)), 2*(1:(Region-1))+1)   # (A - 31)
  ind44 <- cbind(ind44[rep(seq_len(nrow(ind44)), each=length(tau)*NN), ] , grd)
  ind_total[[4]] <- as.matrix(ind44)
  
  ind55 <- cbind(2*(2:Region)-1, 2*(2:Region)-2)         # (A - 24)
  ind55 <- cbind(ind55[rep(seq_len(nrow(ind55)), each=length(tau)*NN), ] , grd)
  ind_total[[5]] <- as.matrix(ind55)
  
  ind66 <- cbind(2*(2:Region)-1, 2*(2:Region)-1)         # (A - 25)
  ind66 <- cbind(ind66[rep(seq_len(nrow(ind66)), each=length(tau)*NN), ] , grd)
  ind_total[[6]] <- as.matrix(ind66)
  
  ind77 <- cbind(2*(2:Region), 2*(2:Region)-2)           # (A - 28)
  ind77 <- cbind(ind77[rep(seq_len(nrow(ind77)), each=length(tau)*NN), ] , grd)
  ind_total[[7]] <- as.matrix(ind77)
  
  ind88 <- cbind(2*(2:Region), 2*(2:Region)-1)           # (A - 29)
  ind88 <- cbind(ind88[rep(seq_len(nrow(ind88)), each=length(tau)*NN), ] , grd)
  ind_total[[8]] <- as.matrix(ind88)
  
  return(ind_total)
}

if(Region != 1)  ind_total <- prod_ind(tau, NN)

detBTM <- function(bb,s){
  
  b_det <- array(0, c(2*Region, dim(s)))
  b_det_prime <- array(0, c(2*Region, dim(s)))
  
  b_det[1,,] <- bb[1,1,,]
  b_det[2,,] <- bb[1,1,,]*bb[2, 2,,] - bb[1,2,,]*bb[2,1,,]
  b_det_prime[2,,] <- bb[2,1,,]
  
  if(Region != 1){
    for ( i in 3:(2*Region) ){
      if (i %% 2 == 0){
        b_det_prime[i,,] <- bb[i,i-1,,]*b_det[i-2,,] - bb[i,i-2,,]*b_det_prime[i-1,,]
        b_det[i,,] <- bb[i,i,,]*b_det[i-1,,] - bb[i-1,i,,]*b_det_prime[i,,]
      }else{
        b_det_prime[i,,] <- bb[i-1,i,,]*b_det[i-2,,] - bb[i-2,i,,]*b_det_prime[i-1,,]
        b_det[i,,] <- bb[i,i,,]*b_det[i-1,,] - bb[i,i-1,,]*b_det_prime[i,,]
      }
    }
  }
  
  return(b_det[2*Region,,])
}

#### 2 functions ommited due to licencing

I <- Themis.create.Superposition(obsPress, obsRate[which_dd], Response(rep(1, length(tau)), tau, 1, 0, tN))[[1]]
