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

## construct the overlap tensor for the given single-well arguments
Themis.create.Superposition <- function(press, rate, response){
  ## Construct the superposition overlap tensor I as defined in Appx. A
  ## I stores information needed to form the integrals which generate the
  ## coefficients C_ij of the convolution matrix (16). This takes the form of
  ## for every pressure point i, the range of rates j which have a non-zero
  ## contribution to the integral, and then for every (i,j) it determined which
  ## node intervals [tau_k-1, tau_k] have a non-zero contribution and what sub-interval
  ## of that node interval actually contributes to C_ij. It is this contributing
  ## sub-interval that is stored, along with the indices (i,j,k).
  ##
  ## The interval is stored either as its endpoints (cMidpointRadius=FALSE) or
  ## as the interval midpoint and half-width (cMidpointRadius=TRUE).
  ##
  ## Given the values of z and q, it is then relatively straightforward to
  ## compute C.
  ##
  ## NB: All problems with the same node positions, pressure obs times and flow period
  ## start-points have the same I, allowing for efficient construction of C as we
  ## change z and q during fitting.
  
  ##===============================================================================
  ## object sizes
  m <- length(press)
  nq <- length(rate)
  n <- length(response)
  
  ##===============================================================================
  ## Screen valid rate:pressure combinations
  s <- exp(response$Node) # nodes in raw time
  nodes1 <- c(0, s[-n]) # starting nodes in node intervals (0 is implicit)
  nodes2 <- s # terminating nodes in node intervals
  firstNode <- s[1]
  lastNode <- s[n]
  ## Determine the range of valid rates for each pressure_i, such that
  ## the vector (1:validRates[i]) gives the indices of rates prior to pressure_i.
  ## This is useful as rates occurring AFTER pressure_i cannot influence its value,
  ## and so can be completely ignored for any convolution for pressure_i.
  validRates <- rep(0, len=m)
  flowStart <- rate$Start
  flowEnd <- rate$End
  pTime <- press$Time
  #browser()
  validRates <- findInterval(pTime, flowStart)
  nijs <- sum(validRates) # num. valid combinations requiring consideration
  ##.debugPrintTime()
  
  ##===============================================================================
  ## Calculating overlap integral tensor values...
  I <- NULL
  if(m < nq){
    ## If we have fewer pressures than rates, then this loop is quickest
    Il <- vector(length=m, mode="list")
    sz <- 0
    js <- ti <- D1s <- D2s <- As <- Bs <- Ks <- 0
    for(i in 1:m){
      njs <- validRates[i]
      if(njs==0)
        next
      js <- 1:njs
      ## Get pressure time i, and the difference with appropriate rate times
      ti <- pTime[i]
      D1s <- ti - flowEnd[js]
      D2s <- ti - flowStart[js]
      ## [A,B] are now intervals which *may* have a contribution to the convolution
      ## integrals. A should be < B, so if B<A we have an empty interval
      As <- matrix(nodes1, ncol=njs, nrow=n)
      As <- sweep(As, 2, D1s, FUN=pmax)
      Bs <- matrix(nodes2, ncol=njs, nrow=n)
      Bs <- sweep(Bs, 2, D2s, FUN=pmin)
      ## Ks[k,j] informs us (TRUE) if we need to consider node k in conjunction
      ## with rate j (and pressure i), or if we can ignore it (FALSE)
      Ks <- (As < Bs)
      if(!any(Ks)){
        ## unique observation at the start of a flow period; ignore
        if(sum(As == Bs) == 1)
          next
        .debugWarning("No valid Ks.")
        if(config(".DEBUG")){
          browser()
          browser()
        }
      }
      As <- log(As[Ks])
      Bs <- log(Bs[Ks])
      js2 <- matrix(js, ncol=njs, nrow=n, byrow=TRUE)[Ks] # valid j indices in Ks[k,j]
      ks <- matrix(1:n, ncol=njs, nrow=n)[Ks] # valid k indices in Ks[k,j]
      nvals <- length(As)
      ## We now have all the (rate_j:node_k) combinations with a non-zero
      ## contribution to pressure_i, these are the non-empty elements of I_ijk.
      ## For each of these we need information on the associated interval.
      vals <- cbind(i, js2, ks, As, Bs)
      ## Found all of the intervals, stick in the list
      Il[[i]] <- vals
      sz <- sz + nvals
    }
    ## trash big temporary objects
    rm(As, Bs, Ks, D1s, D2s, vals)
    ## Populate the tensor object...
    I <- matrix(0, nrow=sz, ncol=5)
    pos <- 1
    for(i in 1:length(Il)){
      if(length(Il[[i]]) == 0)
        next
      nr <- nrow(Il[[i]])
      I[pos + 0:(nr - 1),] <- Il[[i]]
      pos <- pos + nr
    }
    rm(Il)
  } else {
    ##===============================================================================
    ## If we have fewer rates than pressures (most common), then this loop is quickest
    ## JAC: This is my own modification which deviates from the standard approach
    Il <- vector(mode="list", length=nq)
    sz <- 0
    js <- ti <- D1s <- D2s <- As <- Bs <- Ks <- 0
    for(j in 1:nq){
      validPress <- which(pTime >= flowStart[j])
      nps <- length(validPress)
      if(nps==0)
        next
      D1s <- pTime[validPress] - flowEnd[j]
      D2s <- pTime[validPress] - flowStart[j]
      As <- matrix(nodes1, ncol=nps, nrow=n)
      As <- sweep(As, 2, D1s, FUN=pmax)
      Bs <- matrix(nodes2, ncol=nps, nrow=n)
      Bs <- sweep(Bs, 2, D2s, FUN=pmin)
      Ks <- (As < Bs)
      As <- log(As[Ks])
      Bs <- log(Bs[Ks])
      is <- matrix(validPress, ncol=nps, nrow=n, byrow=TRUE)[Ks]
      ks <- matrix(1:n, ncol=nps, nrow=n)[Ks]
      nvals <- length(As)
      ## We now have all the (rate_j:node_k) combinations with a non-zero
      ## contribution to pressure_i, these are the non-empty elements of I_ijk.
      ## For each of these we need information on the associated interval.
      vals <- cbind(is, j, ks, As, Bs)
      Il[[j]] <- vals
      sz <- sz + nvals
    }
    ## trash big temporary objects
    rm(As, Bs, Ks, D1s, D2s, vals)
    ## Populate the tensor object...
    I <- matrix(0, nrow=sz, ncol=5)
    pos <- 1
    for(i in 1:length(Il)){
      if(length(Il[[i]]) == 0)
        next
      nr <- nrow(Il[[i]])
      I[pos + 0:(nr - 1),] <- Il[[i]]
      pos <- pos + nr
    }
    o <- order((max(I[, 2] + 1) * max(I[, 3] + 1)) * I[, 1] +
                 (max(I[, 3] + 1)) * I[, 2] + I[, 3])
    I <- I[o,,drop=FALSE]
    rm(Il)
  }
  colnames(I) <- c("i", "j", "k", "u", "v")
  ## create utility data structures
  whichInf <- which(I[,"u"]==-Inf)
  whichKs <- lapply(1:n, FUN=function(a) which(I[,"k"] == a))
  whichJs <- lapply(1:nq, FUN=function(a) which(I[,"j"] == a))
  ## create object
  I <- Themis.SuperposItem(I, m, nq, n, whichInf, whichKs, whichJs)
  return(I)
}

Themis.SuperposItem <- function(I, NP, NQ, NZ, whichInf, whichNodes, whichRates){
  stopifnot(ncol(I)==5)
  if(!inherits(I,"Matrix"))
    I <- Matrix(I)
  if(is.null(colnames(I)))
    colnames(I) <- c("i", "j", "k", "u", "v")
  stopifnot(length(whichNodes) <= NZ,
            length(whichRates) <= NQ)
  stopifnot(all(sapply(whichNodes, is.numeric)),
            all(sapply(whichRates, is.numeric)))
  return(new("SuperposItem",
             Core=I,
             Value=4,
             Pos=matrix(1, 1, 1),
             Len=matrix(nrow(I), 1, 1),
             Count=1,
             Symmetric=FALSE,
             NPress=NP, NRate=NQ, NNode=NZ,
             WhichInf=whichInf,
             WhichNodes=whichNodes,
             WhichRates=whichRates))
}

##===============================================================================
## Find the C matrix corresponding to the (i,j)th response function. 'dims' gives
## the dimensions of C which can't be inferred from I or Z as it is #press x #rate
Themis.calc.CMatrix <- function(response, I, dims){

  beta <- response$Slopes
  alpha <- response$Intercepts

  alpha <- exp(alpha) # alpha is only ever used when exponentiated

  if(any(!is.finite(alpha))){
    .debugWarning("Non-finite values of e^alpha.")
  }
  
  nrowI = nrow(I)
  ci <- numeric(nrowI)
  k <- I$k
  alphak <- alpha[k]
  betak <- beta[k]
  ## Calculate the contribution of rate j to well i due to the response function
  ## between nodes (k-1,k). This is then summed over k to obtain Cij
  ## JAC: Now uses an endpoint representation method (faster)
  ## I now contains the columns u, v which represent the interval [u,v] for which
  ## the integral component (i,j,k) is non-zero. u can be -Inf for contributions
  ## from before the first node.
  u <- I$u
  v <- I$v
  #case1 <- I@WhichInf # := which(I$u==-Inf); semi-infinite case
  case1 <- which(u==-Inf)
  if(length(case1) > 0){
    ak <- alphak[case1]
    bk <- betak[case1]
    ci[case1] <- ak * exp(v[case1] * bk) / bk
  }
  ## CASE 2
  ## flat case: 2 exp{alpha_k} delta_k
  case2 <- which(betak == 0) # beta=0 case
  if(length(case2) > 0){
    ak <- alphak[case2]
    ci[case2] <- ak * (v[case2] - u[case2])
  }
  ## finite, not-flat:
  ##  2 exp{alpha_k} exp{mu_k beta_k} sinh(beta_k delta_k) / beta_k
  done <- c(case1,case2)
  case3 <- (1:nrowI)
  if(length(done)>0)
    case3 <- case3[-done] # otherwise
  if(length(case3) > 0){
    ak <- alphak[case3]
    bk <- betak[case3]
    ci[case3] <- ak * (exp(v[case3] * bk) - exp(u[case3] * bk)) / bk
  }
  
  if(any(is.infinite(ci))){
    .debugWarning("Non-finite contributions to C matrix.")
  }
  
  if(length(case1)==0 && length(case2)==0 && length(case3)==0){
    .debugWarning("C matrix case failure")
  }
  
  Cmat <- sparseMatrix(i=I$i, j=I$j, x=ci, dims=dims, giveCsparse=FALSE)
  return(Cmat)
}

I <- Themis.create.Superposition(obsPress, obsRate[which_dd], Response(rep(1, length(tau)), tau, 1, 0, tN))[[1]]
