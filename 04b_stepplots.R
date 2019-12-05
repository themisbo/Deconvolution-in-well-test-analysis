
RD_dat <- dat %>% dplyr::select(., starts_with("R"))
M_dat <- dat %>% dplyr::select(., starts_with("M"))
eta_dat <- dat %>% dplyr::select(., starts_with("eta"))

MAP_point <- MAP(out)$parametersMAP

if(Region == 1){
  RD1_test <- as.matrix(RD_dat)
}else{
  RD1_test <- as.matrix(t(apply(exp(RD_dat), 1, cumsum)))
  RD1_test <- as.matrix( log(RD1_test) )
}

if(Region == 1){
  M1_test <- as.matrix(M_dat)
}else{
  M1_test <- as.matrix(t(apply(exp(M_dat), 1, cumprod)))
  M1_test <- as.matrix( log(M1_test) )
}

if(Region == 1){
  eta1_test <- as.matrix(eta_dat)
}else{
  eta1_test <- as.matrix(eta_dat)
}

RDi_list <- list(minRDi <- apply(RD1_test, 2, min),
                 mnRDi <- apply(RD1_test, 2, mean),
                 maxRDi <- apply(RD1_test, 2, max),
                 RD1_test,
                 RDi95 <- apply(RD1_test, 2, function(x){quantile(x,probs=c(.025,.975))}),
                 RDi99 <- apply(RD1_test, 2, function(x){quantile(x,probs=c(.005,.995))}),
                 MAPRDi <- cumsum( exp( MAP_point[(1:Region) + 3] )  ))

M1i_list <- list(minM1i <- apply(M1_test, 2, min),
                 mnM1i <- apply(M1_test, 2, mean),
                 maxM1i <- apply(M1_test, 2, max),
                 M1_test,
                 M1i95 <- apply(M1_test, 2, function(x){quantile(x,probs=c(.025,.975))}),
                 M1i99 <- apply(M1_test, 2, function(x){quantile(x,probs=c(.005,.995))}),
                 MAPM1i <- cumprod(exp( MAP_point[2*Region + 1:Region + 3]  )))

eta1i_list <- list(mineta1i <- apply(eta1_test, 2, min),
                   mneta1i <- apply(eta1_test, 2, mean),
                   maxeta1i <- apply(eta1_test, 2, max),
                   eta1_test,
                   eta1i95 <- apply(eta1_test, 2, function(x){quantile(x, probs=c(.025,.975))}),
                   eta1i99 <- apply(eta1_test, 2, function(x){quantile(x, probs=c(.005,.995))}),
                   MAPeta1i <- exp( MAP_point[Region + 1:Region + 3] ))

Uncty_plot <- function(xx, yy, type, zoom, plotminmax, plot99, plot95, ylab_string){
  
  minx <- unlist(xx[[1]])
  mnx  <- unlist(xx[[2]])
  maxx <- unlist(xx[[3]])
  x    <- unlist(xx[[4]])
  x95  <- unlist(xx[[5]])
  x99  <- unlist(xx[[6]])
  MAPx <- unlist(xx[[7]])
  
  if (Region==1) x95 = cbind(x95, x95)
  if (Region==1) x99 = cbind(x99, x99)
  
  miny <- unlist(yy[[1]])
  mny  <- unlist(yy[[2]])
  maxy <- unlist(yy[[3]])
  y    <- unlist(yy[[4]])
  y95  <- unlist(yy[[5]])
  y99  <- unlist(yy[[6]])
  MAPy <- unlist(yy[[7]])
  
  if (Region==1) y95 = cbind(y95, y95)
  if (Region==1) y99 = cbind(y99, y99)

  plot(c(0, mnx, maxx[Region], maxx[Region]+10), c(0, mny, mny[Region], mny[Region]),
       ty = "s", col = "gold4", lwd = 2, xlim=c(0, maxx[Region]+2),
       ylim=c(min(min(miny) - 1,-1), max(maxy, 2)) , xlab = "RD", ylab = ylab_string,cex.lab = 1.5, cex.axis = 2)
  
  if(plotminmax == TRUE){
    rect(minx, c(0,miny[1:(Region-1)]), maxx, maxy, col = "cyan4", border = NA)
    rect(minx, miny, c(maxx[2:Region], maxx[Region]), maxy, col = "cyan4", border = NA)
  }
  
  if(plot99 == TRUE){
    rect(x99[1,], c(0, y99[1,1:(Region-1)]), x99[2,], y99[2,], col = "cyan3", border = NA)
    rect(x99[1,], y99[1,], c(x99[2,2:Region], maxx[Region]), y99[2,], col = "cyan3", border = NA)
  }
  
  if(plot95 == TRUE){
    rect(x95[1,], c(0, y95[1,1:(Region-1)]), x95[2,], y95[2,], col = "cyan", border = NA)
    rect(x95[1,], y95[1,], c(x95[2,2:Region], maxx[Region]), y95[2,], col = "cyan", border = NA)
  }
  
  lines(c(min(RD1_test, 0), mnx, maxx[Region]), c(0, mny, mny[Region]), ty = "s", col = "gold4", lwd = 3)
  lines(c(min(RD1_test, 0), mnx, maxx[Region]), c(0, mny, mny[Region]), ty = "s", col = "gold3", lty = 2, lwd = 3)
  points(mnx, mny, col = 2:(Region+1), pch = 19)
  
}
#X11()
Uncty_plot(RDi_list, M1i_list, zoom = "mnx",
           plotminmax = T, plot99 = T, plot95 = T, ylab_string = "M" )

Uncty_plot(RDi_list, eta1i_list, zoom = "mnx",
           plotminmax = T, plot99 = T, plot95 = T,ylab_string = expression(eta))
