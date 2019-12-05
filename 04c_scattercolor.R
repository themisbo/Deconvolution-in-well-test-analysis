dat_full <- dat

if(Region == 1){
ggplot(dat_full, aes(eta1, M1)) +
  geom_point(aes(color = RD1)) + xlab(paste0(sprintf('\u03B7'), 1))
}else if(Region == 2){
  gg1 <- ggplot(dat_full, aes(eta1, M1)) +
    geom_point(aes(color = RD1)) + xlab(paste0(sprintf('\u03B7'), 1))
  gg2 <- ggplot(dat_full, aes(eta2, M2)) +
    geom_point(aes(color = RD2)) + xlab(paste0(sprintf('\u03B7'), 2))
  grid.arrange(gg1, gg2 ,ncol=1)
}else if(Region == 3){
  gg1 <- ggplot(dat_full, aes(eta1, M1)) +
    geom_point(aes(color = RD1)) + xlab(paste0(sprintf('\u03B7'), 1))
  gg2 <- ggplot(dat_full, aes(eta2, M2)) +
    geom_point(aes(color = RD2)) + xlab(paste0(sprintf('\u03B7'), 2))
  gg3 <- ggplot(dat_full, aes(eta3, M3)) +
    geom_point(aes(color = RD3)) + xlab(paste0(sprintf('\u03B7'), 3))
  grid.arrange(gg1, gg2, gg3 , layout_matrix = rbind(c(1,2),c(3,3)))
}else if(Region == 4){
  gg1 <- ggplot(dat_full, aes(eta1, M1)) +
    geom_point(aes(color = RD1)) + xlab(paste0(sprintf('\u03B7'), 1))
  gg2 <- ggplot(dat_full, aes(eta2, M2)) +
    geom_point(aes(color = RD2)) + xlab(paste0(sprintf('\u03B7'), 2))
  gg3 <- ggplot(dat_full, aes(eta3, M3)) +
    geom_point(aes(color = RD3)) + xlab(paste0(sprintf('\u03B7'), 3))
  gg4 <- ggplot(dat_full, aes(eta4, M4)) +
    geom_point(aes(color = RD4)) + xlab(paste0(sprintf('\u03B7'), 4))
  grid.arrange(gg1, gg2, gg3, gg4, ncol=2)
}


library(gg3D)
theta=0
phi=20

ggplot(dat_full, aes(x=eta1, y=RD1, z=M1, col = RD1) ) + 
  theme_void() +
  axes_3D(theta=theta, phi=phi) +
  stat_3D(theta=theta, phi=phi) +
  labs_3D(theta=0, phi=20, hjust=c(1,0,0), vjust=c(1.5,1,-.2),
          labs=c(paste0(sprintf('\u03B7'), 1), "RD1", "M1"))


