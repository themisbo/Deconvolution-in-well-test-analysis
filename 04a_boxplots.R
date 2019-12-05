
RD_dat <- dat %>% dplyr::select(., starts_with("RD"))

if(Region == 1){
  RD1_test <- RD_dat
  
  ggplot(RD1_test, aes(1, RD1)) +
    geom_boxplot() +
    labs(y = "value", x = "Transition")
  
}else{
  RD_dat2 <- t(apply(exp(RD_dat), 1, cumsum))
  RD_dat2 <- as_tibble( log(RD_dat2) )
  
  dataRD2_long <- gather(RD_dat2, factor_key=TRUE)
  
  ggplot(dataRD2_long, aes(key, value)) +
    geom_boxplot() +
    labs(y = "value", x = "Transition")
}

M_dat <- dat %>% dplyr::select(., starts_with("M"))

if(Region == 1){
  M1_test <- M_dat
  
  ggplot(M1_test, aes(1,M1)) +
    geom_boxplot() +
    labs(y = "value", x = "Transition")
  
}else{
  M_dat2 <- t(apply(exp(M_dat), 1, cumprod))
  M_dat2 <- as_tibble( log(M_dat2) )
  
  dataM2_long <- gather(M_dat2, factor_key=TRUE)
  
  ggplot(dataM2_long, aes(key, value)) +
    geom_boxplot() +
    labs(y = "value", x = "Transition")
}

eta_dat <- dat %>% dplyr::select(., starts_with("eta"))

if(Region == 1){
  eta1_test <- eta_dat
  
  ggplot(eta1_test, aes(1,eta1)) +
    geom_boxplot() +
    labs(y = "value", x = "Transition")
  
}else{
  eta_dat2 <- eta_dat
  
  dataeta2_long <- gather(eta_dat2, factor_key=TRUE)
  
  ggplot(dataeta2_long, aes(key, value)) +
    geom_boxplot() +
    labs(y = "value", x = "Transition")
}



