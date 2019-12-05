load("oil1r_DEzs.RData")
#load("oil2r_DEzs.RData")
#load("oil3r_DEzs.RData")
#load("oil4r_DEzs.RData")

library(BayesianTools)
library(tidyverse)
library(gridExtra)
library(GGally)
library(coda)
library(xtable)

load("real_TLS.RData")
obsPress <- real_TLS$swPress
obsRate    <- real_TLS$swRate
trueResp    <- real_TLS$TlsResp
which_dd <- which(obsRate$Rate != 0)

p0_value <- max(obsPress$Press)

start_from = 30001
Posterior <- getSample(out, start = start_from, coda = FALSE)

colnames(Posterior) <- c("TM", "PM", "CDe2Skin", paste0("RD", 1:Region), paste0("eta", 1:Region), 
                         paste0("M", 1:Region), "sigma_p" )

colnames(Posterior) <- c("T", "P", "W", paste0("R", 1:Region), paste0("eta", 1:Region), 
                         paste0("M", 1:Region), "sigma_p" )

Region = (ncol(Posterior) - 4)/3

dat <- as_tibble(Posterior)

p1a <- dat %>%
  ggplot(., aes(TM)) + 
  geom_histogram()

p2a <- dat %>%
  ggplot(., aes(PM)) + 
  geom_histogram()

p3a <- dat %>%
  ggplot(., aes(CDe2Skin)) + 
  geom_histogram()

p4a <- dat %>%
  ggplot(., aes(RD1)) + 
  geom_histogram()

p5a <- dat %>%
  ggplot(., aes(eta1)) + 
  geom_histogram()+ xlab(paste0(sprintf('\u03B7'), 1))

p6a <- dat %>%
  ggplot(., aes(M1)) + 
  geom_histogram()

p1b <- dat %>%
  ggplot(., aes(x=as.numeric(rownames(.)), y = TM)) + 
  geom_line() + labs(x = "Iterations")

p2b <- dat %>%
  ggplot(., aes(x=as.numeric(rownames(.)), y = PM)) + 
  geom_line() + labs(x = "Iterations")

p3b <- dat %>%
  ggplot(., aes(x=as.numeric(rownames(.)), y = CDe2Skin)) + 
  geom_line() + labs(x = "Iterations")

p4b <- dat %>%
  ggplot(., aes(x=as.numeric(rownames(.)), y = RD1)) + 
  geom_line() + labs(x = "Iterations")

p5b <- dat %>%
  ggplot(., aes(x=as.numeric(rownames(.)), y = eta1)) + 
  geom_line() + labs(x = "Iterations") + ylab(paste0(sprintf('\u03B7'), 1))

p6b <- dat %>%
  ggplot(., aes(x=as.numeric(rownames(.)), y = M1)) + 
  geom_line() + labs(x = "Iterations")

grid.arrange(p1a, p2a, p3a, p4a, p5a, p6a, p1b, p2b, p3b, p4b, p5b, p6b,ncol=6)

autocorPosterior <- getSample(out, start = start_from, coda = TRUE, whichParameters = 1:6 )[[1]]

colnames(autocorPosterior) <- c("TM", "PM", "CDe2Skin", paste0("RD", 1), paste0(sprintf('\u03B7'), 1), 
                                paste0("M", 1) )
autocorr.plot(autocorPosterior, lag.max = 100)

dat %>%
  dplyr::select(starts_with("sigma")) %>%
  gather(variable, value) %>%
  ggplot(., aes(x = value)) + 
  facet_wrap(~ variable, scales = "free", ncol = 4) + 
  geom_histogram()

data_long <- gather(dat, factor_key=TRUE)

print(xtable(data_long%>% group_by(key)%>%
               summarise(min = min(value), "25%"=quantile(value, probs=0.25), mean= mean(value),
                         "75%"=quantile(value, probs=0.75), max = max(value), sd= sd(value))), include.rownames=FALSE)

ggpairs(dat)

load("oil1r_DEzsRATES.RData")
#load("oil2r_DEzsRATES.RData")
#load("oil3r_DEzsRATES.RData")
#load("oil4r_DEzsRATES.RData")

colnames(p0q_res) <- c("p0", paste0("q", 1:N) )

dat_p0q <- as_tibble(p0q_res)

dat_p0q %>%
  dplyr::select(starts_with("p0")) %>%
  gather(variable, value) %>%
  ggplot(., aes(x = value)) + 
  facet_wrap(~ variable, scales = "free", ncol = 4) + 
  geom_histogram() +
  geom_vline(xintercept = p0_value, col = "red")

dat_p0q %>%
  dplyr::select(starts_with("p0")) %>%
  gather(variable, value) %>%
  ggplot(., aes(x = value)) + 
  facet_wrap(~ variable, scales = "free", ncol = 4) + 
  geom_histogram() +
  theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_vline(xintercept = p0_value, col = "red")
