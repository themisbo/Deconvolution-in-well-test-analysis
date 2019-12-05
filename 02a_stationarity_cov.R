
library(coda)
library(factoextra)

Posterior_coda <- getSample(out, start = start_from, coda = TRUE) 

post_coda <- Posterior_coda[[1]][,c(1:4,Region +4, 2*Region+4)]
colnames(post_coda) <- c("TM", "PM", "CDe2Skin", paste0("RD", 1), paste0(sprintf('\u03B7'), 1), 
                         paste0("M", 1) )

geweke.plot(post_coda)

heidel.diag(post_coda)

effectiveSize(Posterior_coda)

crosscorr.plot(post_coda)

cumuplot(post_coda)

all_cov = array(dim = c(dim(Posterior)[1]-1, dim(Posterior)[2], dim(Posterior)[2]))
for(i in 1:(dim(Posterior)[1]-1)){
  all_cov[i,,] = cov(Posterior[1:(i+1),])
}

plot(1:(dim(Posterior)[1]-1), all_cov[,1,1], ty = "l", xlab = "Iteration", ylab = "Cov(1,1)")

cors <- cov2cor(all_cov[dim(Posterior)[1]-1,,])

rownames(cors) <- colnames(cors) <- c("TM", "PM", "CDe2Skin", paste0("RD", 1:Region), paste0("eta", 1:Region), 
                                      paste0("M", 1:Region), "sigma_p" )

pca <- prcomp(cors)

fviz_eig(pca)
fviz_pca_var(pca)
