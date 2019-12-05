
 Posterior_modified <- getSample(out, coda = F, start = start_from)


if(Region == 1){
  Posterior_modified[, (1:Region) + 3] <- exp(Posterior_modified[, (1:Region) + 3])
  Posterior_modified[, 2*Region + 1:Region + 3] <- exp(Posterior_modified[, 2*Region + 1:Region + 3])
}else{
  Posterior_modified[, (1:Region) + 3] <- t(apply(exp(Posterior_modified[, (1:Region) + 3]), 1, cumsum))
  Posterior_modified[, 2*Region + 1:Region + 3] <- t(apply(exp(Posterior_modified[, 2*Region + 1:Region + 3]), 1, cumprod))
}

Posterior_modified[, Region + 1:Region + 3] <- exp( Posterior_modified[, Region + 1:Region + 3])

all_cov = array(dim = c(dim(Posterior_modified)[1]-1, dim(Posterior_modified)[2], dim(Posterior_modified)[2]))
for(i in 1:(dim(Posterior_modified)[1]-1)){
  all_cov[i,,] = cov(Posterior_modified[1:(i+1),])
}

cors_modified <- cov2cor(all_cov[dim(Posterior_modified)[1]-1,,])

rownames(cors_modified) <- colnames(cors_modified) <- c("TM", "PM", "CDe2Skin", paste0("RD", 1:Region), paste0("eta", 1:Region), 
                                                        paste0("M", 1:Region), "sigma_p"  )

heatmap(cors_modified, symm = TRUE, keep.dendro = TRUE)
