auto.corr <- readRDS('./results/simulation_study_output_autocorrelation.RData')

results <- vector(mode = "list", length = length(auto.corr))

for (zz in 1:length(auto.corr)){
  
  print(paste0(zz, "/", length(auto.corr)))
  
  out <- auto.corr[[zz]]
  
  q.hat <- q.hat.lower <- q.hat.upper <- matrix(NA, nrow = length(out), ncol = 5)
  sigma2 <- r <- rep(NA, length(out))
  
  for (i in 1:length(out)){
    
    these.mcmc <- out[[i]]
    
    to.use <- sapply(these.mcmc, length) > 1
    these.mcmc <- these.mcmc[to.use]
    
    if (length(these.mcmc) == 0){ next }
    
    all.betas <- do.call(rbind, lapply(these.mcmc, function(X) apply(X, 2, mean)[c(16:20)]))
    q.hat[i, ] <- apply(matrix(all.betas, ncol = 5), 2, mean) * 3.6
    
    all.lower <- do.call(rbind, lapply(these.mcmc, function(X) apply(X, 2, function(X) quantile(X, probs = 0.025))[c(16:20)]))
    q.hat.lower[i, ] <- apply(matrix(all.lower, ncol = 5), 2, mean) * 3.6
    
    all.upper <- do.call(rbind, lapply(these.mcmc, function(X) apply(X, 2, function(X) quantile(X, probs = 0.975))[c(16:20)]))
    q.hat.upper[i, ] <- apply(matrix(all.upper, ncol = 5), 2, mean) * 3.6
    
    all.sigma2 <- do.call(rbind, lapply(these.mcmc, function(X) apply(X, 2, mean)[c(21)]))
    sigma2[i] <- apply(all.sigma2, 2, mean) 
    
    all.r <- do.call(rbind, lapply(these.mcmc, function(X) apply(X, 2, mean)[c(22)]))
    r[i] <- apply(all.r, 2, mean) 
    
  }
  
  results[[zz]] <- list(q.hat = q.hat,
                        q.hat.lower = q.hat.lower,
                        q.hat.upper = q.hat.upper,
                        sigma2 = sigma2,
                        r = r)
  
}

saveRDS(results, './results/simulation_study_output_autocorrelation_condensed.RData')

rm(auto.corr)
gc()

spike.misalignment <- readRDS('./results/simulation_study_output_spike_misalignment.RData')

results <- vector(mode = "list", length = length(spike.misalignment))

for (zz in 1:length(spike.misalignment)){
  
  print(paste0(zz, "/", length(spike.misalignment)))
  
  out <- spike.misalignment[[zz]]
  
  q.hat <- q.hat.lower <- q.hat.upper <- matrix(NA, nrow = length(out), ncol = 5)
  sigma2 <- r <- rep(NA, length(out))
  
  for (i in 1:length(out)){
    
    these.mcmc <- out[[i]]
    
    to.use <- sapply(these.mcmc, length) > 1
    these.mcmc <- these.mcmc[to.use]
    
    if (length(these.mcmc) == 0){ next }
    
    all.betas <- do.call(rbind, lapply(these.mcmc, function(X) apply(X, 2, mean)[c(16:20)]))
    q.hat[i, ] <- apply(matrix(all.betas, ncol = 5), 2, mean) * 3.6
    
    all.lower <- do.call(rbind, lapply(these.mcmc, function(X) apply(X, 2, function(X) quantile(X, probs = 0.025))[c(16:20)]))
    q.hat.lower[i, ] <- apply(matrix(all.lower, ncol = 5), 2, mean) * 3.6
    
    all.upper <- do.call(rbind, lapply(these.mcmc, function(X) apply(X, 2, function(X) quantile(X, probs = 0.975))[c(16:20)]))
    q.hat.upper[i, ] <- apply(matrix(all.upper, ncol = 5), 2, mean) * 3.6
    
    all.sigma2 <- do.call(rbind, lapply(these.mcmc, function(X) apply(X, 2, mean)[c(21)]))
    sigma2[i] <- apply(all.sigma2, 2, mean) 
    
    all.r <- do.call(rbind, lapply(these.mcmc, function(X) apply(X, 2, mean)[c(22)]))
    r[i] <- apply(all.r, 2, mean) 
    
  }
  
  results[[zz]] <- list(q.hat = q.hat,
                        q.hat.lower = q.hat.lower,
                        q.hat.upper = q.hat.upper,
                        sigma2 = sigma2,
                        r = r)
  
}

saveRDS(results, './results/simulation_study_output_spike_misalignment_condensed.RData')