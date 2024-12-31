# a.vec = rep(1, ncol(X))
# b.vec = rep(1, ncol(X))
# c.vec = rep(1, ncol(X))
# d.vec = rep(1, ncol(X))
# k = l = 1
# n.samples = 1500
# n.burn.in = round(n.samples / 4, 0)

ss.regress <- function(y, X, 
                       a.vec = rep(1, ncol(X)), b.vec = rep(1, ncol(X)), # Hyper priors for theta_i's
                       c.vec = rep(1, ncol(X)), d.vec = rep(1, ncol(X)), # Hyper priors for tau2_i's
                       k = 1, l = 1, # hyper priors on nu
                       n.samples, n.burn.in = round(n.samples / 4, 0),
                       plot.trace = F) {
  
  library(MASS)
  
  p <- ncol(X)
  n <- nrow(X)
  
  # res is where we store the posterior samples
  res <- matrix(NA, nrow = n.samples, ncol = 4*p + 3)
  
  colnames(res) <- c(paste0('pi', seq(p)),
                     paste0('tau2.', seq(p)),
                     paste0('theta', seq(p)),
                     paste0('beta', seq(p)),
                     'sigma2', 'r', "nu")
  
  # take the MLE estimate as the values for the first sample
  m <- lm(y ~ X - 1)
  sigma.orig <- var(predict(m) - y)
  beta.orig <- ifelse(coef(m) > 0, coef(m), 0.25)
  beta.orig <- ifelse(beta.orig < 100, beta.orig, 10)
  res[1, ] <- c(rep(0, p), rep(1, p), rep(0.5, p), beta.orig, sigma.orig, 0, 2)
  res[1, ] <- ifelse(is.na(res[1,]), 0.25, res[1,])
  
  # compute only once
  XtX <- t(X) %*% X
  Xty <- t(X) %*% y
  
  accepted <- vector(length = n.samples)
  
  # we start running the Gibbs sampler
  for (i in seq(2, n.samples)) {
    
    # print(paste0( i, "/", n.samples))
    
    # first, get all the values of the previous time point
    pi.prev <- res[i-1, seq(1, p)]
    tau2.prev <- res[i-1,seq(p + 1, 2*p)]
    theta.prev <- res[i-1, seq(2*p + 1, 3*p)]
    beta.prev <- res[i-1, seq(3*p + 1, 4*p)]
    sigma2.prev <- res[i-1, ncol(res)-2]
    r.prev <- res[i-1, ncol(res)-1]
    nu.prev <- res[i-1, ncol(res)]
    
    ## Start sampling from the conditional posterior distributions
    ##############################################################
    
    # sample theta from a Beta
    theta.new <- vector(length = p)
    for (j in sample(seq(p))){
      theta.new[j] <- rbeta(1, a.vec[j] + pi.prev[j], 1 - pi.prev[j] + b.vec[j])
    }
    
    R.coef <- 1 / (1-r.prev^2)
    R.mat <- matrix(0, nrow = n, ncol = n)
    diag(R.mat) <- 1 + r.prev^2
    R.mat[1,1] <- R.mat[n,n] <- 1
    R.mat[row(R.mat) == (col(R.mat)-1)] <- -r.prev
    R.mat[row(R.mat) == (col(R.mat)+1)] <- -r.prev
    
    R.inv <- R.coef * R.mat
    
    err <- y - X %*% beta.prev
    
    # sample sigma2 from an Inverse-Gamma
    sigma2.new <- 1 / rgamma(1, n/2 + nu.prev/2, t(err) %*% R.inv %*% err / 2 + nu.prev/2)
    
    nu.target <- function(nu){
      
      part1 <- (nu/2) * log(nu/2)
      part2 <- -lgamma(nu/2)
      part3 <- (-(nu/2)-1) * log(sigma2.new)
      part4 <- (-k-1) * log(nu)
      part5 <- -(nu/2) / sigma2.new
      part6 <- -l / nu
      
      return(part1 + part2 + part3 + part4 + part5 + part6)
    }
    
    counter <- 1
    while(T){
      proposed.nu <- rnorm(1, mean = nu.prev, sd = 10)
      counter <- counter + 1
      if (proposed.nu > 0){ break }
      if (counter > 10000){
        print("broke on nu")
        return(NA)
      }
    }
    
    accept.prob <- exp(nu.target(proposed.nu) - nu.target(nu.prev))
    
    if(runif(1) <= accept.prob) {
      nu.new <- proposed.nu
      # accepted[i] <- T
    } else {
      nu.new <- nu.prev
      # accepted[i] <- F
    }
    
    
    #R <- qr.solve(R.inv)
    
    r.target <- function(r){
      
      R.coef <- 1 / (1-r^2)
      R.mat <- matrix(0, nrow = n, ncol = n)
      diag(R.mat) <- 1 + r^2
      R.mat[1,1] <- R.mat[n,n] <- 1
      R.mat[row(R.mat) == (col(R.mat)-1)] <- -r
      R.mat[row(R.mat) == (col(R.mat)+1)] <- -r
      
      R.inv <- R.coef * R.mat
      
      part1 <- log( (1-r^2)^(-(n-1)/2) )
      part2 <- (-1/(2*sigma2.new)) * t(err) %*% R.inv %*% err
      
      return(part1 + part2)
    }
    
    counter <- 1
    while(T){
      proposed.r <- rnorm(1, mean = r.prev, sd = 0.2)
      counter <- counter + 1
      if (proposed.r > 0 & proposed.r <= 0.97){ break }
      if (counter > 10000){
        print("broke on r")
        return(NA)
      }
    }
    
    accept.prob <- exp(r.target(proposed.r) - r.target(r.prev))
    
    if(runif(1) <= accept.prob) {
      r.new <- proposed.r
      accepted[i] <- T
    } else {
      r.new <- r.prev
      accepted[i] <- F
    }
    
    
    # sample tau2 from an Inverse Gamma
    tau2.new <- vector(length = p)
    for (j in sample(seq(p))){
      tau2.new[j] <- 1 / rgamma(1,
                                c.vec[j] + pi.prev[j],
                                d.vec[j] + beta.prev[j]/sigma2.new)
    }
    
    XtRX <- t(X) %*% R.inv %*% X
    XtRy <- t(X) %*% R.inv %*% y
    
    
    # sample beta from multivariate Gaussian
    # if(all(is.na(
    #   tryCatch({ forceSymmetric(qr.solve( (1/sigma2.new) * XtRX )) },
    #            error = function(msg){ return(NA) })
    # ))){
    #   print("broke on beta cov")
    #   return(NA)
    # } else {
    #   beta.cov <- forceSymmetric(qr.solve( (1/sigma2.new) * XtRX ))
    # }
    
    beta.cov <- qr.solve( (1/sigma2.new) * XtRX )
    
    beta.new <- vector(length = p)
    
    
    
    for (j in sample(seq(p))){
      e <- rep(0, p)
      e[j] <- 1
      
      beta.mean <- beta.cov %*% ((XtRy/sigma2.new) - matrix(e/(tau2.new * sigma2.new), nrow = p))
      
      
      counter <- 1
      while(T){
        this.beta <- mvrnorm(n = 1, mu = as.numeric(beta.mean), Sigma = beta.cov)[j]
        counter <- counter + 1
        if (this.beta >= 0) {break}
        if (counter > 10000) {
          print("broke on betas")
          return(NA)
        }
      }
      
      beta.new[j] <- this.beta
    }
    
    # sample each pi_j in random order
    for (j in sample(seq(p))) {
      
      # get the betas for which beta_j is zero
      pi0 <- pi.prev
      pi0[j] <- 0
      bp0 <- matrix(beta.new * pi0, nrow = p)
      
      # compute the z variables
      xj <- X[, j]
      xj.star <- R.inv %*% X[,j]
      sum.x2 <- sum(xj*xj.star)
      z <- y - X %*% bp0
      z.star <- R.inv %*% y - R.inv %*% X %*% bp0
      
      # compute chance parameter of the conditional posterior of pi_j (Bernoulli)
      l0 <- log(1 - theta.new[j])
      
      l1 <- log(theta.new[j] / (2 * tau2.new[j] * sigma2.new)) + 
        (sum(xj.star * z + xj * z.star) - 2/tau2.new[j])^2 / (4 * sigma2.new * sum.x2) + 
        0.5 * log( 2 * pi * sigma2.new / sum.x2) 
      
      
      # sample pi_j from a Bernoulli
      pi.prev[j] <- rbinom(1, 1, 1 - (exp(l0) / (exp(l0) + exp(l1))))
    }
    
    pi.new <- pi.prev
    
    # add new samples
    res[i, ] <- c(pi.new, tau2.new, theta.new, beta.new*pi.new, sigma2.new, r.new, nu.new)
    
  } # End Gibbs sampler
  
  out <- as.data.frame(res[-seq(n.burn.in), ])
  
  if (plot.trace){
    
    # png(save.dir, res = 100, pointsize = 24, width = 1920, height = 1080)
    
    par(mfrow = c(ceiling(ncol(out) / floor(ncol(out)/3)), 
                  floor(ncol(out)/3)))
    par(mar = c(2,2,2,2))
    
    for (i in 1:ncol(out)){
      plot(out[,i], type = "l", main = colnames(out)[i])
    }
    
    # dev.off()
  }
  
  # remove the first n.burnin number of samples
  return(out)
}


