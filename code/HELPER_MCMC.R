
inv.transform <- function(sim.order, p, beta.mean, beta.chol, j){
  # print(paste0( i, "/", n.samples, " - sim order: ", sim.order))
  num.samples <- round(ifelse(25^log(sim.order) < 500, 500, 25^log(sim.order)))
  these.betas <- vector(length= num.samples)
  for (o in 1:length(these.betas)){
    u <- runif(p, min = 1-10^(-sim.order), max = 1)
    Z <- qnorm(u)
    these.betas[o] <- (beta.mean + beta.chol %*% Z)[j]
  }
  these.betas[is.infinite(these.betas)] <- NA
  # hist(these.betas, xlim = c(-0.1, 0.01))
  # hist(these.betas, xlim = c(-0.06,0.01), breaks = seq(-0.06, 0.02, by = 0.01/2))
  # hist(these.betas, xlim = c(-0.06,0.01), breaks = seq(-0.06, 0.02, by = 0.01/2), col = alpha("green", 0.1), add = T)
  this.beta <- NA
  if (any(these.betas > 0, na.rm = T)){
    this.beta <- these.betas[these.betas > 0 & !is.na(these.betas)][1]
    return(this.beta)
  } else {
    return(NA)
  }
  
}

mills.approx <- function(b, p, beta.mean, beta.chol, j){
  # print(paste0( i, "/", n.samples, " - b value: ", b))
  num.samples <- 500
  these.betas <- vector(length= num.samples)
  for (o in 1:length(these.betas)){
    U <- runif(p)
    Z <- sqrt(b^2 - 2 * log(U/b))
    these.betas[o] <- (beta.mean + beta.chol %*% Z)[j]
  }
  
  if (any(these.betas > 0)){
    this.beta <- these.betas[these.betas > 0][1]
    return(list(this.beta, mean(these.betas)))
  } else {
    return(list(NA, mean(these.betas)))
  }
  
}


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
    
    beta.cov <- qr.solve( (1/sigma2.new) * XtRX )
    beta.chol <- chol(beta.cov)
    beta.new <- vector(length = p)
    
    for (j in sample(seq(p))){
      e <- rep(0, p)
      e[j] <- 1
      
      beta.mean <- beta.cov %*% ((XtRy/sigma2.new) - matrix(e/(tau2.new * sigma2.new), nrow = p))
      
      found.beta <- F
      
      # Check the first inv transform
      this.beta <- inv.transform(sim.order = 0, p, beta.mean, beta.chol, j)
      if (!is.na(this.beta)){ found.beta <- T }
      
      # If first inv transform didn't work, do other stuff
      if (!found.beta){
        
        # First check the last inv transform
        this.beta <- inv.transform(sim.order = 10, p, beta.mean, beta.chol, j)
        if (!is.na(this.beta)){ found.beta <- T }
        
        # If the last inv transform found something, check the intermediate inv transforms
        if (found.beta){
          for (repeater in 1:10){
            for (sim.order in seq(0,10, by = 0.2)){
              this.beta <- inv.transform(sim.order, p, beta.mean, beta.chol, j)
              if (!is.na(this.beta)){ break }
            }
            if (!is.na(this.beta)){ break }
          }
          if (!is.na(this.beta)){
            beta.new[j] <- this.beta  
          } else {
            found.beta <- F
          }
        }
        
        # If the last inv transform didn't find anything, move onto the mills approx
        if (!found.beta){
          b.seq <- seq(1, 10000, by = 100)
          avg.betas <- vector(length = length(b.seq))
          count <- 1
          for (b in b.seq){
            tmp <- mills.approx(b, p, beta.mean, beta.chol, j)
            this.beta <- tmp[[1]]
            avg.betas[count] <- tmp[[2]]
            count <- count + 1
            if (!is.na(this.beta)){ 
              found.beta <- T
              break 
            }
          }
          
          # If the mills approx found something, zoom in the b grid
          if (found.beta){
            
            if (b == 1){
              b.seq.zoom <- seq(1,100, by = 1)
            } else {
              # NOTE: the extra plus 2 is for when a beta was found at the very upper
              # end of the range. A resample from that range might not find it again,
              # so extend the range slightly so that a beta can still be found
              b.seq.zoom <- seq(b.seq[count-2], b.seq[count-1+2], by = 1)
            }
            
            avg.betas <- vector(length = length(b.seq.zoom))
            for (b in b.seq.zoom){
              tmp <- mills.approx(b, p, beta.mean, beta.chol, j)
              this.beta <- tmp[[1]]
              if (!is.na(this.beta)){ 
                break 
              }
            }
            
            beta.new[j] <- this.beta  
            
            # If the mills approx didnt' find something, diagnose why not
          } else {
            
            this.fit <- lm(avg.betas ~ b.seq)
            if (coef(this.fit)[2] < 0){
              print('decreasing betas')
              return('decreasing betas')
            } else {
              print('b did not go high enough')
              return('b did not go high enough')
            }
            
          } # End check if the mills approx found something
        } # End check if the last inv transform found something
        
        # If first inv transform worked, save beta
      } else {
        beta.new[j] <- this.beta  
      }
      
    } # end loop through the betas
    
    
    
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
    
    par(mfrow = c(ceiling(ncol(out) / floor(ncol(out)/3)), 
                  floor(ncol(out)/3)))
    par(mar = c(2,2,2,2))
    
    for (i in 1:ncol(out)){
      plot(out[,i], type = "l", main = colnames(out)[i])
    }
    
  }
  
  # remove the first n.burnin number of samples
  return(out)
}


