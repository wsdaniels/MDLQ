
run.mdlq.mcmc <- function(y, X,
                          a.vec = rep(1, ncol(X)), b.vec = rep(1, ncol(X)), # Hyper priors for theta_i's
                          c.vec = rep(1, ncol(X)), d.vec = rep(1, ncol(X)), # Hyper priors for tau2_i's
                          alpha1 = 1, alpha2 = 1, # Hyper priors on b
                          p = 0.001,
                          n.samples, n.burn.in = round(n.samples / 4, 0),
                          plot.trace = F) {
  
  
  k <- ncol(X)
  n <- nrow(X)
  
  # res is where we store the posterior samples
  res <- matrix(NA, nrow = n.samples, ncol = 4*k + 1)
  
  colnames(res) <- c(paste0('pi', seq(k)),
                     paste0('s', seq(k)),
                     paste0('theta', seq(k)),
                     paste0('beta', seq(k)),
                     'b')
  
  # take the MLE estimate as the values for the first sample
  m <- lm(y ~ X - 1)
  sigma.orig <- var(predict(m) - y)
  beta.orig <- ifelse(coef(m) > 0, coef(m), 0.25)
  beta.orig <- ifelse(beta.orig < 100, beta.orig, 10)
  res[1, ] <- c(rep(0, k), rep(1, k), rep(0.5, k), beta.orig, sigma.orig)
  res[1, ] <- ifelse(is.na(res[1,]), 0.25, res[1,])
  
  # we start running the Gibbs sampler
  for (i in seq(2, n.samples)) {
    
    # first, get all the values of the previous time point
    pi.prev <- res[i-1, seq(1, k)]
    s.prev <- res[i-1,seq(k + 1, 2*k)]
    theta.prev <- res[i-1, seq(2*k + 1, 3*k)]
    beta.prev <- res[i-1, seq(3*k + 1, 4*k)]
    b.prev <- res[i-1, ncol(res)]
    
    ## Start sampling from the conditional posterior distributions
    ##############################################################
    
    theta.new <- vector(length = k)
    for (j in sample(seq(k))){
      theta.new[j] <- rbeta(1, a.vec[j] + pi.prev[j], 1 - pi.prev[j] + b.vec[j])
    }
    
    b.target <- function(b){
      part1 <- n*log(1/(2*b)) + (alpha1+1) * log(1/b) 
      part2 <- -sum( abs( y - X %*% beta.prev ) ) / b
      part3 <- -alpha2 / b
      return(part1 + part2 + part3)
    }
    
    counter <- 1
    while(T){
      proposed.b <- rnorm(1, mean = b.prev, sd = 1)
      counter <- counter + 1
      if (proposed.b > 0){ break }
      if (counter > 10000){
        return( NA )
      }
    }
    
    accept.prob <- exp(b.target(proposed.b) - b.target(b.prev))
    if(runif(1) <= accept.prob) {
      b.new <- proposed.b
    } else {
      b.new <- b.prev
    }
    
    for (j in sample(seq(k))){
      if (pi.prev[j] == 1){
        s.target <- function(s.vec){
          part1 <- log(pi.prev[j] * (1/s.vec[j]))
          part2 <- -beta.prev[j]/s.vec[j]
          part3 <- log( (1/s.vec[j])^(c.vec[j]+1) )
          part4 <- -d.vec[j]/s.vec[j]
          return(part1 + part2 + part3 + part4)
        }
      } else {
        s.target <- function(s.vec){
          part1 <- log( (1-pi.prev[j]) * (1/p) )
          part2 <- -beta.prev[j]/p
          part3 <- log( (1/s.vec[j])^(c.vec[j]+1) )
          part4 <- -d.vec[j]/s.vec[j]
          return(part1 + part2 + part3 + part4)
        }
      }
      
      counter <- 1
      while(T){
        proposed.s <- rnorm(1, mean = s.prev[j], sd = 2)
        counter <- counter + 1
        if (proposed.s > 0){ break }
        if (counter > 10000){
          return(NA)
        }
      }
      
      proposed.s.vec <- s.prev
      proposed.s.vec[j] <- proposed.s
      
      accept.prob <- exp(s.target(proposed.s.vec) - s.target(s.prev))
      if(runif(1) <= accept.prob) {
        s.prev[j] <- proposed.s
      } else {
        s.prev[j] <- s.prev[j]
      }
    }
    
    s.new <- s.prev
    for (j in sample(seq(k))){
      if (pi.prev[j] == 1){
        beta.target <- function(beta.vec){
          part1 <- n * log( 1/(2*b.new) )
          part2 <- -sum( abs( y - X %*% beta.vec ) ) / b.new
          part3 <- (s.new[j]^-1)
          part4 <- -beta.vec[j] / s.new[j]
          return(part1 + part2 + part3 + part4)
        }
      } else {
        beta.target <- function(beta.vec){
          part1 <- n * log( 1/(2*b.new) )
          part2 <- -sum( abs( y - X %*% beta.vec ) ) / b.new
          part3 <- (p^-1)
          part4 <- -beta.vec[j] / p
          return(part1 + part2 + part3 + part4)
        }
      }
      
      counter <- 1
      while(T){
        proposed.beta <- rnorm(1, mean = beta.prev[j], sd = 0.1)
        counter <- counter + 1
        if (proposed.beta > 0){ break }
        if (counter > 10000){
          return(NA)
        }
      }
      
      proposed.beta.vec <- beta.prev
      proposed.beta.vec[j] <- proposed.beta
      
      accept.prob <- exp(beta.target(proposed.beta.vec) - beta.target(beta.prev))
      if(runif(1) <= accept.prob) {
        beta.prev[j] <- proposed.beta
      } else {
        beta.prev[j] <- beta.prev[j]
      }
    }
    
    beta.new <- beta.prev
    for (j in sample(seq(k))) {
      pi.target <- function(pi.vec){
        p0 <- (p^-1) * (1-theta.new[j]) * exp(-beta.new[j] / p)
        p1 <- (s.new[j]^-1) * theta.new[j] * exp(-beta.new[j] / s.new[j])
        psi <- exp(log(p0) - log(p0 + p1))
        return((1-psi)^pi.vec[j] * (psi)^(1-pi.vec[j]))
      }
      
      proposed.pi <- rbinom(1,1,0.5)
      proposed.pi.vec <- pi.prev
      proposed.pi.vec[j] <- proposed.pi
      
      accept.prob <- pi.target(proposed.pi.vec) / pi.target(pi.prev)
      if(runif(1) <= accept.prob) {
        pi.prev[j] <- proposed.pi
      } else {
        pi.prev[j] <- pi.prev[j]
      }
    }
    
    pi.new <- pi.prev
    
    res[i, ] <- c(pi.new, s.new, theta.new, beta.new, b.new)
    
  } # End Gibbs sampler
  
  out <- as.data.frame(res[-seq(n.burn.in), ])
  
  if (F){
    
    par(mfrow = c(ceiling(ncol(out) / floor(ncol(out)/3)),
                  floor(ncol(out)/3)))
    par(mar = c(2,2,2,2))
    for (i in 1:ncol(out)){
      plot(out[,i], type = "l", main = colnames(out)[i])
    }
    par(mfrow = c(1,1))
  }
  
  return(out)
  
}
