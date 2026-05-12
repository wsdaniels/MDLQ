
run.mdlq.mcmc <- function(y, X, # response vector and covariate matrix
                          a.vec = rep(1, ncol(X)), b.vec = rep(1, ncol(X)), # hyper priors for theta_i's
                          c.vec = rep(1, ncol(X)), d.vec = rep(1, ncol(X)), # hyper priors for tau2_i's
                          alpha1 = 1, alpha2 = 1, # hyper priors on nu
                          n.samples, n.burn.in = round(n.samples / 4, 0), # number of iterations to run the sampler and amount to burn in
                          plot.trace = F) { # flag to plot traces
  
  library(MASS)
  library(tmvtnorm)
  library(TruncatedNormal)
  library(R.utils)
  
  # Function to force symmetry and add a tiny ridge until chol succeeds (if needed)
  make.spd <- function(Sig, start.ridge = 0, max.ridge = 1e-4){
    d <- nrow(Sig)
    Sig <- 0.5 * (Sig + t(Sig))
    ridge <- start.ridge
    while (TRUE) {
      Sig.try <- if (ridge > 0) Sig + diag(ridge, d) else Sig
      ok <- try(chol(Sig.try), silent = TRUE)
      if (!inherits(ok, "try-error")) return(Sig.try)
      ridge <- if (ridge == 0) 1e-10 else min(ridge * 10, max.ridge)
      if (ridge >= max.ridge && inherits(ok, "try-error")) return(Sig.try) # last resort
    }
  }
  
  # Function to call TruncatedNormal::mvrandn and do error handling
  safe.mvrandn <- function(n, mu, Sig, lower, upper,
                           per.attempt.timeout.sec = 2, 
                           max_attempts = 5){
    d <- length(mu)
    stopifnot(length(lower) == d, length(upper) == d)
    Sig <- make.spd(Sig)
    
    attempt <- 1
    ridge   <- 0
    last.err <- NULL
    
    repeat {
      # strengthen ridge slightly across attempts, in case of borderline SPD
      Sig.use <- if (ridge > 0) Sig + diag(ridge, d) else Sig
      
      res <- try(
        withTimeout(
          expr = TruncatedNormal::mvrandn(n = n, l = lower, u = upper, mu = mu, Sig = Sig.use),
          timeout = per.attempt.timeout.sec,
          onTimeout = "error"  
        ),
        silent = TRUE
      )
      
      if (!inherits(res, "try-error") && is.matrix(res)) {
        if (nrow(res) == d && ncol(res) == n && all(is.finite(res))) return(res)
        last.err <- "nonfinite_or_bad_shape"
      } else {
        last.err <- attr(res, "condition")
      }
      
      attempt <- attempt + 1
      if (attempt > max_attempts) {
        stop(sprintf("mvrandn timeout or error after %d attempts (%s)", max_attempts, as.character(last.err)))
      }
      # increase ridge a bit for numerical stability before retrying
      ridge <- if (ridge == 0) 1e-10 else ridge * 10
    }
  }
  
  # Get dimensions
  p <- ncol(X)
  n <- nrow(X)
  
  # res is where we store the posterior chains
  res <- matrix(NA, nrow = n.samples, ncol = 4*p + 3)
  colnames(res) <- c(paste0('z', seq(p)),
                     paste0('tau2.', seq(p)),
                     paste0('theta', seq(p)),
                     paste0('beta', seq(p)),
                     'sigma2', 'r', "nu")
  
  # take the MLE estimate as the values for the first sample
  m <- lm(y ~ X - 1)
  sigma.orig <- var(predict(m) - y)
  beta.orig <- ifelse(coef(m) > 0, coef(m), 0.25)
  beta.orig <- ifelse(beta.orig < 100, beta.orig, 10)
  
  # Initialize res as before
  res[1, ] <- c(rep(0, p), rep(1, p), rep(0.5, p), beta.orig, sigma.orig, 0, 2)
  res[1, ] <- ifelse(is.na(res[1,]), 0.25, res[1,])
  
  # compute only once
  XtX <- t(X) %*% X
  Xty <- t(X) %*% y
  
  # Run the Gibbs sampler
  for (i in seq(2, n.samples)) {
    
    print(paste0(i, "/", n.samples))
    
    # Get the parameter values from the previous iteration
    z.prev <- res[i-1, seq(1, p)]
    tau2.prev <- res[i-1,seq(p + 1, 2*p)]
    theta.prev <- res[i-1, seq(2*p + 1, 3*p)]
    beta.prev <- res[i-1, seq(3*p + 1, 4*p)]
    sigma2.prev <- res[i-1, ncol(res)-2]
    r.prev <- res[i-1, ncol(res)-1]
    nu.prev <- res[i-1, ncol(res)]
    
    
    #------------------------------------------------------------
    #--- Sample the probability of an emission: theta
    #------------------------------------------------------------
    theta.new <- vector(length = p)
    for (j in sample(seq(p))){
      theta.new[j] <- rbeta(1, a.vec[j] + z.prev[j], 1 - z.prev[j] + b.vec[j])
    }
    
    # Construct the correlation matrix for the errors
    R.coef <- 1 / (1-r.prev^2)
    R.mat <- matrix(0, nrow = n, ncol = n)
    diag(R.mat) <- 1 + r.prev^2
    R.mat[1,1] <- R.mat[n,n] <- 1
    R.mat[row(R.mat) == (col(R.mat)-1)] <- -r.prev
    R.mat[row(R.mat) == (col(R.mat)+1)] <- -r.prev
    
    # Compute inverse of R
    R.inv <- R.coef * R.mat
    
    # Compute residuals based on previous rate estimates
    err <- y - X %*% beta.prev
    
    #------------------------------------------------------------
    #--- Sample the error variance: sigma2
    #------------------------------------------------------------
    sigma2.new <- 1 / rgamma(1, n/2 + nu.prev/2, t(err) %*% R.inv %*% err / 2 + nu.prev/2)
    
    
    #------------------------------------------------------------
    #--- Sample the prior belief about degrees of freedom: nu
    #------------------------------------------------------------
    nu.target <- function(nu){
      part1 <- (nu/2) * log(nu/2)
      part2 <- -lgamma(nu/2)
      part3 <- (-(nu/2)-1) * log(sigma2.new)
      part4 <- (-alpha1-1) * log(nu)
      part5 <- -(nu/2) / sigma2.new
      part6 <- -alpha2 / nu
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
    } else {
      nu.new <- nu.prev
    }
    
    
    #------------------------------------------------------------
    #--- Sample the autocorrelation coefficient: r
    #------------------------------------------------------------
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
    } else {
      r.new <- r.prev
    }
    
    
    #------------------------------------------------------------
    #--- Sample the emission rate scale parameter: tau2
    #------------------------------------------------------------
    tau2.new <- vector(length = p)
    for (j in sample(seq(p))){
      tau2.new[j] <- 1 / rgamma(1,
                                c.vec[j] + z.prev[j],
                                d.vec[j] + beta.prev[j])
    }
    
    
    # Compute once
    XtRX <- t(X) %*% R.inv %*% X
    XtRy <- t(X) %*% R.inv %*% y    
    
    #------------------------------------------------------------
    #--- Sample the emission rates: beta
    #------------------------------------------------------------
    beta.cov <- qr.solve( (1/sigma2.new) * XtRX )
    beta.cov <- as.matrix(nearPD(beta.cov)$mat)
    beta.new <- vector(length = p)
    
    for (j in sample(seq(p))){
      e <- rep(0, p)
      e[j] <- 1
      
      # Construct the mean of the update distribution
      beta.mean <- as.vector(beta.cov %*% ((XtRy/sigma2.new) - matrix(e/(tau2.new), nrow = p)))
      
      # bounds for the orthant constraint beta >= 0
      lower.vec <- rep(0, length(beta.mean))
      upper.vec <- rep(Inf, length(beta.mean))
      
      # Try to get 100 samples; if it times out or errors, throw "broke on betas" error
      samples.dxN <- NULL
      try({
        samples.dxN <- safe.mvrandn(
          n = 100,
          mu = beta.mean,
          Sig = beta.cov,
          lower = lower.vec,
          upper = upper.vec,
          per.attempt.timeout.sec = 60,  # hard cap per attempt
          max_attempts = 5               # a few retries with adaptive ridge
        )
      }, silent = TRUE)
      
      if (is.null(samples.dxN)) {
        print("broke on betas")
        return("broke on betas")
      }
      
      # Convert to (N x d) 
      beta.samples <- t(samples.dxN)                      
      beta.samples <- matrix(beta.samples, ncol = ncol(X))  
      
      # Trim out first few samples
      beta.samples.trimmed <- beta.samples[50:100, , drop = FALSE]
      beta.samples.trimmed[is.infinite(beta.samples.trimmed)] <- NA
      
      # Choose the last fully finite row; if none, treat as failure
      is.row.ok <- apply(beta.samples.trimmed, 1, function(z) all(is.finite(z)))
      if (!any(is.row.ok)) {
        print("broke on betas")
        return("broke on betas")
      }
      last.row <- max(which(is.row.ok))
      
      beta.new[j] <- beta.samples.trimmed[last.row, j]
      if (is.na(beta.new[j])) {
        print("broke on betas")
        return("broke on betas")
      }
      
    } # end loop through the betas
    
    
    #------------------------------------------------------------
    #--- Sample the spike-slab indicator: z
    #------------------------------------------------------------
    for (j in sample(seq(p))) {
      
      # get the betas for which beta_j is zero
      z0 <- z.prev
      z0[j] <- 0
      bz0 <- matrix(beta.new * z0, nrow = p)
      
      # compute the w term
      xj <- X[, j]
      xj.star <- R.inv %*% X[,j]
      sum.x2 <- sum(xj*xj.star)
      w <- y - X %*% bz0
      w.star <- R.inv %*% y - R.inv %*% X %*% bz0
      
      # compute chance parameter of the conditional posterior of z_j (Bernoulli)
      # Break it up into the z=0 and z=1 cases
      l0 <- log(1 - theta.new[j])
      l1 <- log(theta.new[j] / (2 * tau2.new[j])) + 
        (sum(xj.star * w + xj * w.star) - 2*sigma2.new/tau2.new[j])^2 / (4 * sigma2.new * sum.x2) + 
        0.5 * log( 2 * pi * sigma2.new / sum.x2) 
      
      # sample z_j from a Bernoulli
      z.prev[j] <- rbinom(1, 1, 1 - (exp(l0) / (exp(l0) + exp(l1))))
    }
    
    z.new <- z.prev
    
    # save the new parameter values from this iteration of the sampler
    res[i, ] <- c(z.new, tau2.new, theta.new, beta.new*z.new, sigma2.new, r.new, nu.new)
  } # End Gibbs sampler
  
  # Remove the burn in samples
  out <- as.data.frame(res[-seq(n.burn.in), ])
  
  if (plot.trace){
    par(mfrow = c(ceiling(ncol(out) / floor(ncol(out)/3)), 
                  floor(ncol(out)/3)))
    par(mar = c(2,2,2,2))
    
    for (i in 1:ncol(out)){
      plot(out[,i], type = "l", main = colnames(out)[i])
    }
  }
  
  return(out)
}


