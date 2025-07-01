# Description: Performs emission event detection, localization, and quantification
#              using predictions from the Gaussian puff simulation
# Author: William Daniels (wdaniels@mines.edu)
# Last Updated: November 10, 2022

# Clear environment
if(!is.null(dev.list())){dev.off()}
rm(list = ls())

# Import necessary libraries
library(lubridate)
library(foreach)
library(doParallel)
library(zoo)


# START USER INPUT
#---------------------------------------------------------------------------

num.cores.to.use <- 8

forward.model.path         <- '../input_data/forward_model_output_ADED2024.RData'
output.file.path           <- './results/simulation_study_output_spike_misalignment.RData'
spike.detection.alg.path   <- 'https://raw.github.com/wsdaniels/DLQ/master/code/HELPER_spike_detection_algorithm.R'
spike.slab.regression.path <- '../code/HELPER_spike_slab_regression.R'
helper.function.path       <- '../code/HELPER_functions.R'
truth.path                 <- '../output_data/truth_ADED2024_30min.RData'

# END OF USER INPUT - NO MODIFICATION NECESSARY BELOW THIS POINT
#---------------------------------------------------------------------------


# Source helper files which contain helper functions 
source(spike.detection.alg.path)
source(spike.slab.regression.path)
source(helper.function.path)

# Read in simulation data
data <- readRDS(forward.model.path)

# Read in true emission data
truth <- readRDS(truth.path)

# Pull out sensor observations and replace NA's that are not on edge of 
# the time series with interpolated values
obs <- na.approx(data$obs, na.rm = F)

# Number of sensors
n.r <- ncol(obs)

# Pull out time stamps of observations and simulations
times <- data$times

# Pull out the simulation predictions
sims <- data[5:length(data)]

# Grab source info
n.s <- length(sims) 
source.names <- names(sims)




# Quantification interval length [minutes]
interval.length <- 30

# Number of minutes to advance the quantification interval 
step.size <- 10

# Compute total number of intervals
num.intervals <- (length(times)/step.size) - (interval.length/step.size) + 1
num.intervals <- floor(num.intervals)


# True error variance to use in simulated observations
sigma2.true <- 1

# True spike misalignment percents
misalignment.vals <- c(0, 12.5, 25, 37.5, 50)/100

# Initialize variables to hold results
exp.results <- vector(mode = "list", length = length(misalignment.vals))
names(exp.results) <- misalignment.vals



# STEP 4: LOCALIZATION AND QUANTIFICATION 
# ---------------------------------------------------------------------------

for (zz in 1:length(misalignment.vals)){
  
  print(paste0(zz, "/", length(misalignment.vals)))
  
  # Fire up the parallel cluster
  cl <- makeCluster(num.cores.to.use)
  registerDoParallel(cl)
  
  big.out <- foreach(a = 1:num.intervals) %dopar% {
    
    q.hat <- q.hat.lower <- q.hat.upper <- sigma2 <- r.vec <- vector(length = n.s)
    
    print(paste0(a, "/", num.intervals))
    
    # Mask in this quantification interval
    sub.mask <- seq((a-1)*step.size + 1,
                    (a-1)*step.size + interval.length)
    
    # Number of time steps in this quantification interval
    N <- length(sub.mask)
    
    # Create matrix that holds dispersion model output
    X <- matrix(NA, nrow = N*n.r, ncol = n.s)
    colnames(X) <- source.names
    
    # Fill X matrix with simulation output
    for (s in 1:n.s){
      this.source <- c()
      for (r in 1:n.r){
        this.source <- c(this.source, sims[[s]][sub.mask, r])
      }
      X[,s] <- this.source
    }
    
    # Get observation data
    y <- c()
    for (r in 1:n.r){
      y <- c(y, obs[sub.mask, r])
    }
    
    betas <- as.matrix(truth[a,])
    
    y.no.error <- X %*% betas
    y <- y.no.error + rnorm(nrow(X), 0, sqrt(sigma2.true))
    
    # plot(y, type = "l")
    
    y.spikes <- find.spikes(times= 1:length(y.no.error), obs = y.no.error, 
                            amp.threshold = 0.5, make.plot = F)
    event.nums <- na.omit(unique(y.spikes$events))
    n.spikes <- length(event.nums)
    
    n.spikes.to.move <- round(n.spikes * misalignment.vals[zz])
    spike.ind.to.move <- sample.int(n.spikes, size = n.spikes.to.move, replace = F)
    
    if (n.spikes.to.move > 0){
      for (k in spike.ind.to.move){
        
        this.spike.mask <- which(y.spikes$events == event.nums[k])
        y[this.spike.mask] <- rnorm(length(this.spike.mask), 0, sqrt(sigma2.true))
      }
    }
    
    # lines(y, col = "red")
    
    if (n.spikes.to.move > 0){
      
      for (k in spike.ind.to.move){
        
        this.spike.mask <- y.spikes$events == event.nums[k]
        this.y.to.add <- na.omit(y.no.error[this.spike.mask])
        this.start.ind <- round(runif(1, min = 1, max = length(y) - length(this.y.to.add) - 1))
        to.add.mask <- seq(this.start.ind, this.start.ind + length(this.y.to.add) -1)
        y[to.add.mask] <- y[to.add.mask] + this.y.to.add
      }
      
    }
    
    # lines(y, col = "blue")
    
    info.mask <- apply(X, 2, function(this.col) sum(this.col > 0.5) > 4)
    
    q.hat[!info.mask] <- 
      q.hat.lower[!info.mask] <- q.hat.upper[!info.mask] <-
      sigma2[!info.mask] <- r.vec[!info.mask] <- "no info"
    
    if (all(y == 0)){
      q.hat[info.mask] <- 
        q.hat.lower[info.mask] <- q.hat.upper[info.mask] <-
        sigma2[info.mask] <- r.vec[info.mask] <- 0
      
    } else {
      
      X <- matrix(X[, info.mask], nrow = length(y))
      
      # Run the MDLQ model
      out <- tryCatch(
        { out <- ss.regress(y=y, X=X,
                            n.samples = 1200,
                            n.burn.in = 200)
        }, error = function(msg){
          out <- "DNC"
          return(out)
        }
      )
      
      if (length(out) == 1){
        q.hat[info.mask] <- 
          q.hat.lower[info.mask] <- q.hat.upper[info.mask] <-
          sigma2[info.mask] <- r.vec[info.mask] <- out
        
      } else {
        these.rates <- out[, colnames(out) %in% paste0("beta", 1:ncol(X))]
        these.sigma2 <- out[, colnames(out) == "sigma2"]
        these.r <- out[, colnames(out) == "r"]
        
        these.rates <- as.matrix(these.rates, ncol = ncol(X))
        
        q.hat[info.mask] <- apply(these.rates, 2, mean) 
        q.hat.lower[info.mask] <- apply(these.rates, 2, function(X) quantile(X, probs = 0.025)) 
        q.hat.upper[info.mask] <- apply(these.rates, 2, function(X) quantile(X, probs = 0.975)) 
        
        sigma2[info.mask] <- mean(these.sigma2)
        
        r.vec[info.mask] <- mean(these.r)
      }
    }
    
    # Save output
    to.save <- list(q.hat = q.hat, q.hat.lower = q.hat.lower, q.hat.upper = q.hat.upper,
                    sigma2 = sigma2, r.vec = r.vec)
    
    to.save
    
  } # End loop through inversion windows
  
  stopCluster(cl)
  
  exp.results[[zz]] <- big.out
  
} # end loop through ac.vals


# Package up MDLQ results
to.save <- list(times = times, 
                obs = obs,
                source.names = source.names,
                WD = data$WD,
                WS = data$WS,
                out = exp.results,
                Wellhead.West = data$Wellhead.West,
                Wellhead.East = data$Wellhead.East,
                Tanks = data$Tanks,
                Separator.West = data$Separator.West,
                Separator.East = data$Separator.East)

# Save results
saveRDS(to.save, output.file.path)

