# Description: Performs multi-source emission event detection, localization,
#              and quantification using output from the Gaussian puff 
#              atmospheric dispersion model and observations from point sensor
#              networks.
# Author: William Daniels (wdaniels@mines.edu)
# Last Updated: December 30, 2024

# Clear environment
if(!is.null(dev.list())){dev.off()}
rm(list = ls())

# Import necessary libraries
library(lubridate)
library(foreach)
library(doParallel)
library(zoo)

# Start code timer
code.start.time <- Sys.time()

# START USER INPUT
#---------------------------------------------------------------------------

# Number of cores to use. Set equal to 1 to run in serial.
num.cores.to.use <- 4

# Length of inversion window [minutes]
interval.length <- 30

# Number of minutes to advance the inversion window
step.size <- 30

# Path to output from atmospheric dispersion model
forward.model.path <- '../input_data/forward_model_output_ADED2024.RData'

# Location to save MDLQ output
output.file.path <- '../output_data/MDLQ_output_ADED2024_30min_interval_30min_step.RData'

# Path to helper file that contains the Gibbs updates for the MDLQ model
spike.slab.regression.path <- '../code/HELPER_MCMC.R'

# Path to helper file that contains functions from: https://doi.org/10.1525/elementa.2023.00110
# These functions remove background concentrations and perform event detection
helper.function.path <- '../code/HELPER_functions.R'

# Index (starting at 1) of the first simulation file in the "data" object
first.sim.ind <- 5


# END OF USER INPUT - NO MODIFICATION NECESSARY BELOW THIS POINT
#---------------------------------------------------------------------------


# STEP 1: READ IN DATA AND PRELIMINARIES
#---------------------------------------------------------------------------

# Source helper files 
source('https://raw.github.com/wsdaniels/DLQ/master/code/HELPER_spike_detection_algorithm.R')
source(spike.slab.regression.path)
source(helper.function.path)

# Read in simulation data
data <- readRDS(forward.model.path)

# Trim data so that they start and end at either the hour or half hour. 
# This makes it so that all 30-minute intervals at aligned 
# This step was written by Spencer Kidd and added June 11, 2025
first.clean.time <- min(which(minute(data$times) %in% c(0,30)))
last.clean.time <- max(which(minute(data$times) %in% c(0,30)))
to.keep <- first.clean.time:last.clean.time
data$times <- data$times[to.keep]
data$WD <- data$WD[to.keep]
data$WS <- data$WS[to.keep]
data$obs <- data$obs[to.keep, ]
data[first.sim.ind:length(data)] <- lapply(data[first.sim.ind:length(data)], function(X) X[to.keep, ])


# Pull out sensor observations and replace NA's that are not on edge of 
# the time series with interpolated values
obs <- na.approx(data$obs, na.rm = F, maxgap = 10)

# Number of sensors
n.r <- ncol(obs)

# Pull out time stamps of observations and simulations
times <- data$times

# Pull out the simulation predictions
sims <- data[first.sim.ind:length(data)]

# Grab source info
n.s <- length(sims) 
source.names <- names(sims)



# STEP 2: REMOVE BACKGROUND FROM OBSERVATIONS
#---------------------------------------------------------------------------

# Remove background using the method from: https://doi.org/10.1525/elementa.2023.00110
obs <- remove.background(obs, gap.time = 30)


# STEP 4: METHANE EMISSION SOURCE APPORTIONMENT
# ---------------------------------------------------------------------------

# Compute total number of intervals
num.intervals <- (length(times)/step.size) - (interval.length/step.size) + 1
num.intervals <- floor(num.intervals)

# Fire up the parallel cluster
cl <- makeCluster(num.cores.to.use)
registerDoParallel(cl)

# Loop through 30-minute intervals
big.out <- foreach(a = 1:num.intervals) %dopar% {
  
  # Initialize vectors to hold output
  q.hat <- q.hat.lower <- q.hat.upper <- pis <- q.hat.median <- vector(length = n.s)
  these.rates <- q.max <- NA
  
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
  
  # Check for NA y values
  if (any(is.na(y))){ 
    print("NA y")
    
  } else {
    
    if(T){
      plot(y, type = "l", ylim = c(0,50), lwd = 3)
      for (p in 1:length(q.hat)){
        if (!is.na(as.numeric(q.hat[p]))){
          lines(X[,p], col = p+1)    
        }}}
    
    # Determine which sources have downwind sensors (i.e., which sources have information)
    info.mask <- apply(X, 2, function(this.col) sum(this.col > 0.25) > 4)
    
    # Save information / no information mask
    q.hat[!info.mask] <-
      q.hat.lower[!info.mask] <- q.hat.upper[!info.mask] <-
      pis[!info.mask] <- q.hat.median[!info.mask] <- "no info"
    
    # Subset to just the sources with information
    X <- matrix(X[, info.mask], nrow = length(y))
    
    # If mostly zeros, set emission rate for sources with information to zero
    # OLD VERSION THAT WAS TOO STRICT: if (all(y == 0)){
    if (sum(y > 0) < 5 & all(y[y > 0] < 1)){
      
      q.hat[info.mask] <-
        q.hat.lower[info.mask] <- q.hat.upper[info.mask] <-
        pis[info.mask] <- q.hat.median[info.mask] <- 0
      
    } else {
      
      # Determine which sources have overlap between simulated enhancements and enhancements in observations
      dot <- as.vector(crossprod(y, X))
      norm2 <- colSums(X^2, na.rm = T)
      q.max <- 2 * dot / norm2
      overlap.mask <- q.max > (0.1 / 3.6) # q of at least 0.1 kg/hr
      
      # Save overlap / no overlap mask
      q.hat[info.mask][!overlap.mask] <-
        q.hat.lower[info.mask][!overlap.mask] <- q.hat.upper[info.mask][!overlap.mask] <-
        pis[info.mask][!overlap.mask] <- q.hat.median[info.mask][!overlap.mask] <- "no overlap"
      
      # Subset to just the sources with overlap
      X <- matrix(X[, overlap.mask], nrow = length(y))
      
      # Run the MDLQ model if any sources are left
      if (ncol(X) > 0){
        
        # q.hat[info.mask] <-
        #   q.hat.lower[info.mask] <- q.hat.upper[info.mask] <-
        #   pis[info.mask] <- q.hat.median[info.mask] <- "invert"
        
        out <- tryCatch(
          { out <- run.mdlq.mcmc(y=y, X=X,
                                 n.samples = 750,
                                 n.burn.in = 250)
          }, error = function(msg){
            out <- "DNC"
            return(out)
          }
        )
        
        # catch any errors
        if (length(out) == 1){
          q.hat[info.mask][overlap.mask] <-
            q.hat.lower[info.mask][overlap.mask] <-
            q.hat.upper[info.mask][overlap.mask] <-
            pis[info.mask][overlap.mask] <-
            q.hat.median[info.mask][overlap.mask] <- out
          
          # Save output
        } else {
          these.rates <- out[, colnames(out) %in% paste0("beta", 1:ncol(X))]
          these.pis   <- out[, colnames(out) %in% paste0("pi",   1:ncol(X))]
          
          these.rates <- as.matrix(these.rates, ncol = ncol(X))
          these.pis   <- as.matrix(these.pis,   ncol = ncol(X))
          
          q.hat[info.mask][overlap.mask]       <- apply(these.rates, 2, mean) * 3.6
          q.hat.lower[info.mask][overlap.mask] <- apply(these.rates, 2, function(X) quantile(X, probs = 0.025)) * 3.6
          q.hat.upper[info.mask][overlap.mask] <- apply(these.rates, 2, function(X) quantile(X, probs = 0.975)) * 3.6
          
          pis[info.mask][overlap.mask] <- apply(these.pis, 2, mean)
          
          q.hat.median[info.mask][overlap.mask] <- apply(these.rates, 2, median) * 3.6
          
        } # End if to check for MCMC convergence
      } # End if to check if there are any columns of X left
    } # End if to check for enhancements in observations
    
    if(F){
      plot(y, type = "l", ylim = c(0,30), lwd = 3)
      counter <- 1
      for (p in 1:length(q.hat)){
        if (!is.na(as.numeric(q.hat[p]))){
          lines(X[,counter] * as.numeric(q.hat[p])/3.6, col = p+1)    
          counter <- counter + 1
        }}}
    
    # Save output
    to.save <- list(q.hat = q.hat, q.hat.lower = q.hat.lower, q.hat.upper = q.hat.upper,
                    pis = pis, q.hat.median = q.hat.median, rates = these.rates * 3.6, q.max = q.max)
    
    to.save
    
  } # End check for NA y values
  
} # End loop through 30-minute intervals

# Stop the parallel cluster
stopCluster(cl)



# Package up MDLQ results
to.save <- c(list(times = times,
                  obs = obs,
                  source.names = source.names,
                  WD = data$WD,
                  WS = data$WS,
                  out = big.out),
             sims)

# Save results
saveRDS(to.save, output.file.path)


# End code timer
code.stop.time <- Sys.time() 

# Print wall clock execution time
difftime(code.stop.time, code.start.time, units = "mins")