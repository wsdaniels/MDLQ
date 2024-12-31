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


# START USER INPUT
#---------------------------------------------------------------------------

# When removing background, spikes that are separated by gap.time minutes or less
# will be combined into one event.
gap.time <- 30

# Number of cores to use. Set equal to 1 to run in serial.
num.cores.to.use <- 6

# Set how to run the MDLQ. Options include:
#  'run.on.releases' - runs the MDLQ using the true start and end times of the controlled releases. USE THIS TO REPRODUCE MANUSCRIPT RESULTS.
#  'event.detection' - runs the MDLQ using estimated start and end times of the controlled releases.
#  '30.min' - runs the MDLQ on overlapping 30-minute intervals. 
run.mode <- 'event.detection'

# Path to output from atmospheric dispersion model
forward.model.path <- '../input_data/forward_model_output_ADED2024.RData'

# Location to save MDLQ output
output.file.path <- '../output_data/MDLQ_output_ADED2024.RData'

# Path to helper file that contains the Gibbs updates for the MDLQ model
spike.slab.regression.path <- '../code/HELPER_spike_slab_regression.R'

# Path to helper file that contains functions from: https://doi.org/10.1525/elementa.2023.00110
# These functions remove background concentrations and perform event detection
helper.function.path <- '../code/HELPER_functions.R'

# Path to controlled release data. Necessary if run.mode = run.on.releases
leak.data.path <- '../input_data/leak_data.RData'

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


# STEP 2: REMOVE BACKGROUND FROM OBSERVATIONS
#---------------------------------------------------------------------------

# Remove background using the method from: https://doi.org/10.1525/elementa.2023.00110
obs <- remove.background(obs, gap.time)




# STEP 3: DEFINE EMISSION EVENTS
#---------------------------------------------------------------------------

# If running on exact timing of controlled releases, get controlled release data
if (run.mode == 'run.on.releases'){
  
  # Read in controlled release data
  leak.data <- readRDS(leak.data.path)
  
  # NA's mean no emissions, so replace with a 0 kg/hr value
  leak.data[is.na(leak.data)] <- 0
  
  # Reorganize columns
  leak.data <- leak.data[,c(1,2,3,4,6,7,5)]
  
  # Get controlled release durations
  durations <- as.numeric(difftime(leak.data$end, leak.data$start, units = "min"))
  
  # Filter out events less than 30 minutes
  to.keep <- durations >= 30
  leak.data <- leak.data[to.keep, ]
  
  # Get just the emission rates (not times)
  leak.data.mat <- leak.data[,3:7]
  
  # Remove controlled releases that occurred after CMS data stops
  to.remove <- leak.data$end > max(data$times)
  leak.data <- leak.data[!to.remove, ]
  
  # Number of emission events
  n.ints <- nrow(leak.data)
  
  # If estimating start and end times, run event detection method from: https://doi.org/10.1525/elementa.2023.00110 
} else if (run.mode == 'event.detection'){
  
  # Perform event detection 
  spikes <- perform.event.detection(obs, gap.time, length.threshold = 15)
  
  # Grab event number again after filtering by length
  event.nums <- na.omit(unique(spikes$events))
  
  # Number of emission events
  n.ints <- length(event.nums)
  
}




# STEP 4: METHANE EMISSION SOURCE APPORTIONMENT
# ---------------------------------------------------------------------------


if (run.mode == '30.min'){
  
  # Quantification interval length [minutes]
  interval.length <- 30
  
  # Number of minutes to advance the quantification interval 
  step.size <- 10
  
  # Compute total number of intervals
  num.intervals <- (length(times)/step.size) - (interval.length/step.size) + 1
  num.intervals <- floor(num.intervals)
  
  # Fire up the parallel cluster
  cl <- makeCluster(num.cores.to.use)
  registerDoParallel(cl)
  
  # Loop through 30-minute intervals
  big.out <- foreach(a = 1:num.intervals) %dopar% {
    
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
    
    # plot(y, type = "l", ylim = c(0,15))
    # lines(X[,1], col = "red")
    # lines(X[,2], col = "blue")
    # lines(X[,3], col = "green")
    # lines(X[,4], col = "orange")
    # lines(X[,5], col = "purple")
    
    
    code.start.time <- Sys.time() 
    
    # Run the MDLQ model
    out <- tryCatch(
      { out <- ss.regress(y=y, X=X,
                          n.samples = 3000,
                          n.burn.in = 750)
      }, error = function(msg){
        out <- NA
        return(out)
      }
    )
    
    code.end.time <- Sys.time() 
    difftime(code.end.time, code.start.time, units = "mins")
    
    # Save outputs and inputs
    big.to.save <- list(out, list(y = y, X = X))
    big.to.save
    
  } # End loop through 30-minute intervals
  
  # Stop the parallel cluster
  stopCluster(cl)
  
} else {
  
  # Fire up the parallel cluster
  cl <- makeCluster(num.cores.to.use)
  registerDoParallel(cl)
  
  # Loop through emission events
  big.out <- foreach(t = 1:n.ints) %dopar% {
    
    # Mask in the time steps of this event
    if (run.mode == 'run.on.releases'){
      this.mask <- which(data$times >= leak.data$start[t] & data$times <= leak.data$end[t])
      
    } else if (run.mode == 'event.detection'){
      this.mask <- seq(min(which(spikes$events == event.nums[t])),
                       max(which(spikes$events == event.nums[t])))
    }
    
    # Quantification interval length [minutes]
    interval.length <- 30
    
    # Number of minutes to advance the quantification interval 
    step.size <- 10
    
    # Compute total number of intervals
    num.intervals <- length(this.mask) / step.size - interval.length / step.size + 1
    num.intervals <- floor(num.intervals)
    
    # Initialize variables to hold MCMC output and input data
    mcmc.to.save <- data.to.save <- vector(mode = "list", length = num.intervals)
    
    # Loop through the number of quantification intervals within this event
    for (a in 1:num.intervals){
      
      print(paste0(a, "/", num.intervals))
      
      # Mask in this quantification interval
      time.mask <- seq((a-1)*step.size + 1,
                       (a-1)*step.size + interval.length)
      
      # Create sub-mask of the "this.mask" variable corresponding to this quantification interval
      if (run.mode == 'run.on.releases'){
        sub.mask <- this.mask[time.mask]
        
      } else if (run.mode == 'event.detection'){
        if (length(time.mask) > length(this.mask)){
          sub.mask <- this.mask
        } else {
          sub.mask <- this.mask[time.mask]  
        }
      }
      
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
      
      # plot(y, type = "l", ylim = c(0,15))
      # lines(X[,1], col = "red")
      # lines(X[,2], col = "blue")
      # lines(X[,3], col = "green")
      # lines(X[,4], col = "orange")
      # lines(X[,5], col = "purple")
      
      
      code.start.time <- Sys.time() 
      
      # Run the MDLQ model
      out <- tryCatch(
        { out <- ss.regress(y=y, X=X,
                            n.samples = 3000,
                            n.burn.in = 750)
        }, error = function(msg){
          out <- NA
          return(out)
        }
      )
      
      code.end.time <- Sys.time() 
      difftime(code.end.time, code.start.time, units = "mins")
      
      # Save outputs and inputs
      mcmc.to.save[[a]] <- out
      data.to.save[[a]] <- list(y = y, X = X)
      
    } # End loop through 30-minute intervals
    
    big.to.save <- list(mcmc.to.save, data.to.save)
    big.to.save 
    
  } # End loop through events
  
  # Stop the parallel cluster
  stopCluster(cl)
  
}

# Package up MDLQ results
to.save <- list(times = times, 
                obs = obs,
                source.names = source.names,
                WD = data$WD,
                WS = data$WS,
                out = big.out,
                Wellhead.West = data$Wellhead.West,
                Wellhead.East = data$Wellhead.East,
                Tanks = data$Tanks,
                Separator.West = data$Separator.West,
                Separator.East = data$Separator.East)

# Save results
saveRDS(to.save, output.file.path)
