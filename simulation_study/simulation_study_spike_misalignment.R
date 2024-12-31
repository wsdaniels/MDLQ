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

num.cores.to.use <- 4

forward.model.path         <- '../output_data/forward_model_output_ADED2024.RData'
output.file.path           <- '../output_data/simulation_study_output_spike_misalignment.RData'
spike.detection.alg.path   <- 'https://raw.github.com/wsdaniels/DLQ/master/code/HELPER_spike_detection_algorithm.R'
spike.slab.regression.path <- '../code/HELPER_spike_slab_regression.R'
helper.function.path       <- '../code/HELPER_functions.R'
leak.data.path             <- '../input_data/leak_data_clean.RData'

# END OF USER INPUT - NO MODIFICATION NECESSARY BELOW THIS POINT
#---------------------------------------------------------------------------


# Source helper files which contain helper functions 
source(spike.detection.alg.path)
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


# Read in controlled release data and organize
leak.data <- readRDS(leak.data.path)
leak.data[is.na(leak.data)] <- 0
leak.data <- leak.data[,c(1,2,3,4,6,7,5)]

# Filter out durations < 30
durations <- as.numeric(difftime(leak.data$end, leak.data$start, units = "min"))
to.keep <- durations >= 30
leak.data <- leak.data[to.keep, ]

# Filter out releases that end after data stops
to.remove <- leak.data$end > max(data$times)
leak.data <- leak.data[!to.remove, ]

# Get just the true release rates
leak.data.mat <- leak.data[,3:7]

# Number of releases
n.ints <- nrow(leak.data)

# True error variance to use in simulated observations
sigma2.true <- 2

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
  
  big.out <- foreach(t = 1:n.ints) %dopar% {
    
    this.mask <- which(data$times >= leak.data$start[t] & data$times <= leak.data$end[t])
    
    interval.length <- 30
    step.size <- 10
    num.intervals <- length(this.mask) / step.size - interval.length / step.size + 1
    num.intervals <- floor(num.intervals)
    mcmc.to.save <- data.to.save <- vector(mode = "list", length = num.intervals)
    
    betas <- t(as.matrix(leak.data.mat[t, ]))
    
    for (a in 1:num.intervals){
      
      time.mask <- seq((a-1)*step.size + 1,
                       (a-1)*step.size + interval.length)
      
      sub.mask <- this.mask[time.mask]
      
      N <- length(sub.mask)
      
      X <- matrix(NA, nrow = N*n.r, ncol = n.s)
      colnames(X) <- source.names
      
      for (s in 1:n.s){
        
        this.source <- c()
        for (r in 1:n.r){
          this.source <- c(this.source, sims[[s]][sub.mask, r])
        }
        
        X[,s] <- this.source
      }
      
      
      
      y.no.error <- X %*% betas
      y <- y.no.error + rnorm(nrow(X), 0, sqrt(sigma2.true))
      
      plot(y, type = "l")
      
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
      
      lines(y, col = "red")
      
      if (n.spikes.to.move > 0){
        
        for (k in spike.ind.to.move){
          
          this.spike.mask <- y.spikes$events == event.nums[k]
          
          this.y.to.add <- na.omit(y.no.error[this.spike.mask])
          
          this.start.ind <- round(runif(1, min = 1, max = length(y) - length(this.y.to.add) - 1))
          
          to.add.mask <- seq(this.start.ind, this.start.ind + length(this.y.to.add) -1)
          
          y[to.add.mask] <- y[to.add.mask] + this.y.to.add
        }
        
      }
      
      lines(y, col = "blue")
      
      out <- tryCatch(
        { out <- ss.regress(y=y, X=X,
                            n.samples = 3000)
        }, error = function(msg){
          out <- NA
          return(out)
        }
      )
      
      
      mcmc.to.save[[a]] <- out
      data.to.save[[a]] <- list(y = y, X = X)
      
    } # End loop through 30-minute intervals
    
    big.to.save <- list(mcmc.to.save, data.to.save)
    
    big.to.save 
    
  } # End loop through events
  
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

