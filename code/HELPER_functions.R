
remove.background <- function(obs, gap.time){
  
  # Skip sensors that have only NA values
  to.use <- which(apply(obs, 2, function(X) !all(is.na(X))))
  
  # Loop through sensor units
  for (j in to.use){
    
    # Grab observations from this sensor
    this.raw.obs <- obs[,j]
    
    # Remove the NA's at the beginning and end of the time series
    # (some sensors start late or end early)
    to.keep <- !is.na(this.raw.obs)
    trimmed.obs <- this.raw.obs[to.keep]
    trimmed.times <- times[to.keep]
    
    # Flag spikes
    spikes <- find.spikes(obs = trimmed.obs,
                          times = trimmed.times,
                          amp.threshold = 0.1,
                          make.plot = F)
    
    # Add points immediately before and after the spike to the spike mask
    # This better captures the full event, as there is some delay between the
    # start of an emission and when the sensors first see an enhancement
    for (i in na.omit(unique(spikes$events))){
      min.ind <- min(which(spikes$events == i))
      max.ind <- max(which(spikes$events == i))
      spikes$events[c(max(c(min.ind-1, 1)), min(c(max.ind+1,nrow(spikes))))] <- i
    }
    
    # Pull out integers that uniquely identify an event
    event.nums <- na.omit(unique(spikes$events))
    
    # If there are at least two events
    if (length(event.nums) > 1){
      
      # Loop through spikes and combine spikes that are separated by less than gap.time minutes
      for (i in 2:length(event.nums)){
        
        # Mask in this spike and last spike
        this.spike <- spikes$events == event.nums[i]
        previous.spike <- spikes$events == event.nums[i-1]
        
        # Clean up
        this.spike[is.na(this.spike)] <- F
        previous.spike[is.na(previous.spike)] <- F
        
        # Get start time of current spike and end time of previous spike
        this.spike.start.time <- spikes$time[this.spike][1]
        previous.spike.end.time <- spikes$time[previous.spike][length(spikes$time[previous.spike])]
        
        # Compute time difference
        time.diff <- difftime(this.spike.start.time, previous.spike.end.time, units = "mins")
        
        # Check gap
        if (minutes(time.diff) < minutes(gap.time)){
          spikes$events[this.spike] <- event.nums[i-1]
          event.nums[i] <- event.nums[i-1]
        }
      }
    }
    
    # Grab the new event numbers after combining events in the previous for loop
    event.nums <- na.omit(unique(spikes$events))
    
    # If there are any events
    if (length(event.nums) > 0){
      
      # Loop through events and (1) fill in any gaps between events with the correct
      # event number, and (2) estimate background as mean of first and last spike
      # points, which occur before the sharp increase and after the spike has
      # returned to return.threshold percent of the max value
      for (i in 1:length(event.nums)){
        
        # Fill in gaps
        first.ob <- min(which(spikes$events == event.nums[i]))
        last.ob <- max(which(spikes$events == event.nums[i]))
        this.mask <- first.ob:last.ob
        spikes$events[this.mask] <- event.nums[i]
        
        # Estimate background using observations before and after spike
        # This is just the first and last observation within the spike mask,
        # since we already added the points before and after the spike to the
        # spike mask earlier on
        b.left <- trimmed.obs[first.ob]
        b.right <- trimmed.obs[last.ob]
        b <- mean(c(b.left, b.right))
        
        # Remove background from this spike
        trimmed.obs[this.mask] <- trimmed.obs[this.mask] - b
      }
      
    }
    
    # Remove background from all non-spike data. Since by definition these points
    # are not in an event, their concentration value is directly taken as the
    # background estimate. Hence removing background is just setting the value to zero
    trimmed.obs[is.na(spikes$events)] <- 0
    
    # Set any negative values to zero. Negative value would arise if the
    # background estimate was too large (greater than actual concentration value)
    trimmed.obs[trimmed.obs < 0] <- 0
    
    # Save background removed data
    obs[to.keep, j] <- trimmed.obs
  }
  
  return(obs)
  
}





perform.event.detection <- function(obs, gap.time, length.threshold){
  
  # Compute maximum concentration value across sensors for each minute
  # Note that since all non-spike observations have been set to zero, any 
  # non-zero value in the max.obs time series should be considered an event 
  # (and hence no need to run the spike detection algorithm again)
  max.obs <- apply(obs, 1, max, na.rm = T)
  
  # Create data frame with time steps and event mask
  spikes <- data.frame(time = times, events = max.obs > 0)
  
  # Find gaps between events that are shorter than gap.time and turn them into events
  #   to.replace holds the indices that need to be switched from F to T
  #   first.gap is an indicator for the first gap (which should not be replaced,
  #   regardless of length)
  #   false.seq holds the indices of each sequence of FALSEs (non-events)
  to.replace <- c()
  first.gap <- T
  false.seq <- c()
  
  # Loop through times
  for (i in 1:length(times)){
    
    # If not a spike, add index to false.sequence
    if (!spikes$events[i]){
      false.seq <- c(false.seq, i)
      
      # Otherwise, check length of false sequence, if greater than gap time,
      # save those indices to replace later
    } else if (spikes$events[i]){
      
      if (length(false.seq) <= gap.time & !first.gap){
        to.replace <- c(to.replace, false.seq)
      }
      
      first.gap <- F
      false.seq <- c()
    }
  }
  
  # Replace gaps shorter than gap.time with T (meaning they are in an event)
  spikes$events[to.replace] <- T
  
  # Now we replace the T/F with an integer to distinguish between events
  # Start by replacing F with NA and T with zero
  spikes$events[!spikes$events] <- NA
  spikes$events[spikes$events] <- 0
  
  # Get indices of spikes
  spike.points <- which(!is.na(spikes$events))
  count <- 0
  
  # Loop through spike points, if last point was not a spike, increase counter
  for (i in spike.points){
    
    if (is.na(spikes$events[i-1])){
      count <- count + 1
      spikes$events[i] <- count
    } else {
      spikes$events[i] <- count
    }
  }
  
  # Get integers that uniquely define the different events
  event.nums <- na.omit(unique(spikes$events))
  
  # Filter events by the length threshold
  for (i in 1:length(event.nums)){
    this.spike <- which(spikes$events == event.nums[i])
    
    if (length(this.spike) < length.threshold){
      spikes$events[this.spike] <- NA
    }
  }
  
  return(spikes)
}