# Description: 
# Summarizes event detection, localization, and quantification results and 
# reproduces paper figures.
# Author: William Daniels (wdaniels@mines.edu)
# Last Updated: December 2024

# Clear environment
if(!is.null(dev.list())){dev.off()}
rm(list = ls())
gc()

set.seed(1234)

# Import necessary libraries
library(ald)
library(fields)
library(scales)
library(lubridate)

if (commandArgs()[1] == "RStudio"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}


# START USER INPUT
#---------------------------------------------------------------------------

# Size of the inversion window used to run the MDLQ
interval.length <- 30

# Amount of time [minutes] that the inversion windows are moved forward
step.size <- 30

# Location of MDLQ output 
data <- readRDS('../output_data/MDLQ_output_ADED2024_30min_interval_30min_step.RData')

# Location of the METEC controlled release ground truth data
leak.data <- readRDS('../input_data/leak_data.RData')

# Generate ground truth data on 30-minute intervals (T), or read in a previously generated file (F)
generate.truth <- F
truth.location <- '../output_data/truth_ADED2024_30min_interval_30min_step.RData'

# Generate info mask on 30-minute intervals (T), or read in a previously generated file (F)
generate.info.mask <- F
info.mask.location <- '../output_data/info_mask_ADED2024_30min_interval_30min_step.RData'

# Path to the forward model output
forward.model.path <- '../input_data/forward_model_output_ADED2024.RData'

# END OF USER INPUT - NO MODIFICATION NECESSARY BELOW THIS POINT
#---------------------------------------------------------------------------




# STEP 1: PARSE OUT RESULTS
#---------------------------------------------------------------------------

# Read in the forward model output
orig.data <- readRDS(forward.model.path)

# Parse out source names
source.names <- data$source.names

q.hat <- lapply(data$out, function(X) X[[1]])
q.hat <- do.call(rbind, q.hat)

q.hat.lower <- lapply(data$out, function(X) X[[2]])
q.hat.lower <- do.call(rbind, q.hat.lower)

q.hat.upper <- lapply(data$out, function(X) X[[3]])
q.hat.upper <- do.call(rbind, q.hat.upper)

pis <- lapply(data$out, function(X) X[[4]])
pis <- do.call(rbind, pis)

q.hat.median <- lapply(data$out, function(X) X[[5]])
q.hat.median <- do.call(rbind, q.hat.median)


sum(q.hat == "no info")/(nrow(q.hat)*ncol(q.hat))
sum(q.hat == "broke on betas")/(nrow(q.hat)*ncol(q.hat))
sum(q.hat == "decreasing betas")/(nrow(q.hat)*ncol(q.hat))
sum(q.hat == "no overlap")/(nrow(q.hat)*ncol(q.hat))


ni.mask <- q.hat == "no info"
dnc.mask <- q.hat == "broke on betas" | q.hat == "no overlap"

q.hat.lower[q.hat.lower == "no info"] <- NA
q.hat.lower[q.hat.lower == "broke on betas"] <- NA
q.hat.lower[q.hat.lower == "decreasing betas"] <- NA
q.hat.lower[q.hat.lower == "no overlap"] <- NA
q.hat.lower <- matrix(as.numeric(q.hat.lower), ncol = ncol(q.hat.lower))

q.hat.upper[q.hat.upper == "no info"] <- NA
q.hat.upper[q.hat.upper == "broke on betas"] <- NA
q.hat.upper[q.hat.upper == "decreasing betas"] <- NA
q.hat.upper[q.hat.upper == "no overlap"] <- NA
q.hat.upper <- matrix(as.numeric(q.hat.upper), ncol = ncol(q.hat.upper))

pis[pis == "no info"] <- NA
pis[pis == "broke on betas"] <- NA
pis[pis == "decreasing betas"] <- NA
pis[pis == "no overlap"] <- NA
pis <- matrix(as.numeric(pis), ncol = ncol(pis))

q.hat[q.hat == "no info"] <- NA
q.hat[q.hat == "broke on betas"] <- NA
q.hat[q.hat == "decreasing betas"] <- NA
q.hat[q.hat == "no overlap"] <- NA
q.hat <- matrix(as.numeric(q.hat), ncol = ncol(q.hat))

q.hat.median[q.hat.median == "no info"] <- NA
q.hat.median[q.hat.median == "broke on betas"] <- NA
q.hat.median[q.hat.median == "decreasing betas"] <- NA
q.hat.median[q.hat.median == "no overlap"] <- NA
q.hat.median <- matrix(as.numeric(q.hat.median), ncol = ncol(q.hat.median))

mcmc <- lapply(data$out, function(X) X[[6]])

# Some intervals subset to only the sources with downwind sensors.
# This results in matrices of different dimensions.
# Here we get the MCMC output matrices to all have the same dimension, inputting NA
# values when a given source was omitted.
mcmc.correct.dim <- mcmc
for (i in 1:length(mcmc)){
  if (length(mcmc[[i]]) > 1){
    tmp <- matrix(NA, nrow = nrow(mcmc[[i]]), ncol = length(source.names))
    this.im <- !ni.mask[i, ]
    tmp[, this.im] <- mcmc[[i]]
    mcmc.correct.dim[[i]] <- tmp
  } 
}


times <- data$times
sims <- data[7:11]

num.intervals <- (length(times)/step.size) - (interval.length/step.size) + 1
num.intervals <- floor(num.intervals)

source('https://raw.github.com/wsdaniels/DLQ/master/code/HELPER_spike_detection_algorithm.R')
source('https://raw.github.com/wsdaniels/CMS-durations/main/code/helper_functions.R')



# STEP 2: READ IN CONTROLLED RELEASE DATA
#---------------------------------------------------------------------------

# NA's mean no emissions, so replace with a 0 kg/hr value
leak.data[is.na(leak.data)] <- 0

# Rename leak columns to match data$source.names
original.name <- c("EG4W", "EG4S", "EG5S", "EG4T", "EG5W")
matching.name <- c("Wellhead.West", "Separator.West", "Separator.East", "Tanks", "Wellhead.East")
colnames(leak.data)[3:7] <- matching.name[match(colnames(leak.data)[3:7], original.name)]

this.order <- match(data$source.names, colnames(leak.data)[3:7])

# Reorganize columns of data matrix to match order of sources in data$source.names
leak.data <- leak.data[ , c(1,2,this.order+2)]

# Get controlled release start and end times
leak.data$start <- as_datetime(leak.data$start, tz = "America/Denver")
leak.data$end   <- as_datetime(leak.data$end, tz = "America/Denver")

# Get true site-level release rates
leak.total <- apply(leak.data[, 3:7], 1, sum)
leak.data.mat <- as.matrix(leak.data[,c(3:7)])




# STEP 3: GENERATE INFORMATION MASK
#---------------------------------------------------------------------------

if (generate.info.mask){
  for (i in 1:num.intervals){
    print(paste0(i, "/", num.intervals))
    # scale.factors <- q.hat[i, ]/3.6
    time.mask <- seq((i-1)*step.size + 1,
                     (i-1)*step.size + interval.length)
    for (j in 1:length(source.names)){
      # sims[[j]][time.mask, ] <- sims[[j]][time.mask, ] * scale.factors[j]
      sims[[j]][time.mask, ] <- sims[[j]][time.mask, ] 
    }
  }
  info.mask <- create.info.mask(times, sims)
  saveRDS(info.mask, info.mask.location)
} else {
  info.mask <- readRDS(info.mask.location)
}

# Translate information mask to the 30-minute intervals
info <- matrix(NA, nrow = num.intervals, ncol = length(source.names))
for (j in 1:length(source.names)){
  this.info.mask <- info.mask[[which(names(sims) == source.names[j])]]
  for (i in 1:num.intervals){
    time.mask <- seq((i-1)*step.size + 1,
                     (i-1)*step.size + interval.length)
    info[i,j] <- sum(!is.na(this.info.mask$events[time.mask]))/interval.length
  }
}

# Mask in only the 30-minute windows that have full information for all sources
# info.to.use <- apply(info, 1, function(X) all(X == 1))
info.to.use <- apply(info, 1, function(X) all(X > 0.95))



# STEP 4: AVERAGE MDLQ OUTPUT WHEN 30-MINUTE WINDOWS OVERLAP
#---------------------------------------------------------------------------

num.to.avg <- (interval.length/step.size) - 1
q.hat.avg <- q.hat.lower.avg <- q.hat.upper.avg <- matrix(NA, nrow = num.intervals, ncol = length(source.names))
q.hat.avg[1, ] <- q.hat[1, ]
q.hat.lower.avg[1, ] <- q.hat.lower[1, ]
q.hat.upper.avg[1, ] <- q.hat.upper[1, ]

if (interval.length > 10){
  for (i in 2:num.intervals){
    ind.to.avg <- seq(max(1, i-num.to.avg), i)
    q.hat.avg[i, ] <- apply(matrix(q.hat[ind.to.avg, ], ncol = length(source.names)), 2, mean, na.rm = T)
    q.hat.lower.avg[i, ] <- apply(matrix(q.hat.lower[ind.to.avg, ], ncol = length(source.names)), 2, mean, na.rm = T)
    q.hat.upper.avg[i, ] <- apply(matrix(q.hat.upper[ind.to.avg, ], ncol = length(source.names)), 2, mean, na.rm = T)
  }
  
} else {
  for (i in 2:num.intervals){
    ind.to.avg <- seq(max(1, i-num.to.avg), i)
    q.hat.avg[i, ] <- q.hat[ind.to.avg, ]
    q.hat.lower.avg[i, ] <- q.hat.lower[ind.to.avg, ]
    q.hat.upper.avg[i, ] <- q.hat.upper[ind.to.avg, ]
  }
}



# STEP 5: GET GROUND TRUTH ON THE 30-MINUTE INTERVALS
#---------------------------------------------------------------------------

get.times <- function(a, step.size){
  sub.mask <- seq((a-1)*step.size + 1,
                  (a-1)*step.size + step.size)
  return(times[sub.mask])
}

if (generate.truth){
  truth <- matrix(NA, nrow = num.intervals, ncol = length(source.names))
  
  for (i in 1:num.intervals){
    print(paste0(i, "/", num.intervals))
    
    int.times <- get.times(i, step.size)
    which.leak <- rep(NA, length(int.times))
    
    for (j in 1:nrow(leak.data)){
      leak.interval <- interval(leak.data$start[j], leak.data$end[j])
      overlap.mask <- int.times %within% leak.interval
      
      which.leak[overlap.mask] <- j
    }
    
    vals.to.avg <- matrix(NA, nrow = length(which.leak), ncol = length(source.names))
    for (k in 1:nrow(vals.to.avg)){
      if (is.na(which.leak[k])){
        vals.to.avg[k, ] <- rep(0, length(source.names))
      } else {
        vals.to.avg[k, ] <- leak.data.mat[which.leak[k], ]
      }
    }
    truth[i, ] <- apply(vals.to.avg, 2, mean)
  }
  saveRDS(truth, truth.location)
  
} else {
  truth <- readRDS(truth.location)
}


# STEP 6: COMPUTE SITE-LEVEL EMISSION RATE ESTIMATES AND COMPUTE CONVERAGES
#---------------------------------------------------------------------------

# Compute site-level emission rate estimates
q.hat.total <- apply(q.hat.avg, 1, sum)
q.hat.total.lower <- apply(q.hat.lower.avg, 1, sum)
q.hat.total.upper <- apply(q.hat.upper.avg, 1, sum)
truth.total <- apply(truth, 1, sum)

# Compute site-level emission rate coverage
coverage <- sum(truth.total >= q.hat.total.lower & truth.total <= q.hat.total.upper, na.rm = T) / sum(!is.na(q.hat.total))

# Compute source-level emission rate coverage
source.level.coverage <- vector(length = length(source.names))
for (i in 1:length(source.names)){
  source.level.coverage[i] <- sum(truth[,i] >= q.hat.lower.avg[,i] & truth[,i] <= q.hat.upper.avg[,i], na.rm = T) / sum(!is.na(q.hat.lower.avg[,i]))
}




# STEP 7: CREATE RESULT TIME SERIES AND DATA EXAMPLE FIGURES
#---------------------------------------------------------------------------

# Set colors for plots
tank.color <- "#3062CF" #blue
wellhead.east.color <- "#9147B8" #purple
separator.east.color <- "#7EAD52" #green
separator.west.color <- "#F1C30E" #gold
wellhead.west.color <- "#C7383C" #red

cols <- c(wellhead.west.color, separator.west.color, tank.color, wellhead.east.color, separator.east.color)

interval.times <- vector(length = num.intervals)
for (i in 1:num.intervals){
  interval.times[i] <- get.times(i, step.size)[1]
}
interval.times <- as_datetime(interval.times, tz = "America/Denver")


png('../figures/result_time_series.png',
    res = 100, pointsize = 24, width = 1920, height = 1080*0.9)

par(mfrow = c(2,1))
par(mar = c(1.5,0,0,0))
par(oma = c(0.5,3.5,0.5,0.5))
par(mgp = c(2.5, 0.75, 0))

truth.cum <- matrix(NA, nrow = nrow(truth), ncol = ncol(truth))
for (i in 1:nrow(truth)){
  truth.cum[i, ] <- cumsum(truth[i, ])
}

q.hat.cum <- matrix(NA, nrow = nrow(q.hat.avg), ncol = ncol(q.hat.avg))
for (i in 1:nrow(q.hat.cum)){
  q.hat.cum[i, ] <- cumsum(q.hat.avg[i, ])
}

# this.mask <- 1040:1770
this.mask <- 347:592

ylim.max <- 10

truth.alpha.val <- 0.8
plot(interval.times[this.mask], truth.cum[this.mask,1], type = "l", ylim = c(0,ylim.max), col = NA, xaxt= "n")
envelopePlot(x1 = interval.times[this.mask], y1 = rep(0, length(this.mask)), y2 = truth.cum[this.mask, 1],
             lineCol = NA, col = alpha(cols[1], truth.alpha.val))
envelopePlot(x1 = interval.times[this.mask], y1 = truth.cum[this.mask, 1], y2 = truth.cum[this.mask, 2],
             lineCol = NA, col = alpha(cols[2], truth.alpha.val))
envelopePlot(x1 = interval.times[this.mask], y1 = truth.cum[this.mask, 2], y2 = truth.cum[this.mask, 3],
             lineCol = NA, col = alpha(cols[3], truth.alpha.val))
envelopePlot(x1 = interval.times[this.mask], y1 = truth.cum[this.mask, 3], y2 = truth.cum[this.mask, 4],
             lineCol = NA, col = alpha(cols[4], truth.alpha.val))
envelopePlot(x1 = interval.times[this.mask], y1 = truth.cum[this.mask, 4], y2 = truth.cum[this.mask, 5],
             lineCol = NA, col = alpha(cols[5], truth.alpha.val))
lines(interval.times[this.mask], truth.total[this.mask], lwd = 3)
mtext("Emission Rate [kg/hr]", side = 2, line = 2.25)

alpha.val <- 0.8
plot(interval.times[this.mask], q.hat.cum[this.mask,1], col = NA, type = "l", ylim = c(0,ylim.max), xaxt = "n")
envelopePlot(x1 = interval.times[this.mask], y1 = rep(0, length(this.mask)),
             y2 = ifelse(is.na(q.hat.cum[this.mask, 1]), 0, q.hat.cum[this.mask, 1]),
             lineCol = NA, col = alpha(cols[1], alpha.val))
envelopePlot(x1 = interval.times[this.mask], y1 = ifelse(is.na(q.hat.cum[this.mask, 1]), 0, q.hat.cum[this.mask, 1]),
             y2 = ifelse(is.na(q.hat.cum[this.mask, 2]), 0, q.hat.cum[this.mask, 2]),
             lineCol = NA, col = alpha(cols[2], alpha.val))
envelopePlot(x1 = interval.times[this.mask], y1 = ifelse(is.na(q.hat.cum[this.mask, 2]), 0, q.hat.cum[this.mask, 2]),
             y2 = ifelse(is.na(q.hat.cum[this.mask, 3]), 0, q.hat.cum[this.mask, 3]),
             lineCol = NA, col = alpha(cols[3], alpha.val))
envelopePlot(x1 = interval.times[this.mask], y1 = ifelse(is.na(q.hat.cum[this.mask, 3]), 0, q.hat.cum[this.mask, 3]),
             y2 = ifelse(is.na(q.hat.cum[this.mask, 4]), 0, q.hat.cum[this.mask, 4]),
             lineCol = NA, col = alpha(cols[4], alpha.val))
envelopePlot(x1 = interval.times[this.mask], y1 = ifelse(is.na(q.hat.cum[this.mask, 4]), 0, q.hat.cum[this.mask, 4]),
             y2 = ifelse(is.na(q.hat.cum[this.mask, 5]), 0, q.hat.cum[this.mask, 5]),
             lineCol = NA, col = alpha(cols[5], alpha.val))
lines(interval.times[this.mask], truth.total[this.mask], lwd = 3)

mtext("Emission Rate [kg/hr]", side = 2, line = 2.25)

date.seq <- seq(round_date(interval.times[1], unit = "days"), 
                round_date(interval.times[length(interval.times)], unit = "days"),
                by = "1 day")

axis(side = 1, 
     at = date.seq,
     labels = paste0(month.abb[month(date.seq)]," ", day(date.seq)))

dev.off()



leak.cum <- matrix(NA, nrow = nrow(leak.data.mat), ncol = ncol(leak.data.mat))
for (i in 1:nrow(leak.data.mat)){
  leak.cum[i, ] <- cumsum(leak.data.mat[i, ])
}

start.time <- as_datetime("2024-02-12T00:00:00", tz = "America/Denver")
end.time <- as_datetime("2024-02-12T23:50:00", tz = "America/Denver")

this.int <- interval(start.time, end.time)

min.freq <- seq(start.time, end.time, by = "1 min")

leak.cum.min.freq <- matrix(nrow = length(min.freq), ncol = ncol(leak.cum))

leak.data.subset <- leak.data[leak.data$start %within% this.int, ]
leak.cum.subset <- leak.cum[leak.data$start %within% this.int, ]

for (i in 1:nrow(leak.cum.min.freq)){
  
  print(paste0(i, "/", nrow(leak.cum.min.freq)))
  
  for (j in 1:nrow(leak.data.subset)){
    
    this.leak.int <- interval(leak.data.subset$start[j],
                              leak.data.subset$end[j])
    
    if (min.freq[i] %within% this.leak.int){
      
      leak.cum.min.freq[i, ] <- leak.cum.subset[j, ]
    }
  }
}

leak.cum.min.freq[is.na(leak.cum.min.freq)] <- 0

high.res.mask <- orig.data$times %within% this.int
low.res.mask <- interval.times %within% this.int

png('../figures/data_example.png',
    res = 100, pointsize = 24, width = 1920, height = 1080*0.75)

par(mfcol = c(2,2))
par(mar = c(0.5,3.5,1.5,0.5))
par(oma = c(1.5,0,0,0))
par(mgp = c(2, 0.75, 0))

ylim.max <- 8

truth.alpha.val <- 0.8
plot(min.freq, leak.cum.min.freq[,1], type = "l", ylim = c(0,ylim.max), col = NA,
     ylab = "Emission rate [kg/hr]", yaxt= "n", xpd = NA, xaxt = "n", xlab = "")
envelopePlot(x1 = min.freq, y1 = rep(0, length(min.freq)), y2 = leak.cum.min.freq[, 1],
             lineCol = NA, col = alpha(cols[1], truth.alpha.val))
envelopePlot(x1 = min.freq, y1 = leak.cum.min.freq[, 1], y2 = leak.cum.min.freq[, 2],
             lineCol = NA, col = alpha(cols[2], truth.alpha.val))
envelopePlot(x1 = min.freq, y1 = leak.cum.min.freq[, 2], y2 = leak.cum.min.freq[, 3],
             lineCol = NA, col = alpha(cols[3], truth.alpha.val))
envelopePlot(x1 = min.freq, y1 = leak.cum.min.freq[, 3], y2 = leak.cum.min.freq[, 4],
             lineCol = NA, col = alpha(cols[4], truth.alpha.val))
envelopePlot(x1 = min.freq, y1 = leak.cum.min.freq[, 4], y2 = leak.cum.min.freq[, 5],
             lineCol = NA, col = alpha(cols[5], truth.alpha.val))
lines(min.freq, leak.cum.min.freq[,5], lwd = 3)

axis(side = 2, at = seq(0,8, by = 2))

lwd.val <- 2

plot(times[high.res.mask], apply(orig.data$obs[high.res.mask, ], 1, max), type = "l",
     ylab = "Methane concentration [ppm]", xaxt = "n", lwd = lwd.val,
     xpd = NA, xlab = "", yaxt = "n")
axis(side = 2, at = seq(0,60, by= 15))

date.seq <- seq(start.time, end.time+hours(4), by = "4 hours")
axis(side = 1, at = date.seq, labels = paste0(hour(date.seq), ":00"))

plot(times[high.res.mask], orig.data$WS[high.res.mask], type = "l",
     ylab = "Wind speed [m/s]", xaxt = "n", lwd = lwd.val, xpd = NA, xlab = "")

plot(times[high.res.mask], orig.data$WD[high.res.mask], type = "l",
     ylab = "Wind direction", lwd = lwd.val, xpd = NA,
     yaxt = "n", xaxt = "n")

axis(side = 2, at = c(0, pi/2, pi, 3*pi/2, 2*pi),
     labels = c("W", "S", "E", "N", "W"))

date.seq <- seq(start.time, end.time+hours(4), by = "4 hours")
axis(side = 1, at = date.seq, labels = paste0(hour(date.seq), ":00"))

dev.off()



png('../figures/data_example_full.png',
    res = 100, pointsize = 24, width = 1920/1.25, height = 1080*1.5)

par(mfrow = c(11,1))
par(mar = c(0.25,1.75,0.25,1.75))
par(oma = c(1.5,0,0,0))
par(mgp = c(2, 0.75, 0))

ylim.max <- 8

truth.alpha.val <- 0.8
plot(min.freq, leak.cum.min.freq[,1], type = "l", ylim = c(0,ylim.max), col = NA,
     ylab = "", yaxt= "n", xpd = NA, xaxt = "n", xlab = "")
envelopePlot(x1 = min.freq, y1 = rep(0, length(min.freq)), y2 = leak.cum.min.freq[, 1],
             lineCol = NA, col = alpha(cols[1], truth.alpha.val))
envelopePlot(x1 = min.freq, y1 = leak.cum.min.freq[, 1], y2 = leak.cum.min.freq[, 2],
             lineCol = NA, col = alpha(cols[2], truth.alpha.val))
envelopePlot(x1 = min.freq, y1 = leak.cum.min.freq[, 2], y2 = leak.cum.min.freq[, 3],
             lineCol = NA, col = alpha(cols[3], truth.alpha.val))
envelopePlot(x1 = min.freq, y1 = leak.cum.min.freq[, 3], y2 = leak.cum.min.freq[, 4],
             lineCol = NA, col = alpha(cols[4], truth.alpha.val))
envelopePlot(x1 = min.freq, y1 = leak.cum.min.freq[, 4], y2 = leak.cum.min.freq[, 5],
             lineCol = NA, col = alpha(cols[5], truth.alpha.val))
lines(min.freq, leak.cum.min.freq[,5], lwd = 3)

axis(side = 2, at = seq(0,8, by = 2))

mtext("True emission state [kg/hr]", side =3 , line = -1.15, cex = 0.65,
      col = "black", adj = 0.025)

lwd.val <- 2

tmp.obs <- orig.data$obs[high.res.mask, ]
colnames(tmp.obs) <- c("East", "East-Central", "Northeast", "North", "South", "Southeast", "Southwest", "West", "West-Central", "Northwest")

for (i in 1:ncol(tmp.obs)){
  plot(times[high.res.mask], tmp.obs[,i], type = "l", ylim = c(0,65),
       ylab = "", xaxt = "n", lwd = lwd.val,
       xpd = NA, xlab = "", yaxt = "n")
  axis(side = 4, at = seq(0,60, by= 15))
  
  mtext("Methane concentration\n[ppm]", side =3 , line = -2, cex = 0.65,
        col = "black", adj = 0.985)
  mtext(paste0("Sensor: ", colnames(tmp.obs[i])), side =3 , line = -1.15, cex = 0.65,
        col = "black")
}

date.seq <- seq(start.time, end.time+hours(4), by = "4 hours")
axis(side = 1, at = date.seq, labels = paste0(hour(date.seq), ":00"))

dev.off()





# STEP 8: CREATE INVENTORY AND ALERT RESULTS FIGURE ON ENTIRE DATASET
#---------------------------------------------------------------------------

# Take samples from the posterior of the betas to create distribution of inventory estimates
n.samples <- 4000
samples <- array(NA, dim = c(nrow(q.hat), ncol(q.hat), n.samples))

for (i in 1:ncol(q.hat)){
  for (j in 1:nrow(q.hat)){
    if (is.na(q.hat[j,i])){
      samples[j,i, ] <- sample(na.omit(q.hat[,i]), n.samples, replace = T)
    } else if (q.hat[j,i] == 0){
      samples[j,i, ] <- 0
    } else {
      samples[j,i, ] <- sample(mcmc.correct.dim[[j]][,i], n.samples, replace = T)
    }
  }
}

# Sum over the inversion windows
inventory.samples <- t(apply(samples, c(2,3), sum))*step.size/60/1000

equip.sum <- apply(inventory.samples, 2, mean)
equip.sum.lower <- apply(inventory.samples, 2, function(X) quantile(X, probs = 0.025))
equip.sum.upper <- apply(inventory.samples, 2, function(X) quantile(X, probs = 0.975))

site.total <- sum(equip.sum)
site.total.lower <- sum(equip.sum.lower)
site.total.upper <- sum(equip.sum.upper)

# OLD VERSION THAT IS NOT BASED ON SAMPLES
if (F){
  
  to.remove <- rep(F, length(q.hat.total))
  
  q.hat.avg.outlier.removed <- q.hat.avg[!to.remove,]*step.size/60/1000
  q.hat.lower.avg.outlier.removed <- q.hat.lower.avg[!to.remove,]*step.size/60/1000
  q.hat.upper.avg.outlier.removed <- q.hat.upper.avg[!to.remove,]*step.size/60/1000
  
  for (i in 1:length(source.names)){
    to.add.avg <- is.na(q.hat.avg.outlier.removed[,i])
    q.hat.avg.outlier.removed[to.add.avg, i] <- mean(q.hat.avg.outlier.removed[, i], na.rm = T)
    
    to.add.avg <- is.na(q.hat.lower.avg.outlier.removed[,i])
    q.hat.lower.avg.outlier.removed[to.add.avg, i] <- mean(q.hat.lower.avg.outlier.removed[, i], na.rm = T)
    
    to.add.avg <- is.na(q.hat.upper.avg.outlier.removed[,i])
    q.hat.upper.avg.outlier.removed[to.add.avg, i] <- mean(q.hat.upper.avg.outlier.removed[, i], na.rm = T)
  }
  
  equip.sum <- apply(q.hat.avg.outlier.removed, 2, sum)
  equip.sum.lower <- apply(q.hat.lower.avg.outlier.removed, 2, sum)
  equip.sum.upper <- apply(q.hat.upper.avg.outlier.removed, 2, sum)
  
  site.total <- sum(equip.sum)
  site.total.lower <- sum(equip.sum.lower)
  site.total.upper <- sum(equip.sum.upper)
  
}

equip.sum.truth <- apply(truth, 2, sum)*step.size/60/1000
site.total.truth <- sum(equip.sum.truth)

to.plot <- c(equip.sum, site.total)
to.plot.truth <- c(equip.sum.truth, site.total.truth)
to.plot.lower <- c(equip.sum.lower, site.total.lower)
to.plot.upper <- c(equip.sum.upper, site.total.upper)


png('../figures/loc_quant_results.png',
    res = 100, pointsize = 30, width = 1920, height = 1080*0.625)

par(mfrow = c(1,3))

par(mar = c(3,3,2,1))

par(mfrow = c(1,3))
par(mgp = c(2, 0.75, 0))

this.order <- order(to.plot, decreasing = T)

alpha.val <- 0.5
bar.cols <- c(wellhead.west.color,
              separator.west.color,
              tank.color, 
              wellhead.east.color, 
              separator.east.color,
              "black")

total.error.abs <- to.plot - to.plot.truth
total.error <- round(100* (to.plot - to.plot.truth) / to.plot.truth, 1)

b <- barplot(to.plot[this.order], col = alpha(bar.cols[this.order], alpha.val), 
             ylim = c(0,3),
             border = NA, ylab = "Total emissions [metric tons]")

adj.val <- 0.5
segments(x0 = b-adj.val, x1 = b+adj.val, y0 = to.plot.truth[this.order],
         col = bar.cols[this.order], lwd = 6)

segments(x0 = b, y0 = to.plot.lower[this.order], y1 = to.plot.upper[this.order], lwd = 2)

legend("right", c("West Wellhead", "West Separator", "Tank", "East Wellhead", "East Separator", "Site Total")[this.order],
       fill = bar.cols[this.order], box.lwd = NA)

legend("topright", c("Truth (line)", "Estimate (box)"),
       lwd = c(6, NA), fill = c(NA, alpha("black", alpha.val)),
       box.col = "white", border = c("white", "black"))

text(x = b-0.25,
     y = -0.2,
     labels = trimws(paste0(format(total.error[this.order], nsmall = 1), "%")),
     offset = 5,
     srt = 25,
     cex = 0.95,
     xpd = NA)

line.col <- "steelblue4"

error <- q.hat.total - truth.total
hist(error, xlim = c(-6,6), breaks = seq(-999,999, by = 0.5), xaxt = "n",
     yaxt = "n",
     ylab = "Frequency [thousands]",
     xlab = "",
     ylim= c(0,1500))
axis(side = 1, at = seq(-6,6, by = 2))
segments(x0 = mean(error, na.rm = T), y0 = -999, y1 = 99999, lwd = 4, col = line.col)
segments(x0 = quantile(error, probs = c(0.025, 0.975), na.rm = T), 
         y0 = -999, y1 = 99999, lwd = 4, col = line.col, lty = 2)
axis(side = 2, at = seq(0,1500, by = 500), labels = seq(0,1.5, by = 0.5))



# Determine true emission state (on / off)
is.emitting <- truth > 0

# Determine estimated emission state (on / off)
est.emitting <- matrix(NA, nrow = nrow(is.emitting), ncol = ncol(is.emitting))
for (i in 1:5){
  est.emitting[,i] <- pis[,i] > mean(pis[,i], na.rm = T)
}
est.emitting[q.hat == 0] <- F

# Figure out the number of correct localization estimates per release
num.correct <- matrix(NA, nrow = nrow(is.emitting), ncol = length(source.names))
for (i in 1:nrow(is.emitting)){
  num.correct[i,] <- (is.emitting[i, ] & est.emitting[i,]) | (!is.emitting[i, ] & !est.emitting[i,])
}

# Format as table and percent
num.correct.vec <- apply(num.correct, 1, sum)

percent <- table(num.correct.vec) / sum(!is.na(num.correct.vec))
to.plot <- table(num.correct.vec)

b <- barplot(to.plot,  yaxt = "n",
             col = alpha(mako(7)[1:6], 1),
             border = NA,
             ylim = c(0,1500),
             main = round(mean(num.correct.vec, na.rm = T), 2),
             ylab = "Frequency [thousands]")

axis(side = 1, at = b, labels = NA)
axis(side = 2, at = seq(0,1500, by = 500), labels = seq(0,1.5, by = 0.5))
text(x = b,
     y = to.plot+50,
     labels = paste0(round(100*percent, 1), "%"),
     offset = 5,
     xpd = NA)


dev.off()



source.names
round(total.error.abs, 1)
total.error
round(mean(error, na.rm = T), 2)
round(IQR(error, na.rm = T), 2)
round(coverage, 2)
round(mean(num.correct.vec, na.rm = T), 2)

for (i in 1:ncol(est.emitting)){
  
  source.names[i]
  
  to.remove <- is.na(est.emitting[,i])
  est.emitting2 <- est.emitting[!to.remove, ]
  is.emitting2 <- is.emitting[!to.remove, ]
  
  sum(is.emitting2[,i])
  sum(!is.emitting2[,i])
  sum(est.emitting2[,i])
  sum(!est.emitting2[,i])
  
  tp = sum(is.emitting2[,i] & est.emitting2[,i]) 
  fp = sum(!is.emitting2[,i] & est.emitting2[,i]) 
  
  fn = sum(is.emitting2[,i] & !est.emitting2[,i]) 
  tn = sum(!is.emitting2[,i] & !est.emitting2[,i]) 
  
  tpr = tp / (tp + fn)
  tnr = tn / (fp + tn)
  
  ppv = tp / (tp + fp)
  npv = tn / (tn + fn)
  
  accuracy = (tp + tn) / (tp + fp + fn + tn)
  
  matrix(c(source.names[i],         sum(is.emitting2[,i]),                                     sum(!is.emitting2[,i]),                                    NA, 
           sum(est.emitting2[,i]),  paste0(tp, " (", 100*round(tp/nrow(is.emitting2),3), ")"), paste0(fp, " (", 100*round(fp/nrow(is.emitting2),3), ")"), round(100*ppv,1), 
           sum(!est.emitting2[,i]), paste0(fn, " (", 100*round(fn/nrow(is.emitting2),3), ")"), paste0(tn, " (", 100*round(tn/nrow(is.emitting2),3), ")"), round(100*npv,1), 
           NA,                      round(100*tpr, 1),                                         round(100*tnr,1),                                          round(100*accuracy, 1)),
         nrow = 4, byrow = T)
  
}



site.level.est.emitting <- apply(est.emitting, 1, any, na.rm = T)
site.level.is.emitting <- apply(is.emitting, 1, any, na.rm = T)

sum(site.level.is.emitting)
sum(!site.level.is.emitting)
sum(site.level.est.emitting)
sum(!site.level.est.emitting)

tp = sum(site.level.is.emitting & site.level.est.emitting) 
fp = sum(!site.level.is.emitting & site.level.est.emitting) 

fn = sum(site.level.is.emitting & !site.level.est.emitting) 
tn = sum(!site.level.is.emitting & !site.level.est.emitting) 

tpr = tp / (tp + fn)
tnr = tn / (fp + tn)

ppv = tp / (tp + fp)
npv = tn / (tn + fn)

accuracy = (tp + tn) / (tp + fp + fn + tn)

matrix(c("site-level",                  sum(site.level.is.emitting),                                           sum(!site.level.is.emitting),                                          NA, 
         sum(site.level.est.emitting),  paste0(tp, " (", 100*round(tp/length(site.level.is.emitting),3), ")"), paste0(fp, " (", 100*round(fp/length(site.level.is.emitting),3), ")"), round(100*ppv,1), 
         sum(!site.level.est.emitting), paste0(fn, " (", 100*round(fn/length(site.level.is.emitting),3), ")"), paste0(tn, " (", 100*round(tn/length(site.level.is.emitting),3), ")"), round(100*npv,1), 
         NA,                            round(100*tpr, 1),                                                     round(100*tnr,1),                                                      round(100*accuracy, 1)),
       nrow = 4, byrow = T)



# STEP 5: CREATE LOCALIZATION AND QUANTIFICATION RESULTS FIGURE
#---------------------------------------------------------------------------

to.remove <- !info.to.use

mcmc.correct.dim <- mcmc.correct.dim[!to.remove]

q.hat.tmp <- q.hat[!to.remove, ]

n.samples <- 4000
samples <- array(NA, dim = c(nrow(q.hat.tmp), ncol(q.hat.tmp), n.samples))

for (i in 1:ncol(q.hat.tmp)){
  for (j in 1:nrow(q.hat.tmp)){
    if (is.na(q.hat.tmp[j,i])){
      samples[j,i, ] <- sample(na.omit(q.hat.tmp[,i]), n.samples, replace = T)
    } else if (q.hat.tmp[j,i] == 0){
      samples[j,i, ] <- 0
    } else {
      samples[j,i, ] <- sample(mcmc.correct.dim[[j]][,i], n.samples, replace = T)
    }
  }
}

inventory.samples <- t(apply(samples, c(2,3), sum))*step.size/60/1000

equip.sum <- apply(inventory.samples, 2, mean)
equip.sum.lower <- apply(inventory.samples, 2, function(X) quantile(X, probs = 0.025))
equip.sum.upper <- apply(inventory.samples, 2, function(X) quantile(X, probs = 0.975))

# OLD VERSION THAT IS NOT BASED ON SAMPLING
if (F){
  
  q.hat.avg.outlier.removed <- q.hat.avg[!to.remove,]*step.size/60/1000
  q.hat.lower.avg.outlier.removed <- q.hat.lower.avg[!to.remove,]*step.size/60/1000
  q.hat.upper.avg.outlier.removed <- q.hat.upper.avg[!to.remove,]*step.size/60/1000
  
  for (i in 1:length(source.names)){
    to.add.avg <- is.na(q.hat.avg.outlier.removed[,i])
    q.hat.avg.outlier.removed[to.add.avg, i] <- mean(q.hat.avg.outlier.removed[, i], na.rm = T)
    
    to.add.avg <- is.na(q.hat.lower.avg.outlier.removed[,i])
    q.hat.lower.avg.outlier.removed[to.add.avg, i] <- mean(q.hat.lower.avg.outlier.removed[, i], na.rm = T)
    
    to.add.avg <- is.na(q.hat.upper.avg.outlier.removed[,i])
    q.hat.upper.avg.outlier.removed[to.add.avg, i] <- mean(q.hat.upper.avg.outlier.removed[, i], na.rm = T)
  }
  
  equip.sum <- apply(q.hat.avg.outlier.removed, 2, sum)
  equip.sum.lower <- apply(q.hat.lower.avg.outlier.removed, 2, sum)
  equip.sum.upper <- apply(q.hat.upper.avg.outlier.removed, 2, sum)
  
}

q.hat.total.info.filtered <- apply(q.hat.avg[!to.remove, ], 1, sum)
truth.total.info.filtered <- apply(truth[!to.remove, ], 1, sum)

site.total <- sum(equip.sum)
site.total.lower <- sum(equip.sum.lower)
site.total.upper <- sum(equip.sum.upper)

equip.sum.truth <- apply(truth[!to.remove, ], 2, sum)*step.size/60/1000
site.total.truth <- sum(equip.sum.truth)

to.plot <- c(equip.sum, site.total)
to.plot.truth <- c(equip.sum.truth, site.total.truth)
to.plot.lower <- c(equip.sum.lower, site.total.lower)
to.plot.upper <- c(equip.sum.upper, site.total.upper)


png('../figures/loc_quant_results_info_adjusted.png',
    res = 100, pointsize = 30, width = 1920, height = 1080*0.625)

par(mfrow = c(1,3))

par(mar = c(3,3,2,1))

par(mfrow = c(1,3))
par(mgp = c(2, 0.75, 0))

this.order <- order(to.plot, decreasing = T)

alpha.val <- 0.5
bar.cols <- c(wellhead.west.color,
              separator.west.color,
              tank.color, 
              wellhead.east.color, 
              separator.east.color,
              "black")

total.error <- round(100* (to.plot - to.plot.truth) / to.plot.truth, 1)

b <- barplot(to.plot[this.order], col = alpha(bar.cols[this.order], alpha.val), 
             ylim = c(0,0.3),
             border = NA, ylab = "Total emissions [metric tons]")

adj.val <- 0.5
segments(x0 = b-adj.val, x1 = b+adj.val, y0 = to.plot.truth[this.order],
         col = bar.cols[this.order], lwd = 6)

segments(x0 = b, y0 = to.plot.lower[this.order], y1 = to.plot.upper[this.order], lwd = 2)


legend("right", c("West Wellhead", "West Separator", "Tank", "East Wellhead", "East Separator", "Site Total")[this.order],
       fill = bar.cols[this.order], box.lwd = NA)

legend("topright", c("Truth (line)", "Estimate (box)"),
       lwd = c(6, NA), fill = c(NA, alpha("black", alpha.val)),
       box.col = "white", border = c("white", "black"))

text(x = b-0.25,
     y = -0.02,
     labels = trimws(paste0(format(total.error[this.order], nsmall = 1), "%")),
     offset = 5,
     srt = 25,
     cex = 0.95,
     xpd = NA)

line.col <- "steelblue4"

error <- q.hat.total.info.filtered - truth.total.info.filtered
hist(error, xlim = c(-6,6), breaks = seq(-999,999, by = 0.5), xaxt = "n",
     yaxt = "n",
     ylab = "Frequency",
     xlab = "",
     ylim= c(0,200))
axis(side = 1, at = seq(-6,6, by = 2))
segments(x0 = mean(error, na.rm = T), y0 = -999, y1 = 99999, lwd = 4, col = line.col)
segments(x0 = quantile(error, probs = c(0.025, 0.975), na.rm = T), 
         y0 = -999, y1 = 99999, lwd = 4, col = line.col, lty = 2)
axis(side = 2, at = seq(0,200, by = 50))


# Determine true emission state (on / off)
is.emitting <- truth > 0

# Determine estimated emission state (on / off)
est.emitting <- matrix(NA, nrow = nrow(is.emitting), ncol = ncol(is.emitting))
for (i in 1:5){
  est.emitting[,i] <- pis[,i] > mean(pis[,i], na.rm = T)
}
est.emitting[q.hat == 0] <- F

# Figure out the number of correct localization estimates per release
num.correct <- matrix(NA, nrow = nrow(is.emitting), ncol = length(source.names))
for (i in 1:nrow(is.emitting)){
  num.correct[i,] <- (is.emitting[i, ] & est.emitting[i,]) | (!is.emitting[i, ] & !est.emitting[i,])
}


# Format as table and percent
num.correct.vec <- apply(num.correct, 1, sum)
num.correct.vec <- num.correct.vec[!to.remove]

percent <- table(num.correct.vec) / sum(!is.na(num.correct.vec))
to.plot <- table(num.correct.vec)

b <- barplot(to.plot,  yaxt = "n",
             col = alpha(mako(7)[1:6], 1),
             border = NA,
             ylim = c(0,200),
             main = round(mean(num.correct.vec, na.rm = T), 2),
             ylab = "Frequency")

axis(side = 1, at = b, labels = NA)
axis(side = 2, at = seq(0,200, by = 50))
text(x = b,
     y = to.plot+5,
     labels = paste0(round(100*percent, 1), "%"),
     offset = 5,
     xpd = NA)

dev.off()





# STEP 3: COMPUTE NUMBER OF ESTIMATES NEEDED FOR AVERAGE TO CONVERGE
#---------------------------------------------------------------------------

# Number of MC samples
n.samples <- 5000

# Variable to hold average error across MC samples
avg.error <- matrix(nrow = n.samples, ncol = length(q.hat.total))

# Loop through number of emission events to average
for (i in 1:1000){
  
  print(paste0(i, "/", length(q.hat.total)))
  
  # Loop through MC samples
  for (j in 1:n.samples){
    
    # Sample in events to average
    these.ind <- sample.int(length(q.hat.total), size = i)
    
    # Compute error for each of the sampled events
    this.error <- q.hat.total[these.ind] - truth.total[these.ind]
    
    # Compute average error
    avg.error[j,i] <- mean(this.error, na.rm =T)
    
  }
}


png('../figures/num_to_average.png',
    res = 100, pointsize = 24, width =1920, height = 1080)

par(mfrow = c(1,1))

par(mar = c(4,4,2,2))
par(mgp = c(2.5,1,0))

lwd.val <- 3
cols <- c("#309E78", "#DA6001", "#7570B4", "#E82B8A")

err.quantile <- apply(avg.error, 2, function(X) quantile(X, probs = c(0,0.025, 0.975, 1), na.rm = T))

plot(apply(avg.error, 2, mean), yli = c(-3, 3), col = "white", ylab = "Average error [kg/hr]",
     xlab = "Number of 30-minute inversion windows included in average",
     xlim = c(0,600))

envelopePlot(x1 = 1:1000, y1 = err.quantile[2, 1:1000], y2 = err.quantile[3, 1:1000],
             col = alpha(cols[3], 0.25), lineCol = cols[3])

lines(apply(avg.error, 2, mean), lwd = lwd.val + 1)

lines(1:length(q.hat.total), err.quantile[1,], col = "gray55", lwd = lwd.val)
lines(1:length(q.hat.total), err.quantile[4,], col = "gray55", lwd = lwd.val)

num.to.avg <- max(c(min(which(err.quantile[2, ] > -1)), min(which(err.quantile[3, ] < 1))))

abline(h = c(-1,1), col = cols[2], lty = 2, lwd= lwd.val+1)

abline(v = num.to.avg, col = cols[1],lty = 2, lwd = lwd.val+1)

legend("topright", c("Average of average errors", "Inner 95% of average errors", "Min / max of average errors", "95% probability of average error being within [-1,1] kg/hr"),
       col = c('black', cols[3], "gray55", cols[1]),
       lty= c(1,1,1,2), bty = "n",
       lwd = 5)

dev.off()


