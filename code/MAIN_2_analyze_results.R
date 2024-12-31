# Description: 
# Summarizes event detection, localization, and quantification results and 
# reproduces paper figures.
# Author: William Daniels (wdaniels@mines.edu)
# Last Updated: December 2024

# Clear environment
if(!is.null(dev.list())){dev.off()}
rm(list = ls())
gc()

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

# Set to true if analyzing results using Laplace errors.
# Necessary because the file format of the Laplace error model is different than the 
# normal MDLQ output.
laplace.results <- F

# Read in MDLQ output
if (laplace.results){
  data <- readRDS('../output_data/MDLQ_output_ADED2024_laplace_errors.RData')
} else {
  data <- readRDS('../output_data/MDLQ_output_ADED2024.RData')  
}

# END OF USER INPUT - NO MODIFICATION NECESSARY BELOW THIS POINT
#---------------------------------------------------------------------------





# STEP 1: PARSE OUT RESULTS
#---------------------------------------------------------------------------

# Parse out source names
source.names <- data$source

# Parse out results. Different file format if using laplace errors.
if (laplace.results){
  
  # Get emission rate estimates and 95% credible intervals
  q.hat <- data$q.hat
  q.hat.lower <- data$q.hat.lower
  q.hat.upper <- data$q.hat.upper
  
  # Get emission indicator estimates
  pis <- data$pi
  
  # Get error variance estimates
  sigma2 <- data$sigma2
  
} else {
  
  # Parse out MDLQ output from other saved variables
  out <- data$out
  
  # Initialize variables to hold condensed MDLQ results
  q.hat <- q.hat.lower <- q.hat.upper <- pis <- matrix(NA, nrow = length(out), ncol = 5)
  sigma2 <- vector(length = length(out))
  
  # Loop through emission events
  for (i in 1:length(out)){
    
    # Parse out all of the 30-minute inversion intervals within this event
    these.mcmc <- out[[i]][[1]]
    
    # Keep only the 30-minute intervals where the MDLQ converged
    to.use <- sapply(these.mcmc, length) > 1
    these.mcmc <- these.mcmc[to.use]
    
    # Skip this emission event if no 30-minute intervals converged
    if (length(these.mcmc) == 0){ next }
    
    # Get emission rate estimates from each 30-minute interval within this event
    all.betas <- do.call(rbind, lapply(these.mcmc, function(X) apply(X, 2, mean)[c(16:20)]))
    
    # Remove 30-minute intervals with outlying beta estimates
    to.remove <- rep(F, nrow(all.betas))
    if (nrow(all.betas) > 3){
      to.remove <- matrix(nrow = nrow(all.betas), ncol = ncol(all.betas))
      for (p in 1:ncol(all.betas)){
        these.est <- all.betas[,p]
        for (r in 1:nrow(all.betas)){
          to.remove[r,p] <- these.est[r] / mean(these.est[-seq(max(1,r-1), min(nrow(all.betas), r+1))]) > 10^3
        }
      }
      to.remove <- apply(to.remove, 1, any, na.rm = T)
    }
    
    # Average emission rate estimates across 30-minute intervals within this event
    q.hat[i, ] <- apply(matrix(all.betas[!to.remove, ], ncol = 5), 2, mean) * 3.6
    
    # Average emission rate estimate lower bounds across 30-minute intervals within this event
    all.lower <- do.call(rbind, lapply(these.mcmc, function(X) apply(X, 2, function(X) quantile(X, probs = 0.025))[c(16:20)]))
    q.hat.lower[i, ] <- apply(matrix(all.lower[!to.remove], ncol = 5), 2, mean) * 3.6
    
    # Average emission rate estimate upper bounds across 30-minute intervals within this event
    all.upper <- do.call(rbind, lapply(these.mcmc, function(X) apply(X, 2, function(X) quantile(X, probs = 0.975))[c(16:20)]))
    q.hat.upper[i, ] <- apply(matrix(all.upper[!to.remove], ncol = 5), 2, mean) * 3.6
    
    # Average error variance estimates across 30-minute intervals within this event
    all.sigma2 <- do.call(rbind, lapply(these.mcmc, function(X) apply(X, 2, mean)[c(21)]))
    sigma2[i] <- apply(all.sigma2, 2, mean)
    
    # Average emission indicator estimates across 30-minute intervals within this event
    all.pis <- do.call(rbind, lapply(these.mcmc, function(X) apply(X, 2, mean)[c(1:5)]))
    pis[i, ] <- apply(matrix(all.pis, ncol = 5), 2, mean)
    
  }
  
}

# STEP 2: READ IN CONTROLLED RELEASE DATA
#---------------------------------------------------------------------------

# Read in controlled release data 
leak.data <- readRDS('/Users/wdaniels/Documents/papers/MDLQ/input_data/leak_data_clean.RData')

# NA's mean no emissions, so replace with a 0 kg/hr value
leak.data[is.na(leak.data)] <- 0

# Reorganize columns
leak.data <- leak.data[,c(1,2,3,4,6,7,5)]

# Get controlled release start and end times
leak.data$start <- as_datetime(leak.data$start, tz = "America/Denver")
leak.data$end   <- as_datetime(leak.data$end, tz = "America/Denver")

# Get controlled release durations
durations <- as.numeric(difftime(leak.data$end, leak.data$start, units = "min"))

# Filter out events less than 30 minutes
to.keep <- leak.data$end < max(data$times) & durations >= 30
leak.data <- leak.data[to.keep, ]
durations <- durations[to.keep]

# Get site-level emission rate estimates and 95% credible intervals
q.hat.total <- apply(q.hat, 1, sum)
q.hat.total.lower <- apply(q.hat.lower, 1, sum)
q.hat.total.upper <- apply(q.hat.upper, 1, sum)

# Get true site-level release rates
leak.total <- apply(leak.data[3:7], 1, sum)
leak.data.mat <- as.matrix(leak.data[,c(3:7)])

# Compute site-level emission rate coverage
coverage <- sum(leak.total >= q.hat.total.lower & leak.total <= q.hat.total.upper, na.rm = T) / sum(!is.na(q.hat.total))




# STEP 3: COMPUTE NUMBER OF ESTIMATES NEEDED FOR AVERAGE TO CONVERGE
#---------------------------------------------------------------------------

# Number of MC samples
n.samples <- 5000

# Variable to hold average error across MC samples
avg.error <- matrix(nrow = n.samples, ncol = length(q.hat.total))

# Loop through number of emission events to average
for (i in 1:length(q.hat.total)){
  
  print(paste0(i, "/", length(q.hat.total)))
  
  # Loop through MC samples
  for (j in 1:n.samples){
    
    # Sample in events to average
    these.ind <- sample.int(length(q.hat.total), size = i)
    
    # Compute error for each of the sampled events
    this.error <- q.hat.total[these.ind] - leak.total[these.ind]
    
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
     xlab = "Number of Emission Events Included in Average")

envelopePlot(x1 = 1:337, y1 = err.quantile[2, ], y2 = err.quantile[3,],
             col = alpha(cols[3], 0.25), lineCol = cols[3])

lines(apply(avg.error, 2, mean), lwd = lwd.val + 1)

lines(1:337, err.quantile[1,], col = "gray55", lwd = lwd.val)
lines(1:337, err.quantile[4,], col = "gray55", lwd = lwd.val)

num.to.avg <- max(c(min(which(err.quantile[2, ] > -1)), min(which(err.quantile[3, ] < 1))))

abline(h = c(-1,1), col = cols[2], lty = 2, lwd= lwd.val+1)

abline(v = num.to.avg, col = cols[1],lty = 2, lwd = lwd.val+1)

legend("topright", c("Average of average errors", "Inner 95% of average errors", "Min / max of average errors", "95% probability of average error being within [-1,1] kg/hr"),
       col = c('black', cols[3], "gray55", cols[1]),
       lty= c(1,1,1,2), bty = "n",
       lwd = 5)

dev.off()





# STEP 4: COMPUTE CUMULATIVE EMISSIONS ERROR
#---------------------------------------------------------------------------

# Identify NA emission rate estimates. These correspond to releases where the MDLQ did not converge.
to.remove <- is.na(q.hat.total)

# First compute cumulative error when simply removing these NA events.
q.hat.total.to.compare <- q.hat.total[!to.remove]
leak.total.to.compare <- leak.total[!to.remove]
leak.times.to.compare <- leak.data$start[!to.remove]

# Get total emissions [kg] by multiplying rate by duration of the release in hours
q.hat.total.to.compare <- q.hat.total.to.compare * durations[!to.remove]/60
leak.total.to.compare <- leak.total.to.compare * durations[!to.remove]/60

sum(q.hat.total.to.compare)
sum(leak.total.to.compare)

png('../figures/cumulative_errors.png',
    res = 100, pointsize = 24, width =1920, height = 1080*0.8)

lwd.val <- 5

est.col <- "gray55"
error.col <- mako(10)[7]

par(mfrow = c(1,1))
par(mar = c(2,4,1,1))
par(mgp = c(3, 1, 0))

plot(leak.times.to.compare, cumsum(leak.total.to.compare), type = "l", lwd = lwd.val,
     ylab = "Cumulative Emissions [kg]", xlab = "", ylim = c(0,2500),
     xlim = c(as_datetime("2024-02-1"), as_datetime("2024-5-1")))

abline(h = 0, lwd = lwd.val, lty = 2)

lines(leak.times.to.compare, cumsum(q.hat.total.to.compare), col = est.col, lwd = lwd.val, lty = 1)

cumsum.diff <- cumsum(q.hat.total.to.compare) - cumsum(leak.total.to.compare)

lines(leak.times.to.compare, cumsum.diff, col = error.col, lwd = lwd.val)

100 * cumsum.diff[length(cumsum.diff)] / cumsum(leak.total.to.compare)[length(leak.total.to.compare)]

legend("topleft", legend = c("True Emissions", "Estimated Emissions", "Error (Estimated - True)"),
       lty = c(1,1,1), col = c("black", est.col, error.col),
       lwd = 5)

dev.off()


# Identify NA emission rate estimates. These correspond to releases where the MDLQ did not converge.
to.remove <- is.na(q.hat.total)

# Next compute cumulative error when replacing the NA events with the average emission rate estimate across events.
# Compute average
avg.est <- mean(q.hat.total, na.rm = T)

# Fill in NA events with average, and get the true release rates
q.hat.total.to.compare <- q.hat.total
q.hat.total.to.compare[to.remove] <- avg.est
leak.total.to.compare <- leak.total
leak.times.to.compare <- leak.data$start

# Get total emissions [kg] by multiplying rate by duration of the release in hours
q.hat.total.to.compare <- q.hat.total.to.compare * durations/60
leak.total.to.compare <- leak.total.to.compare * durations/60

sum(q.hat.total.to.compare)
sum(leak.total.to.compare)

png('../figures/cumulative_errors_no_converge_fill_in.png',
    res = 100, pointsize = 24, width =1920, height = 1080*0.8)

lwd.val <- 5

est.col <- "gray55"
error.col <- mako(10)[7]

par(mfrow = c(1,1))
par(mar = c(2,4,1,1))
par(mgp = c(3, 1, 0))

plot(leak.times.to.compare, cumsum(leak.total.to.compare), type = "l", lwd = lwd.val,
     ylab = "Cumulative Emissions [kg]", xlab = "", ylim = c(0,2500),
     xlim = c(as_datetime("2024-02-1"), as_datetime("2024-5-1")))

abline(h = 0, lwd = lwd.val, lty = 2)

lines(leak.times.to.compare, cumsum(q.hat.total.to.compare), col = est.col, lwd = lwd.val, lty = 1)

cumsum.diff <- cumsum(q.hat.total.to.compare) - cumsum(leak.total.to.compare)

lines(leak.times.to.compare, cumsum.diff, col = error.col, lwd = lwd.val)

100 * cumsum.diff[length(cumsum.diff)] / cumsum(leak.total.to.compare)[length(leak.total.to.compare)]

legend("topleft", legend = c("True Emissions", "Estimated Emissions", "Error (Estimated - True)"),
       lty = c(1,1,1), col = c("black", est.col, error.col),
       lwd = 5)

dev.off()






# STEP 5: CREATE LOCALIZATION AND QUANTIFICATION RESULTS FIGURE
#---------------------------------------------------------------------------

# Determine true emission state (on / off)
is.emitting <- leak.data.mat > 0

# Determine estimated emission state (on / off)
est.emitting <- matrix(NA, nrow = nrow(leak.data.mat), ncol = ncol(leak.data.mat))
for (i in 1:5){
  est.emitting[,i] <- pis[,i] > mean(pis[,i], na.rm = T)
}

# Figure out the number of correct localization estimates per release
num.correct <- matrix(NA, nrow = nrow(is.emitting), ncol = 5)
for (i in 1:nrow(is.emitting)){
  num.correct[i,] <- (is.emitting[i, ] & est.emitting[i,]) | (!is.emitting[i, ] & !est.emitting[i,])
}

# Format as table and percent
num.correct.vec <- apply(num.correct, 1, sum)
percent <- table(num.correct.vec) / sum(!is.na(num.correct.vec))
to.plot <- table(num.correct.vec)

# For ADED 2024, can get anywhere from 0 to 5 localization estimates correct per release
# If no estimates in any of these categories (0 through 5), supplement vecotrs with a zero
if (length(to.plot) < 6){
  to.plot <- c(0,as.numeric(to.plot))
  percent <- c(0,as.numeric(percent))
}

png('../figures/controlled_release_results.png',
    res = 100, pointsize = 32, width = 1920, height = 1080*0.625)

par(mar = c(3,3,2,1))

par(mfrow = c(1,3))
par(mgp = c(2, 0.75, 0))

fit <- lm(q.hat.total ~ leak.total)
lim.max <- 10
plot(leak.total, q.hat.total, ylim = c(0,lim.max), xlim = c(0,lim.max), col = "white",
     main = "", xlab = "True emission rate [kg/hr]", asp = 1,
     ylab = "Estimated emission rate [kg/hr]")
mtext(side = 3,
      paste0("Best fit: y = ", round(coef(fit)[2], 2), "x + ", round(coef(fit)[1], 2)),
      line = 0.5)
abline(a = 0, b = 1, lwd = 2)
leak.data.mat <- as.matrix(leak.data[,c(3:7)])
envelopePlot(x1 = seq(0, 20, by= 1), y1 = seq(0, 20/2, by = 1/2), y2 = seq(0,20*2, by = 1*2),
             col = alpha("black", 0.2), lineCol = NA)
points(leak.total, q.hat.total, pch = 19, cex = 0.65, col = alpha("black", 0.5))
segments(x0 = leak.total, y0 = q.hat.total.lower, y1 = q.hat.total.upper, col = alpha("black", 0.2), lwd = 1.25)

abline(fit, col = "darkgoldenrod1", lwd = 4, lty = 1)

total.error <- q.hat.total - leak.total
hist(total.error, xlim = c(-8,8), breaks = seq(-999,999, by = 1), main = "", yaxt = "n", ylim = c(0,150),
     xlab = "Estimated rate - true rate [kg/hr]", xaxt = "n")
axis(side = 1, at = seq(-8,8, by = 2))
segments(x0 = mean(total.error, na.rm = T), y0 = -10, y1 = 140, lwd = 4, col = "brown")
segments(x0 = quantile(total.error, probs = c(0.025, 0.975), na.rm = T), y0 = -10, y1 = 140, lwd = 4, col = "brown", lty = 2)
mtext(side = 3, paste0("Site-level quantification errors\nAverage error = ", round(mean(total.error, na.rm = T),2), " kg/hr"), 
      line = -0.5)
axis(side = 2, at = seq(0,120, by = 30))

sum(total.error<0, na.rm = T)/sum(!is.na(total.error))

names(to.plot) <- 0:5
if (laplace.results){
b <- barplot(to.plot, col = (inferno(6)), yaxt = "n", 
             ylim = c(0,200),
             main = round(mean(num.correct.vec, na.rm = T), 2),
             ylab = "Frequency",
             xlab = "Correct localization estimates per release")
} else {
  b <- barplot(to.plot, col = (inferno(6)), yaxt = "n", 
               ylim = c(0,175),
               main = round(mean(num.correct.vec, na.rm = T), 2),
               ylab = "Frequency",
               xlab = "Correct localization estimates per release")
}

axis(side = 2, at = seq(0,150, by = 50))
text(x = b,
     y = to.plot+12,
     labels = paste0(round(100*percent, 1), "%"),
     offset = 5,
     xpd = NA)

dev.off()

