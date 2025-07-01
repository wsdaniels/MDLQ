rm(list = ls())

library(scales)
library(viridis)
library(lubridate)

data <- readRDS('../input_data/forward_model_output_ADED2024.RData')

truth <- readRDS('../output_data/truth_ADED2024_30min.RData')

# Set colors for plots
tank.color <- "#3062CF" #blue
wellhead.east.color <- "#9147B8" #purple
separator.east.color <- "#7EAD52" #green
separator.west.color <- "#F1C30E" #gold
wellhead.west.color <- "#C7383C" #red


cols <- c(wellhead.west.color, separator.west.color, tank.color, wellhead.east.color, separator.east.color)

cex.val <- 0.75
pt.cex.val <- 1.15 #1.25
lwd.val <- 3
horiz.lwd.val <- 2
length.val <- 0.15
adj.val <- 0
line.val <- 0.1

offset.val <- 0.3

n.exp <- 5
n.beta <- 5
sigma2.true <- 1

ac.vals <- c(0, 0.25, 0.5, 0.75, 0.95)


png('../figures/simulation_study_spike_misalignment.png',
    res = 100, pointsize = 32, width = 1920, height = 1080*1.2)

laytout.mat <- matrix(c(1, 5,
                        2, 6,
                        3, 7,
                        4, 8), nrow = 4, byrow = T)

layout(laytout.mat)

par(mar = c(1.5,2.5,0.75,0))

par(oma = c(1.5,0,0.5,0.25))

par(mgp = c(4, 0.75, 0))

offset <- seq(-offset.val, offset.val, length.out = n.beta+1)
offset[1:5] <- offset[1:5]-0.05
offset[6] <- offset[6] + 0.05

xlim.fix <- 0.2


results <- readRDS('./results/simulation_study_output_spike_misalignment.RData')

plot(1,1, col = "white", xlim = c(0.5+xlim.fix,n.exp+0.5-xlim.fix), ylim = c(0, 300), xaxt = "n", xlab = "",
     yaxt = "n")

axis(side = 2, at = seq(0, 300, length.out = 3))
axis(side = 2, at = seq(0, 300, length.out = 5), labels = NA)

mtext("Error variance [kg/hr]", side = 3, line = line.val, cex = cex.val, adj = adj.val)

rect(xleft  = seq(1,n.exp, by = 2) - 0.5,
     xright = seq(1,n.exp, by = 2) + 0.5,
     ybottom = -999, ytop = 999, 
     col = alpha("black", 0.2),
     border = NA)

abline(h = 2, lty = 2, col = "red", lwd = horiz.lwd.val)


for (i in 1:n.exp){
  
  arrows(x0 = i,
         y0 = quantile(results[[i]]$sigma2, probs = 0.025, na.rm = T),
         y1 = quantile(results[[i]]$sigma2, probs = 0.975, na.rm = T),
         col = "black", lwd = lwd.val,
         code = 3, angle = 90, length = length.val)
  
  points(i, mean(results[[i]]$sigma2, na.rm = T), pch = 19, col = "black", cex = pt.cex.val)
}


plot(1,1, col = "white", xlim = c(0.5+xlim.fix,n.exp+0.5-xlim.fix), ylim = c(0, 1), xaxt = "n", xlab = "",
     yaxt = "n")

axis(side = 2, at = seq(0, 1, length.out = 3))
axis(side = 2, at = seq(0, 1, length.out = 5), labels = NA)

mtext("Autocorrelation coefficient", side = 3, line = line.val, cex = cex.val, adj = adj.val)

rect(xleft  = seq(1,n.exp, by = 2) - 0.5,
     xright = seq(1,n.exp, by = 2) + 0.5,
     ybottom = -999, ytop = 999, 
     col = alpha("black", 0.2),
     border = NA)

abline(h = 0, lty = 2, col = "red", lwd = horiz.lwd.val)


for (i in 1:n.exp){
  
  
  arrows(x0 = i,
         y0 = quantile(results[[i]]$r, probs = 0.025, na.rm = T),
         y1 = quantile(results[[i]]$r, probs = 0.975, na.rm = T),
         col = "black", lwd = lwd.val,
         code = 3, angle = 90, length = length.val)
  
  points(i, mean(results[[i]]$r, na.rm = T), pch = 19, col = "black", cex = pt.cex.val)
}


plot(1,1, col = "white", xlim = c(0.5+xlim.fix,n.exp+0.5-xlim.fix), ylim = c(-2,1), xaxt = "n", xlab = "",
     yaxt = "n")

axis(side = 2, at = seq(-2,1, by = 1))

mtext("Emission rate error: estimated - true [kg/hr]", side = 3, line = line.val, cex = cex.val, adj = adj.val)

rect(xleft  = seq(1,n.exp, by = 2) - 0.5,
     xright = seq(1,n.exp, by = 2) + 0.5,
     ybottom = -999, ytop = 999, 
     col = alpha("black", 0.2),
     border = NA)

abline(h = 0, lty = 2, col = "red", lwd = horiz.lwd.val)



times <- data$times

for (i in 1:n.exp){
  
  interval.length <- 30
  step.size <- 10
  
  num.intervals <- (length(times)/step.size) - (interval.length/step.size) + 1
  num.intervals <- floor(num.intervals)
  
  q.hat <- results[[i]]$q.hat #/ 3.6
  q.hat.lower <- results[[i]]$q.hat.lower #/ 3.6
  q.hat.upper <- results[[i]]$q.hat.upper #/ 3.6
  
  num.to.avg <- (interval.length/step.size) - 1
  q.hat.avg <- q.hat.lower.avg <- q.hat.upper.avg <- matrix(NA, nrow = num.intervals, ncol = 5)
  q.hat.avg[1, ] <- q.hat[1, ]
  q.hat.lower.avg[1, ] <- q.hat.lower[1, ]
  q.hat.upper.avg[1, ] <- q.hat.upper[1, ]
  
  for (aa in 2:num.intervals){
    ind.to.avg <- seq(max(1, aa-num.to.avg), aa)
    q.hat.avg[aa, ] <- apply(q.hat[ind.to.avg, ], 2, mean, na.rm = T)
    q.hat.lower.avg[aa, ] <- apply(q.hat.lower[ind.to.avg, ], 2, mean, na.rm = T)
    q.hat.upper.avg[aa, ] <- apply(q.hat.upper[ind.to.avg, ], 2, mean, na.rm = T)
  }
  
  
  q.hat.total <- apply(q.hat.avg, 1, sum)
  q.hat.total.lower <- apply(q.hat.lower.avg, 1, sum)
  q.hat.total.upper <- apply(q.hat.upper.avg, 1, sum)
  
  truth.total <- apply(truth, 1, sum)
  
  for (j in 1:n.beta){
    
    these.q <- q.hat.avg[,j]
    these.true <- truth[,j] 
    these.error <- these.q - these.true
    
    arrows(x0 = i + offset[j], 
           y0 = quantile(these.error, probs = 0.025, na.rm = T),
           y1 = quantile(these.error, probs = 0.975, na.rm = T),
           col = "black", lwd = lwd.val,
           code = 3, angle = 90, length = length.val)
    
    points(i + offset[j], mean(these.error, na.rm = T), pch = 19, col = cols[j], cex = pt.cex.val)
    
  }
  
  
  these.error <- q.hat.total - truth.total
  
  arrows(x0 = i+offset[length(offset)],
         y0 = quantile(these.error, probs = 0.025, na.rm = T),
         y1 = quantile(these.error, probs = 0.975, na.rm = T),
         col = "black", lwd = lwd.val,
         code = 3, angle = 90, length = length.val)
  
  points(i + offset[length(offset)], mean(these.error, na.rm = T), pch = 15, col = "black", cex = pt.cex.val)
  
}



plot(1,1, col = "white", xlim = c(0.5+xlim.fix,n.exp+0.5-xlim.fix), ylim = c(0.5, 1), xaxt = "n", xlab = "",
     yaxt = "n")

axis(side = 2, at = seq(0.5, 1, length.out = 3))
axis(side = 2, at = seq(0.5, 1, length.out = 5), labels = NA)

mtext("Coverage", side = 3, line = line.val, cex = cex.val, adj = adj.val)

rect(xleft  = seq(1,n.exp, by = 2) - 0.5,
     xright = seq(1,n.exp, by = 2) + 0.5,
     ybottom = -999, ytop = 999, 
     col = alpha("black", 0.2),
     border = NA)

abline(h = 0.95, lty = 2, col = "red", lwd = horiz.lwd.val)

abline(h = 0.5793419, lty = 4, col = "gray40", lwd = horiz.lwd.val)

site.coverage <- vector(length = n.exp)

coverage <- matrix(NA, nrow = n.beta, ncol = n.exp)

for (i in 1:n.exp){
  
  interval.length <- 30
  step.size <- 10
  
  num.intervals <- (length(times)/step.size) - (interval.length/step.size) + 1
  num.intervals <- floor(num.intervals)
  
  q.hat <- results[[i]]$q.hat #/ 3.6
  q.hat.lower <- results[[i]]$q.hat.lower #/ 3.6
  q.hat.upper <- results[[i]]$q.hat.upper #/ 3.6
  
  num.to.avg <- (interval.length/step.size) - 1
  q.hat.avg <- q.hat.lower.avg <- q.hat.upper.avg <- matrix(NA, nrow = num.intervals, ncol = 5)
  q.hat.avg[1, ] <- q.hat[1, ]
  q.hat.lower.avg[1, ] <- q.hat.lower[1, ]
  q.hat.upper.avg[1, ] <- q.hat.upper[1, ]
  
  for (aa in 2:num.intervals){
    ind.to.avg <- seq(max(1, aa-num.to.avg), aa)
    q.hat.avg[aa, ] <- apply(q.hat[ind.to.avg, ], 2, mean, na.rm = T)
    q.hat.lower.avg[aa, ] <- apply(q.hat.lower[ind.to.avg, ], 2, mean, na.rm = T)
    q.hat.upper.avg[aa, ] <- apply(q.hat.upper[ind.to.avg, ], 2, mean, na.rm = T)
  }
  
  
  q.hat.total <- apply(q.hat.avg, 1, sum)
  q.hat.total.lower <- apply(q.hat.lower.avg, 1, sum)
  q.hat.total.upper <- apply(q.hat.upper.avg, 1, sum)
  
  for (j in 1:n.beta){
    
    betas.true <- truth[,j]
    
    these.lower <- q.hat.lower.avg[,j]
    these.upper <- q.hat.upper.avg[,j]
    
    num.within <- betas.true >= these.lower & betas.true <= these.upper
    
    coverage[j,i] <- sum(num.within, na.rm = T) / sum(!is.na(these.upper))
    
  }
  
  points(rep(i, n.beta)+offset[1:5], coverage[,i], pch = 19, col = cols, cex = pt.cex.val)
  
  
  
  truth.total <- apply(truth, 1, sum)
  num.within <- truth.total >= q.hat.total.lower & truth.total <= q.hat.total.upper
  
  site.coverage[i] <- sum(num.within, na.rm = T) / sum(!is.na(q.hat.total.upper))
  
  points(i + offset[length(offset)], site.coverage[i], pch = 15, col = "black", cex = pt.cex.val)
  
}


axis(side = 1, at = 1:n.exp, labels = seq(0,50,length.out = 5))

mtext("Concentration enhancement misalignment [%]", side = 1, line = 1.8, cex = cex.val)


# legend("top", c("West Wellhead", "West Separator", "Tank", "East Wellhead", "East Separator", "Site Total"),
#        pch = c(rep(19, 5), 15), col = c(cols, "black"), pt.cex = 1.2)

# legend("top", c("True value"),
#        lwd = 4, lty = 2, col = "red")

dev.off()









png('../figures/simulation_study_autocorrelation.png',
    res = 100, pointsize = 32, width = 1920, height = 1080*1.2)

results <- readRDS('./results/simulation_study_output_autocorrelation.RData')

laytout.mat <- matrix(c(1, 5,
                        2, 6,
                        3, 7,
                        4, 8), nrow = 4, byrow = T)

layout(laytout.mat)

par(mar = c(1.5,2.5,0.75,0))

par(oma = c(1.5,0,0.5,0.25))

par(mgp = c(4, 0.75, 0))

offset <- seq(-offset.val, offset.val, length.out = n.beta+1)
offset[1:5] <- offset[1:5]-0.05
offset[6] <- offset[6] + 0.05

xlim.fix <- 0.2





plot(1,1, col = "white", xlim = c(0.5+xlim.fix,n.exp+0.5-xlim.fix), ylim = c(0, 12), xaxt = "n", xlab = "",
     yaxt = "n")

axis(side = 2, at = seq(0, 12, by= 4))

mtext("Error variance [kg/hr]", side = 3, line = line.val, cex = cex.val, adj = adj.val)

rect(xleft  = seq(1,n.exp, by = 2) - 0.5,
     xright = seq(1,n.exp, by = 2) + 0.5,
     ybottom = -999, ytop = 999, 
     col = alpha("black", 0.2),
     border = NA)

abline(h = 2, lty = 2, col = "red", lwd = horiz.lwd.val)


for (i in 1:n.exp){
  
  arrows(x0 = i,
         y0 = quantile(results[[i]]$sigma2, probs = 0.025, na.rm = T),
         y1 = quantile(results[[i]]$sigma2, probs = 0.975, na.rm = T),
         col = "black", lwd = lwd.val,
         code = 3, angle = 90, length = length.val)
  
  points(i, mean(results[[i]]$sigma2, na.rm = T), pch = 19, col = "black", cex = pt.cex.val)
}


plot(1,1, col = "white", xlim = c(0.5+xlim.fix,n.exp+0.5-xlim.fix), ylim = c(0, 1), xaxt = "n", xlab = "",
     yaxt = "n")

axis(side = 2, at = seq(0, 1, length.out = 5))

mtext("Autocorrelation coefficient", side = 3, line = line.val, cex = cex.val, adj = adj.val)

rect(xleft  = seq(1,n.exp, by = 2) - 0.5,
     xright = seq(1,n.exp, by = 2) + 0.5,
     ybottom = -999, ytop = 999, 
     col = alpha("black", 0.2),
     border = NA)

abline(h = 2, lty = 2, col = "red", lwd = horiz.lwd.val)


for (i in 1:n.exp){
  
  segments(y0 = ac.vals[i],
           x0 = rep(i, n.beta) - 0.5,
           x1 = rep(i, n.beta) + 0.5,
           col = "red", lty =2, lwd = horiz.lwd.val)
  
  arrows(x0 = i,
         y0 = quantile(results[[i]]$r, probs = 0.025, na.rm = T),
         y1 = quantile(results[[i]]$r, probs = 0.975, na.rm = T),
         col = "black", lwd = lwd.val,
         code = 3, angle = 90, length = length.val)
  
  points(i, mean(results[[i]]$r, na.rm = T), pch = 19, col = "black", cex = pt.cex.val)
}



plot(1,1, col = "white", xlim = c(0.5+xlim.fix,n.exp+0.5-xlim.fix), ylim = c(-1.5,1.5), xaxt = "n", xlab = "",
     yaxt = 'n')

axis(side = 2, seq(-1.5, 1.5, length.out = 3))
axis(side = 2, seq(-1.5, 1.5, length.out = 5), labels = NA)

mtext("Emission rate error: estimated - true [kg/hr]", side = 3, line = line.val, cex = cex.val, adj = adj.val)

rect(xleft  = seq(1,n.exp, by = 2) - 0.5,
     xright = seq(1,n.exp, by = 2) + 0.5,
     ybottom = -999, ytop = 999, 
     col = alpha("black", 0.2),
     border = NA)

abline(h = 0, lty = 2, col = "red", lwd = horiz.lwd.val)



times <- data$times

for (i in 1:n.exp){
  
  interval.length <- 30
  step.size <- 10
  
  num.intervals <- (length(times)/step.size) - (interval.length/step.size) + 1
  num.intervals <- floor(num.intervals)
  
  q.hat <- results[[i]]$q.hat #/ 3.6
  q.hat.lower <- results[[i]]$q.hat.lower #/ 3.6
  q.hat.upper <- results[[i]]$q.hat.upper #/ 3.6
  
  num.to.avg <- (interval.length/step.size) - 1
  q.hat.avg <- q.hat.lower.avg <- q.hat.upper.avg <- matrix(NA, nrow = num.intervals, ncol = 5)
  q.hat.avg[1, ] <- q.hat[1, ]
  q.hat.lower.avg[1, ] <- q.hat.lower[1, ]
  q.hat.upper.avg[1, ] <- q.hat.upper[1, ]
  
  for (aa in 2:num.intervals){
    ind.to.avg <- seq(max(1, aa-num.to.avg), aa)
    q.hat.avg[aa, ] <- apply(q.hat[ind.to.avg, ], 2, mean, na.rm = T)
    q.hat.lower.avg[aa, ] <- apply(q.hat.lower[ind.to.avg, ], 2, mean, na.rm = T)
    q.hat.upper.avg[aa, ] <- apply(q.hat.upper[ind.to.avg, ], 2, mean, na.rm = T)
  }
  
  q.hat.total <- apply(q.hat.avg, 1, sum)
  q.hat.total.lower <- apply(q.hat.lower.avg, 1, sum)
  q.hat.total.upper <- apply(q.hat.upper.avg, 1, sum)
  
  truth.total <- apply(truth, 1, sum)
  
  for (j in 1:n.beta){
    
    these.q <- q.hat.avg[,j]
    these.true <- truth[,j] 
    these.error <- these.q - these.true
    
    arrows(x0 = i + offset[j], 
           y0 = quantile(these.error, probs = 0.025, na.rm = T),
           y1 = quantile(these.error, probs = 0.975, na.rm = T),
           col = "black", lwd = lwd.val,
           code = 3, angle = 90, length = length.val)
    
    points(i + offset[j], mean(these.error, na.rm = T), pch = 19, col = cols[j], cex = pt.cex.val)
    
  }
  
  
  these.error <- q.hat.total - truth.total
  
  arrows(x0 = i + offset[length(offset)],
         y0 = quantile(these.error, probs = 0.025, na.rm = T),
         y1 = quantile(these.error, probs = 0.975, na.rm = T),
         col = "black", lwd = lwd.val,
         code = 3, angle = 90, length = length.val)
  
  points(i + offset[length(offset)], mean(these.error, na.rm = T), pch = 15, col = "black", cex = pt.cex.val)
  
}



plot(1,1, col = "white", xlim = c(0.5+xlim.fix,n.exp+0.5-xlim.fix), ylim = c(0.47, 1), xaxt = "n", xlab = "",
     yaxt = "n")

axis(side = 2, at = seq(0.5, 1, length.out = 3))
axis(side = 2, at = seq(0.5, 1, length.out = 5), labels = NA)


mtext("Coverage", side = 3, line = line.val, cex = cex.val, adj = adj.val)

rect(xleft  = seq(1,n.exp, by = 2) - 0.5,
     xright = seq(1,n.exp, by = 2) + 0.5,
     ybottom = -999, ytop = 999, 
     col = alpha("black", 0.2),
     border = NA)

abline(h = 0.95, lty = 2, col = "red", lwd = horiz.lwd.val)

site.coverage <- vector(length = n.exp)

coverage <- matrix(NA, nrow = n.beta, ncol = n.exp)

for (i in 1:n.exp){
  
  interval.length <- 30
  step.size <- 10
  
  num.intervals <- (length(times)/step.size) - (interval.length/step.size) + 1
  num.intervals <- floor(num.intervals)
  
  q.hat <- results[[i]]$q.hat #/ 3.6
  q.hat.lower <- results[[i]]$q.hat.lower #/ 3.6
  q.hat.upper <- results[[i]]$q.hat.upper #/ 3.6
  
  num.to.avg <- (interval.length/step.size) - 1
  q.hat.avg <- q.hat.lower.avg <- q.hat.upper.avg <- matrix(NA, nrow = num.intervals, ncol = 5)
  q.hat.avg[1, ] <- q.hat[1, ]
  q.hat.lower.avg[1, ] <- q.hat.lower[1, ]
  q.hat.upper.avg[1, ] <- q.hat.upper[1, ]
  
  for (aa in 2:num.intervals){
    ind.to.avg <- seq(max(1, aa-num.to.avg), aa)
    q.hat.avg[aa, ] <- apply(q.hat[ind.to.avg, ], 2, mean, na.rm = T)
    q.hat.lower.avg[aa, ] <- apply(q.hat.lower[ind.to.avg, ], 2, mean, na.rm = T)
    q.hat.upper.avg[aa, ] <- apply(q.hat.upper[ind.to.avg, ], 2, mean, na.rm = T)
  }
  
  
  q.hat.total <- apply(q.hat.avg, 1, sum)
  q.hat.total.lower <- apply(q.hat.lower.avg, 1, sum)
  q.hat.total.upper <- apply(q.hat.upper.avg, 1, sum)
  
  truth.total <- apply(truth, 1, sum)
  
  for (j in 1:n.beta){
    
    betas.true <- truth[,j]
    
    these.lower <- q.hat.lower.avg[,j]
    these.upper <- q.hat.upper.avg[,j]
    
    num.within.source <- betas.true >= these.lower & betas.true <= these.upper
    
    coverage[j,i] <- sum(num.within.source, na.rm = T) / sum(!is.na(these.upper))
    
  }
  
  points(rep(i, n.beta)+offset[1:5], coverage[,i], pch = 19, col = cols, cex = pt.cex.val)
  
  
  
  num.within <- truth.total >= q.hat.total.lower & truth.total <= q.hat.total.upper
  
  # zero.mask <- truth.total > 0
  # sum(num.within[zero.mask], na.rm = T) / sum(!is.na(q.hat.total[zero.mask]))
  
  site.coverage[i] <- sum(num.within, na.rm = T) / sum(!is.na(q.hat.total.upper))
  
  
  points(i + offset[length(offset)], site.coverage[i], pch = 15, col = "black", cex = pt.cex.val)
}


axis(side = 1, at = 1:n.exp, labels = ac.vals)

mtext("Autocorrelation coefficient", side = 1, line = 1.8, cex = cex.val)








results <- readRDS('./results/simulation_study_output_autocorrelation_no_ac_model.RData')



plot(1,1, col = "white", xlim = c(0.5+xlim.fix,n.exp+0.5-xlim.fix), ylim = c(0, 12), xaxt = "n", xlab = "",
     yaxt = "n")

axis(side = 2, at = seq(0, 12, by= 4))

mtext("Error variance [kg/hr]", side = 3, line = line.val, cex = cex.val, adj = adj.val)

rect(xleft  = seq(1,n.exp, by = 2) - 0.5,
     xright = seq(1,n.exp, by = 2) + 0.5,
     ybottom = -999, ytop = 999, 
     col = alpha("black", 0.2),
     border = NA)

abline(h = 2, lty = 2, col = "red", lwd = horiz.lwd.val)


for (i in 1:n.exp){
  
  arrows(x0 = i,
         y0 = quantile(results[[i]]$sigma2, probs = 0.025, na.rm = T),
         y1 = quantile(results[[i]]$sigma2, probs = 0.975, na.rm = T),
         col = "black", lwd = lwd.val,
         code = 3, angle = 90, length = length.val)
  
  points(i, mean(results[[i]]$sigma2, na.rm = T), pch = 19, col = "black", cex = pt.cex.val)
}


plot(1,1, col = "white", xlim = c(0.5+xlim.fix,n.exp+0.5-xlim.fix), ylim = c(0, 1), xaxt = "n", xlab = "",
     yaxt = "n")

axis(side = 2, at = seq(0, 1, length.out = 5))

mtext("Autocorrelation coefficient", side = 3, line = line.val, cex = cex.val, adj = adj.val)

rect(xleft  = seq(1,n.exp, by = 2) - 0.5,
     xright = seq(1,n.exp, by = 2) + 0.5,
     ybottom = -999, ytop = 999, 
     col = alpha("black", 0.2),
     border = NA)

abline(h = 2, lty = 2, col = "red", lwd = horiz.lwd.val)


for (i in 1:n.exp){
  
  segments(y0 = ac.vals[i],
           x0 = rep(i, n.beta) - 0.5,
           x1 = rep(i, n.beta) + 0.5,
           col = "red", lty =2, lwd = horiz.lwd.val)
  
  arrows(x0 = i,
         y0 = quantile(results[[i]]$r, probs = 0.025, na.rm = T),
         y1 = quantile(results[[i]]$r, probs = 0.975, na.rm = T),
         col = "black", lwd = lwd.val,
         code = 3, angle = 90, length = length.val)
  
  points(i, mean(results[[i]]$r, na.rm = T), pch = 19, col = "black", cex = pt.cex.val)
}





plot(1,1, col = "white", xlim = c(0.5+xlim.fix,n.exp+0.5-xlim.fix), ylim = c(-1.5,1.5), xaxt = "n", xlab = "",yaxt = "n")

axis(side = 2, seq(-1.5, 1.5, length.out = 3))
axis(side = 2, seq(-1.5, 1.5, length.out = 5), labels = NA)

mtext("Emission rate error: estimated - true [kg/hr]", side = 3, line = line.val, cex = cex.val, adj = adj.val)

rect(xleft  = seq(1,n.exp, by = 2) - 0.5,
     xright = seq(1,n.exp, by = 2) + 0.5,
     ybottom = -999, ytop = 999, 
     col = alpha("black", 0.2),
     border = NA)

abline(h = 0, lty = 2, col = "red", lwd = horiz.lwd.val)



times <- data$times

for (i in 1:n.exp){
  
  interval.length <- 30
  step.size <- 10
  
  num.intervals <- (length(times)/step.size) - (interval.length/step.size) + 1
  num.intervals <- floor(num.intervals)
  
  q.hat <- results[[i]]$q.hat #/ 3.6
  q.hat.lower <- results[[i]]$q.hat.lower #/ 3.6
  q.hat.upper <- results[[i]]$q.hat.upper #/ 3.6
  
  num.to.avg <- (interval.length/step.size) - 1
  q.hat.avg <- q.hat.lower.avg <- q.hat.upper.avg <- matrix(NA, nrow = num.intervals, ncol = 5)
  q.hat.avg[1, ] <- q.hat[1, ]
  q.hat.lower.avg[1, ] <- q.hat.lower[1, ]
  q.hat.upper.avg[1, ] <- q.hat.upper[1, ]
  
  for (aa in 2:num.intervals){
    ind.to.avg <- seq(max(1, aa-num.to.avg), aa)
    q.hat.avg[aa, ] <- apply(q.hat[ind.to.avg, ], 2, mean, na.rm = T)
    q.hat.lower.avg[aa, ] <- apply(q.hat.lower[ind.to.avg, ], 2, mean, na.rm = T)
    q.hat.upper.avg[aa, ] <- apply(q.hat.upper[ind.to.avg, ], 2, mean, na.rm = T)
  }
  
  q.hat.total <- apply(q.hat.avg, 1, sum)
  q.hat.total.lower <- apply(q.hat.lower.avg, 1, sum)
  q.hat.total.upper <- apply(q.hat.upper.avg, 1, sum)
  
  truth.total <- apply(truth, 1, sum)
  
  for (j in 1:n.beta){
    
    these.q <- q.hat.avg[,j]
    these.true <- truth[,j] 
    these.error <- these.q - these.true
    
    arrows(x0 = i + offset[j], 
           y0 = quantile(these.error, probs = 0.025, na.rm = T),
           y1 = quantile(these.error, probs = 0.975, na.rm = T),
           col = "black", lwd = lwd.val,
           code = 3, angle = 90, length = length.val)
    
    points(i + offset[j], mean(these.error, na.rm = T), pch = 19, col = cols[j], cex = pt.cex.val)
    
  }
  
  
  these.error <- q.hat.total - truth.total
  
  arrows(x0 = i + offset[length(offset)],
         y0 = quantile(these.error, probs = 0.025, na.rm = T),
         y1 = quantile(these.error, probs = 0.975, na.rm = T),
         col = "black", lwd = lwd.val,
         code = 3, angle = 90, length = length.val)
  
  points(i + offset[length(offset)], mean(these.error, na.rm = T), pch = 15, col = "black", cex = pt.cex.val)
  
}






plot(1,1, col = "white", xlim = c(0.5+xlim.fix,n.exp+0.5-xlim.fix), ylim = c(0.47, 1), xaxt = "n", xlab = "",
     yaxt = "n")

axis(side = 2, at = seq(0.5, 1, length.out = 3))
axis(side = 2, at = seq(0.5, 1, length.out = 5), labels = NA)

mtext("Coverage", side = 3, line = line.val, cex = cex.val, adj = adj.val)

rect(xleft  = seq(1,n.exp, by = 2) - 0.5,
     xright = seq(1,n.exp, by = 2) + 0.5,
     ybottom = -999, ytop = 999, 
     col = alpha("black", 0.2),
     border = NA)

abline(h = 0.95, lty = 2, col = "red", lwd = horiz.lwd.val)

site.coverage <- vector(length = n.exp)

coverage <- matrix(NA, nrow = n.beta, ncol = n.exp)

for (i in 1:n.exp){
  
  interval.length <- 30
  step.size <- 10
  
  num.intervals <- (length(times)/step.size) - (interval.length/step.size) + 1
  num.intervals <- floor(num.intervals)
  
  q.hat <- results[[i]]$q.hat #/ 3.6
  q.hat.lower <- results[[i]]$q.hat.lower #/ 3.6
  q.hat.upper <- results[[i]]$q.hat.upper #/ 3.6
  
  num.to.avg <- (interval.length/step.size) - 1
  q.hat.avg <- q.hat.lower.avg <- q.hat.upper.avg <- matrix(NA, nrow = num.intervals, ncol = 5)
  q.hat.avg[1, ] <- q.hat[1, ]
  q.hat.lower.avg[1, ] <- q.hat.lower[1, ]
  q.hat.upper.avg[1, ] <- q.hat.upper[1, ]
  
  for (aa in 2:num.intervals){
    ind.to.avg <- seq(max(1, aa-num.to.avg), aa)
    q.hat.avg[aa, ] <- apply(q.hat[ind.to.avg, ], 2, mean, na.rm = T)
    q.hat.lower.avg[aa, ] <- apply(q.hat.lower[ind.to.avg, ], 2, mean, na.rm = T)
    q.hat.upper.avg[aa, ] <- apply(q.hat.upper[ind.to.avg, ], 2, mean, na.rm = T)
  }
  
  
  q.hat.total <- apply(q.hat.avg, 1, sum)
  q.hat.total.lower <- apply(q.hat.lower.avg, 1, sum)
  q.hat.total.upper <- apply(q.hat.upper.avg, 1, sum)
  
  truth.total <- apply(truth, 1, sum)
  
  for (j in 1:n.beta){
    
    betas.true <- truth[,j]
    
    these.lower <- q.hat.lower.avg[,j]
    these.upper <- q.hat.upper.avg[,j]
    
    num.within.source <- betas.true >= these.lower & betas.true <= these.upper
    
    coverage[j,i] <- sum(num.within.source, na.rm = T) / sum(!is.na(these.upper))
    
  }
  
  points(rep(i, n.beta)+offset[1:5], coverage[,i], pch = 19, col = cols, cex = pt.cex.val)
  
  
  
  num.within <- truth.total >= q.hat.total.lower & truth.total <= q.hat.total.upper
  
  # zero.mask <- truth.total > 0
  # sum(num.within[zero.mask], na.rm = T) / sum(!is.na(q.hat.total[zero.mask]))
  
  site.coverage[i] <- sum(num.within, na.rm = T) / sum(!is.na(q.hat.total.upper))
  
  
  points(i + offset[length(offset)], site.coverage[i], pch = 15, col = "black", cex = pt.cex.val)
}




axis(side = 1, at = 1:n.exp, labels = ac.vals)

mtext("Autocorrelation coefficient", side = 1, line = 1.8, cex = cex.val)

dev.off()



