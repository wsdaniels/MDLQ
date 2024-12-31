rm(list = ls())

library(scales)
library(viridis)
library(lubridate)

cex.val <- 0.75
pt.cex.val <- 1.15 #1.25
lwd.val <- 3
horiz.lwd.val <- 1.5
length.val <- 0.15
adj.val <- 0
line.val <- 0.1

offset.val <- 0.3

n.exp <- 5
n.beta <- 5

cols <- rev(viridis(n.beta))


data <- readRDS('../input_data/forward_model_output_ADED2024.RData')

leak.data <- readRDS('../input_data/leak_data.RData')
leak.data[is.na(leak.data)] <- 0
leak.data <- leak.data[,c(1,2,3,4,6,7,5)]

durations <- as.numeric(difftime(leak.data$end, leak.data$start, units = "min"))
to.keep <- durations >= 30
leak.data <- leak.data[to.keep, ]

to.remove <- leak.data$end > max(data$times)
leak.data <- leak.data[!to.remove, ]

leak.data.mat <- leak.data[,3:7]

n.ints <- nrow(leak.data)

sigma2.true <- 2

ac.vals <- c(0, 0.25, 0.5, 0.75, 0.95)


png('../figures/simulation_study.png',
    res = 100, pointsize = 32, width = 1920, height = 1080*1.2)

results <- readRDS('./results/simulation_study_output_autocorrelation_condensed.RData')

laytout.mat <- matrix(c(1, 5,
                        2, 6,
                        3, 7,
                        4, 8), nrow = 4, byrow = T)

layout(laytout.mat)

par(mar = c(1.5,2.5,0.75,0))

par(oma = c(1.5,0,1.5,0.25))

par(mgp = c(4, 0.75, 0))

offset <- seq(-offset.val, offset.val, length.out = n.beta)

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



plot(1,1, col = "white", xlim = c(0.5+xlim.fix,n.exp+0.5-xlim.fix), ylim = c(-0.25, 0.25), xaxt = "n", xlab = "")

mtext("Emission rate error: estimated - true [kg/hr]", side = 3, line = line.val, cex = cex.val, adj = adj.val)

rect(xleft  = seq(1,n.exp, by = 2) - 0.5,
     xright = seq(1,n.exp, by = 2) + 0.5,
     ybottom = -999, ytop = 999, 
     col = alpha("black", 0.2),
     border = NA)

abline(h = 0, lty = 2, col = "red", lwd = horiz.lwd.val)

for (i in 1:n.exp){
  
  for (j in 1:n.beta){
    
    these.q <- results[[i]]$q.hat[,j] / 3.6
    these.true <- leak.data.mat[,j] 
    these.error <- these.q - these.true
    
    arrows(x0 = i + offset[j], 
           y0 = quantile(these.error, probs = 0.025, na.rm = T),
           y1 = quantile(these.error, probs = 0.975, na.rm = T),
           col = "black", lwd = lwd.val,
           code = 3, angle = 90, length = length.val)
    
    points(i + offset[j], mean(these.error, na.rm = T), pch = 19, col = cols[j], cex = pt.cex.val)
    
  }
}


plot(1,1, col = "white", xlim = c(0.5+xlim.fix,n.exp+0.5-xlim.fix), ylim = c(0.9, 1), xaxt = "n", xlab = "",
     yaxt = "n")

axis(side = 2, at = seq(0.9, 1, by= 0.05))

mtext("Coverage", side = 3, line = line.val, cex = cex.val, adj = adj.val)

rect(xleft  = seq(1,n.exp, by = 2) - 0.5,
     xright = seq(1,n.exp, by = 2) + 0.5,
     ybottom = -999, ytop = 999, 
     col = alpha("black", 0.2),
     border = NA)

abline(h = 0.95, lty = 2, col = "red", lwd = horiz.lwd.val)

coverage <- matrix(NA, nrow = n.beta, ncol = n.exp)

for (i in 1:n.exp){
  
  for (j in 1:n.beta){
    
    betas.true <- leak.data.mat[,j]
    
    these.lower <- results[[i]]$q.hat.lower[,j]/3.6
    these.upper <- results[[i]]$q.hat.upper[,j]/3.6
    
    num.within <- betas.true >= these.lower & betas.true <= these.upper
    
    coverage[j,i] <- sum(num.within, na.rm = T) / sum(!is.na(results[[i]]$q.hat[,j]))
    
  }
  
  points(rep(i, n.beta)+offset, coverage[,i], pch = 19, col = cols, cex = pt.cex.val)
  
}



axis(side = 1, at = 1:n.exp, labels = ac.vals)

mtext("Autocorrelation coefficient", side = 1, line = 1.8, cex = cex.val)



##############################################################
##############################################################
##############################################################
##############################################################

results <- readRDS('./results/simulation_study_output_spike_misalignment_condensed.RData')
outlier.ind <- which(results[[3]]$q.hat[,1] > 500)
results[[3]]$q.hat[outlier.ind,1] <- NA

plot(1,1, col = "white", xlim = c(0.5+xlim.fix,n.exp+0.5-xlim.fix), ylim = c(0, 300), xaxt = "n", xlab = "",
     yaxt = "n")

axis(side = 2, at = seq(0, 300, length.out = 4))

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

abline(h = 0, lty = 2, col = "red", lwd = horiz.lwd.val)


for (i in 1:n.exp){
  
  
  arrows(x0 = i,
         y0 = quantile(results[[i]]$r, probs = 0.025, na.rm = T),
         y1 = quantile(results[[i]]$r, probs = 0.975, na.rm = T),
         col = "black", lwd = lwd.val,
         code = 3, angle = 90, length = length.val)
  
  points(i, mean(results[[i]]$r, na.rm = T), pch = 19, col = "black", cex = pt.cex.val)
}


plot(1,1, col = "white", xlim = c(0.5+xlim.fix,n.exp+0.5-xlim.fix), ylim = c(-3,1), xaxt = "n", xlab = "")

mtext("Emission rate error: estimated - true [kg/hr]", side = 3, line = line.val, cex = cex.val, adj = adj.val)

rect(xleft  = seq(1,n.exp, by = 2) - 0.5,
     xright = seq(1,n.exp, by = 2) + 0.5,
     ybottom = -999, ytop = 999, 
     col = alpha("black", 0.2),
     border = NA)

abline(h = 0, lty = 2, col = "red", lwd = horiz.lwd.val)

for (i in 1:n.exp){
  
  for (j in 1:n.beta){
    
    these.q <- results[[i]]$q.hat[,j] / 3.6
    these.true <- leak.data.mat[,j] 
    these.error <- these.q - these.true
    
    arrows(x0 = i + offset[j], 
           y0 = quantile(these.error, probs = 0.025, na.rm = T),
           y1 = quantile(these.error, probs = 0.975, na.rm = T),
           col = "black", lwd = lwd.val,
           code = 3, angle = 90, length = length.val)
    
    points(i + offset[j], mean(these.error, na.rm = T), pch = 19, col = cols[j], cex = pt.cex.val)
    
  }
}



plot(1,1, col = "white", xlim = c(0.5+xlim.fix,n.exp+0.5-xlim.fix), ylim = c(0.4, 1), xaxt = "n", xlab = "",
     yaxt = "n")

axis(side = 2, at = seq(0.4, 1, length.out = 4))

mtext("Coverage", side = 3, line = line.val, cex = cex.val, adj = adj.val)

rect(xleft  = seq(1,n.exp, by = 2) - 0.5,
     xright = seq(1,n.exp, by = 2) + 0.5,
     ybottom = -999, ytop = 999, 
     col = alpha("black", 0.2),
     border = NA)

abline(h = 0.95, lty = 2, col = "red", lwd = horiz.lwd.val)

coverage <- matrix(NA, nrow = n.beta, ncol = n.exp)

for (i in 1:n.exp){
  
  for (j in 1:n.beta){
    
    betas.true <- leak.data.mat[,j]
    
    these.lower <- results[[i]]$q.hat.lower[,j]/3.6
    these.upper <- results[[i]]$q.hat.upper[,j]/3.6
    
    num.within <- betas.true >= these.lower & betas.true <= these.upper
    
    coverage[j,i] <- sum(num.within, na.rm = T) / sum(!is.na(results[[i]]$q.hat[,j]))
    
  }
  
  points(rep(i, n.beta)+offset, coverage[,i], pch = 19, col = cols, cex = pt.cex.val)
  
}

axis(side = 1, at = 1:n.exp, labels = seq(0,50,length.out = 5))

mtext("Concentration enhancement misalignment [%]", side = 1, line = 1.8, cex = cex.val)


dev.off()






##############################################################
##############################################################
##############################################################
##############################################################
# INDEPENDENT ERRORS


png('../figures/simulation_study_independent_errors.png',
    res = 100, pointsize = 32, width = 1920/2, height = 1080*1.2)


results <- readRDS('./results/simulation_study_output_autocorrelation_independent_errors_condensed.RData')

laytout.mat <- matrix(c(1, 
                        2,
                        3,
                        4), nrow = 4, byrow = T)

layout(laytout.mat)

par(mar = c(1.5,2.5,0.75,0))

par(oma = c(1.5,0,1.5,0.25))

par(mgp = c(4, 0.75, 0))

offset <- seq(-offset.val, offset.val, length.out = n.beta)

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



plot(1,1, col = "white", xlim = c(0.5+xlim.fix,n.exp+0.5-xlim.fix), ylim = c(-0.4, 0.4), xaxt = "n", xlab = "")

mtext("Emission rate error: estimated - true [kg/hr]", side = 3, line = line.val, cex = cex.val, adj = adj.val)

rect(xleft  = seq(1,n.exp, by = 2) - 0.5,
     xright = seq(1,n.exp, by = 2) + 0.5,
     ybottom = -999, ytop = 999, 
     col = alpha("black", 0.2),
     border = NA)

abline(h = 0, lty = 2, col = "red", lwd = horiz.lwd.val)

for (i in 1:n.exp){
  
  for (j in 1:n.beta){
    
    these.q <- results[[i]]$q.hat[,j] / 3.6
    these.true <- leak.data.mat[,j] 
    these.error <- these.q - these.true
    
    arrows(x0 = i + offset[j], 
           y0 = quantile(these.error, probs = 0.025, na.rm = T),
           y1 = quantile(these.error, probs = 0.975, na.rm = T),
           col = "black", lwd = lwd.val,
           code = 3, angle = 90, length = length.val)
    
    points(i + offset[j], mean(these.error, na.rm = T), pch = 19, col = cols[j], cex = pt.cex.val)
    
  }
}


plot(1,1, col = "white", xlim = c(0.5+xlim.fix,n.exp+0.5-xlim.fix), ylim = c(0.4, 1), xaxt = "n", xlab = "",
     yaxt = "n")

axis(side = 2, at = seq(0.4, 1, by= 0.2))

mtext("Coverage", side = 3, line = line.val, cex = cex.val, adj = adj.val)

rect(xleft  = seq(1,n.exp, by = 2) - 0.5,
     xright = seq(1,n.exp, by = 2) + 0.5,
     ybottom = -999, ytop = 999, 
     col = alpha("black", 0.2),
     border = NA)

abline(h = 0.95, lty = 2, col = "red", lwd = horiz.lwd.val)

coverage <- matrix(NA, nrow = n.beta, ncol = n.exp)

for (i in 1:n.exp){
  
  for (j in 1:n.beta){
    
    betas.true <- leak.data.mat[,j]
    
    these.lower <- results[[i]]$q.hat.lower[,j]/3.6
    these.upper <- results[[i]]$q.hat.upper[,j]/3.6
    
    num.within <- betas.true >= these.lower & betas.true <= these.upper
    
    coverage[j,i] <- sum(num.within, na.rm = T) / sum(!is.na(results[[i]]$q.hat[,j]))
    
  }
  
  points(rep(i, n.beta)+offset, coverage[,i], pch = 19, col = cols, cex = pt.cex.val)
  
}



axis(side = 1, at = 1:n.exp, labels = ac.vals)

mtext("Autocorrelation coefficient", side = 1, line = 1.8, cex = cex.val)

dev.off()





##############################################################
##############################################################
##############################################################
##############################################################
# Laplace errors


png('../figures/simulation_study_laplace.png',
    res = 100, pointsize = 32, width = 1920/2, height = 1080*1.2)


results <- readRDS('./results/simulation_study_output_autocorrelation_laplace_errors_condensed.RData')

laytout.mat <- matrix(c(1, 
                        2,
                        3,
                        4), nrow = 4, byrow = T)

layout(laytout.mat)

par(mar = c(1.5,2.5,0.75,0))

par(oma = c(1.5,0,1.5,0.25))

par(mgp = c(4, 0.75, 0))

offset <- seq(-offset.val, offset.val, length.out = n.beta)

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


for (i in 1){
  
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


for (i in 1){
  
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


plot(1,1, col = "white", xlim = c(0.5+xlim.fix,n.exp+0.5-xlim.fix), ylim = c(-0.25, 3), xaxt = "n", xlab = "")

mtext("Emission rate error: estimated - true [kg/hr]", side = 3, line = line.val, cex = cex.val, adj = adj.val)

rect(xleft  = seq(1,n.exp, by = 2) - 0.5,
     xright = seq(1,n.exp, by = 2) + 0.5,
     ybottom = -999, ytop = 999, 
     col = alpha("black", 0.2),
     border = NA)

abline(h = 0, lty = 2, col = "red", lwd = horiz.lwd.val)

for (i in 1){
  
  for (j in 1:n.beta){
    
    these.q <- results[[i]]$q.hat[,j] / 3.6
    these.true <- leak.data.mat[,j] 
    these.error <- these.q - these.true
    
    arrows(x0 = i + offset[j], 
           y0 = quantile(these.error, probs = 0.025, na.rm = T),
           y1 = quantile(these.error, probs = 0.975, na.rm = T),
           col = "black", lwd = lwd.val,
           code = 3, angle = 90, length = length.val)
    
    points(i + offset[j], mean(these.error, na.rm = T), pch = 19, col = cols[j], cex = pt.cex.val)
    
  }
}


plot(1,1, col = "white", xlim = c(0.5+xlim.fix,n.exp+0.5-xlim.fix), ylim = c(0.3, 1), xaxt = "n", xlab = "",
     yaxt = "n")

axis(side = 2, at = seq(0.3, 1, by= 0.3))

mtext("Coverage", side = 3, line = line.val, cex = cex.val, adj = adj.val)

rect(xleft  = seq(1,n.exp, by = 2) - 0.5,
     xright = seq(1,n.exp, by = 2) + 0.5,
     ybottom = -999, ytop = 999, 
     col = alpha("black", 0.2),
     border = NA)

abline(h = 0.95, lty = 2, col = "red", lwd = horiz.lwd.val)

coverage <- matrix(NA, nrow = n.beta, ncol = n.exp)

for (i in 1){
  
  for (j in 1:n.beta){
    
    betas.true <- leak.data.mat[,j]
    
    these.lower <- results[[i]]$q.hat.lower[,j]/3.6
    these.upper <- results[[i]]$q.hat.upper[,j]/3.6
    
    num.within <- betas.true >= these.lower & betas.true <= these.upper
    
    coverage[j,i] <- sum(num.within, na.rm = T) / sum(!is.na(results[[i]]$q.hat[,j]))
    
  }
  
  points(rep(i, n.beta)+offset, coverage[,i], pch = 19, col = cols, cex = pt.cex.val)
  
}

axis(side = 1, at = 1:n.exp, labels = ac.vals)

mtext("Autocorrelation coefficient", side = 1, line = 1.8, cex = cex.val)

dev.off()

