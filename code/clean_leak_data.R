rm(list = ls())

library(lubridate)

data <- read.csv("../input_data/leak_data.csv")

data$tc_ExpStartDatetime <- as_datetime(data$tc_ExpStartDatetime)
data$tc_ExpEndDatetime   <- as_datetime(data$tc_ExpEndDatetime)

data.to.keep <- data$tc_EquipmentGroupID %in% c("4W", "4S", "5S", "4T", "5W")
data <- data[data.to.keep, ]

nr <- nrow(data)

this.order <- order(data$tc_ExpStartDatetime)
data <- data[this.order, ]


combine.start <- combine.end <- matrix(NA, nrow = nr, ncol = nr)

for (i in 1:nr){
  
  print(paste0(i, "/", nr))
  
  event.start <- data$tc_ExpStartDatetime[i]
  event.end   <- data$tc_ExpEndDatetime[i]
  
  for (j in 1:nr){
    
    this.start <- data$tc_ExpStartDatetime[j]
    this.end   <- data$tc_ExpEndDatetime[j]
    
    combine.start[i,j] <- abs(as.numeric(difftime(event.start, this.start, units = "min"))) <= 5
    combine.end[i,j]   <- abs(as.numeric(difftime(event.end,   this.end,   units = "min"))) <= 5
    
  }
}



start.to.combine <- vector(mode = "list")
ind <- 1
b.start <- 1

while(T){
  b.end <- max(which(combine.start[b.start,]))
  
  start.to.combine[[ind]] <- b.start:b.end
  ind <- ind + 1
  
  b.start <- b.end + 1
  if (b.start > nr) {break}
}


end.to.combine <- vector(mode = "list")
ind <- 1
b.start <- 1

while(T){
  b.end <- max(which(combine.end[b.start,]))
  
  end.to.combine[[ind]] <- b.start:b.end
  ind <- ind + 1
  
  b.start <- b.end + 1
  if (b.start > nr) {break}
}


final.keep <- vector(length = length(start.to.combine))
matching.end.ind <- vector(length = length(start.to.combine))

for (i in 1:length(start.to.combine)){
  this.vec <- start.to.combine[[i]]
  
  for (j in 1:length(end.to.combine)){
    cond1 <- all(length(this.vec) == length(end.to.combine[[j]]))
    
    if (cond1){
      cond2 <- all(this.vec == end.to.combine[[j]])
      
      if (cond2){
        final.keep[i] <- T
        matching.end.ind[i] <- j
        break
      }
    }
  }
}

start.to.combine.final <- start.to.combine[final.keep]

end.final.keep <- sort(unique(matching.end.ind))
end.final.keep <- end.final.keep[-(end.final.keep == 0)]

end.to.combine.final <- end.to.combine[end.final.keep]


combine.final <- start.to.combine.final


out <- matrix(NA, nrow = length(combine.final), ncol = 7)
out <- as.data.frame(out)
colnames(out) <- c("start", "end", paste0("EG", unique(data$tc_EquipmentGroupID)))

start.spread <- end.spread <- vector(length = length(combine.final))

for (i in 1:length(combine.final)){
  
  this.sub <- data[combine.final[[i]], ]
  
  out$start[i] <- min(this.sub$tc_ExpStartDatetime)
  
  out$end[i] <- max(this.sub$tc_ExpEndDatetime)
  
  start.spread[i] <- as.numeric(difftime(min(this.sub$tc_ExpStartDatetime), max(this.sub$tc_ExpStartDatetime), units = "min"))
  end.spread[i]   <- as.numeric(difftime(min(this.sub$tc_ExpEndDatetime),   max(this.sub$tc_ExpEndDatetime),   units = "min"))
  
  for (j in 1:nrow(this.sub)){
    
    this.eg <- this.sub$tc_EquipmentGroupID[j]
    
    if (this.eg == "4S"){
      out$EG4S[i] <- this.sub$tc_C1MassFlow[j]/1000
    } else if (this.eg == "4W"){
      out$EG4W[i] <- this.sub$tc_C1MassFlow[j]/1000
    } else if (this.eg == "4T"){
      out$EG4T[i] <- this.sub$tc_C1MassFlow[j]/1000
    } else if (this.eg == "5W"){
      out$EG5W[i] <- this.sub$tc_C1MassFlow[j]/1000
    } else if (this.eg == "5S"){
      out$EG5S[i] <- this.sub$tc_C1MassFlow[j]/1000
    }
    
  }
}


out$start <- as_datetime(out$start)
out$end <- as_datetime(out$end)

this.order <- order(out$start)

out <- out[this.order, ]

saveRDS(out, '../input_data/leak_data.RData')
