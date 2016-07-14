rm(list=ls())
library(reshape)
setwd("~/Dropbox/hedgerow_assembly/analysis/changePoint/cptPeel")

results <- read.table('results_4.txt', sep=' ') 
logs <- read.table('LogLs_4.txt', sep=' ')

## split the year and creates an extra table with year 1 and period
results$V1 <- as.character(results$V1)
results$V2 <- as.character(results$V2)
years1  <-  as.numeric(sapply(strsplit(results$V1, "_"), function(x) x[2]))
years2  <-  as.numeric(sapply(strsplit(results$V2, "_"), function(x) x[2]))
anos <- min(years1):max(years2)
anos <- anos[-which(anos == 2010)]
years <- as.data.frame(matrix(NA, 7,2))
for(n in 1:(length(anos)-1)){
  years[n,2] <- paste(anos[n], anos[n+1], sep="-")  
}
years[,1] <- anos[-length(anos)]
colnames(years) <- c('first', 'period')
## creates othe variables for results 
sites  <- sapply(strsplit(results$V1, "_"), function(x) x[1])
max <- apply(results[,c(3:5)], 1, max)
numberOfOnes <- apply(results[,c(3:5)], 1,
                      function(x) length(which(x==1)))
max.index <- rep(NA, dim(results)[1])
results <- cbind(results,  sites, years1, max, numberOfOnes, max.index)

## checking if there is more than 1 value, and if yes, returns the
## index with the maxlog
for(row in 1:dim(results)[1]){
  if(results$numberOfOnes[row] > 1){## se eh true tem mais de um ai tem
                                  ## que ver o log
    results$max.index[row] <- which.max(logs[row,c(3:5)])
  }
  else {## return the position of the max
    results$max.index[row] <- which(results[row,c(3:5)] ==
                                        results$max[row])[1]
  }
}

## separating only the significative results
## for(i in 1:dim(results)[1]){
##  qual <- which(years$first%in%results$years1[i]==TRUE)
##  results$period[i] <- years[qual,1]
##  }

sigs <- results[results$max>0.949,]
cp <- rep(NA, dim(sigs)[1])
value <- rep(NA, dim(sigs)[1])
sigs <- cbind(sigs, cp, value)

## cp is always the year + next year
## retorna o index que eh a soma de onde ta mais o max.index
for(i in 1:dim(sigs)[1]){
  ## if(i==28) browser()
  y1 <- which(years$first %in% sigs$years1[i] == TRUE)
  print(paste(i,  years$period[y1+(sigs$max.index[i]-1)]))
  sigs$cp[i] <- years$period[y1+(sigs$max.index[i]-1)]
  sigs$value[i] <- sigs[i,c(3:5)][sigs$max.index[i]]
  
}

sigs$value <- unlist(sigs$value)
changing.points <- sigs[,c(6,11,12)]

write.csv(changing.points, "changing_points.csv", row.names=FALSE)
