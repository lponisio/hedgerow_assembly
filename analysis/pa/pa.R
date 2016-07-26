rm(list=ls())
setwd('~/Dropbox/hedgerow_assembly/analysis/pa')

load('../../data/networks/allSpecimens.Rdata')
f.path <- "../changePoint/cptPeel/baci"
load(file=file.path(f.path, "graphs.Rdata"))

## total interactions
plant.tot.int <- lapply(nets, rowSums)
pol.tot.int <- lapply(nets, colSums)

## degree
plant.degree <- lapply(nets, function(x)
                       apply(x, 1, function(x) sum(x >= 1)))
                      
pol.degree <- lapply(nets, function(x)
                       apply(x, 2, function(x) sum(x >= 1)))
## ypr 1 
first.yrs <- spec[spec$ypr == 1,]
first.yrs <- unique(data.frame(Site=first.yrs$Site,
                  Year=first.yrs$Year))
first.yrs <- first.yrs[!apply(first.yrs, 1, function(x) any(is.na(x))),]

ypr1 <- names(nets) %in%
                    paste(first.yrs$Site, first.yrs$Year, sep=".")


sites <- sapply(strsplit(names(nets), "[.]"), function(x) x[1])
years <- sapply(strsplit(names(nets), "[.]"), function(x) x[2])

N <- 2
lapply(unique(sites), function(x){
  this.net <- nets[sites == x]
  pres.sp <- lapply(this.net, function(y){
    y[y > 1] <- 1
    y2 <- y
    y2[rowSums(y) != 0,] <- 1
    y2[,colSums(y) != 0] <- 1
    return(y2)
  })
  browser()
})
