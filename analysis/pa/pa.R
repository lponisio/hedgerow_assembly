rm(list=ls())
setwd('~/Dropbox/hedgerow_assembly/analysis/pa')
library(LaplacesDemon)

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

prob.null.binary <- function(probs, fill, nrow, ncol) {
  ones <- sample(1:length(probs), fill, replace=FALSE,
                 prob=probs)
  ## create resultant matrix
  interact <- matrix(0, nrow=nrow, ncol=ncol)
  interact[ones] <- 1
  return(interact)
}

lapply(unique(sites), function(x){
  this.net <- nets[sites == x]
  pres.sp <- lapply(this.net, function(y){
    ## y[y > 1] <- 1
    degree.pol <- colSums(y) + 1
    degree.plant <- rowSums(y) + 1
    y2 <- y
    y2[rowSums(y) != 0,] <- 1
    y2[,colSums(y) != 0] <- 1
    return(list(y=y,
                fill.y = sum(y),
                pres.mat = y2,
                degree.pol = degree.pol,
                degree.plant = degree.plant))
  })
  these.years <- sapply(strsplit(names(this.net), "[.]"), function(x) x[2])
  lik.years.pa <- numeric(length(these.years)-1)
  lik.years.rand <- numeric(length(these.years)-1)
  for(i in 2:length(these.years)){
      prob.mat <- pres.sp[[i]]$pres.mat*pres.sp[[i-1]]$degree.plant
      lik.years.pa[i-1] <- dmultinom(pres.sp[[i]]$y,
                                     prob=prob.mat)
      lik.years.rand[i-1] <- dmultinom(pres.sp[[i]]$y,
                                     prob=pres.sp[[i]]$pres.mat)
    }
  names(lik.years.pa) <- names(lik.years.rand) <- these.years[-1]
  lik.ratio <- 2*(log(lik.years.pa) - log(lik.years.rand))
})
