rm(list=ls())
library(vegan)
library(bipartite)
library(fields)
setwd('~/Dropbox/hedgerow_assembly/analysis/variability')
source('src/initialize.R')
source('src/commPrep.R')
save.dir.comm <- "saved/communities"
save.dir.nulls <- "saved/nulls"
nnull <- 999

src.dir2 <- '~/Dropbox/hedgerow_network/analysis/beta-div'
source(file.path(src.dir2, 'src/misc.R'))
source(file.path(src.dir2, 'src/vaznull2.R'))

## ************************************************************
## year by species matrices
## ************************************************************
sites <- unique(spec$Site)
comms <- lapply(sites, calcYearBeta,
                species.type="GenusSpecies",
                spec=spec)
names(comms) <- sites
nyears <- sapply(comms, nrow)
comms <- comms[nyears >= 3]
years <- lapply(comms, rownames)

comm.pols <- list(comm=comms,
                  years=years,
                  sites= rep(names(comms),
                    sapply(comms, nrow)))
statuses <- spec$SiteStatus[match(comm.pols$sites,
                                          spec$Site)]
comm.pols$status <- split(statuses,
                          comm.pols$sites)[names(comm.pols$comm)]

save(comm.pols, file=file.path(save.dir.comm, 'pols-abund.Rdata'))

## ************************************************************
## occurrence nulls
## ************************************************************
load(file=file.path(save.dir.comm, 'pols-abund.Rdata'))
nulls.pols <- lapply(comm.pols$comm, function(x){
replicate(nnull, commsimulator(x, method= 'quasiswap'),
          simplify= FALSE)
})
save(nulls.pols, file=file.path(save.dir.nulls, 'pols-occ.Rdata'))

## ************************************************************
## alpha div nulls
## ************************************************************
load(file=file.path(save.dir.comm, 'pols-abund.Rdata'))
nulls.pols <- lapply(comm.pols$comm, vaznull.2, N=nnull)
save(nulls.pols, file=file.path(save.dir.nulls, 'pols-alpha.Rdata'))

## ************************************************************
## individuals nulls
## ************************************************************
load(file=file.path(save.dir.comm, 'pols-abund.Rdata'))
nulls.pols <- lapply(comm.pols$comm, swap.web, N=nnull)
save(nulls.pols, file=file.path(save.dir.nulls, 'pols-indiv.Rdata'))
