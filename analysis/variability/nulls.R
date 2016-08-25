rm(list=ls())
library(vegan)
library(bipartite)
library(fields)
setwd('~/Dropbox/hedgerow_assembly/analysis/variability')
source('src/initialize.R')
source('src/commPrep.R')
save.dir.comm <- "saved/communities"
save.dir.nulls <- "saved/nulls"
nnull <- 9999

src.dir2 <- '~/Dropbox/hedgerow_network/analysis/beta-div'
source(file.path(src.dir2, 'src/misc.R'))
source(file.path(src.dir2, 'src/vaznull2.R'))
sites <- unique(spec$Site)

## ************************************************************
## year by species matrices pollinators!
## ************************************************************

comms <- lapply(sites, calcYearBeta,
                species.type="GenusSpecies",
                spec=spec)
comm.pols <- makePretty(comms, sites, spec)

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

## ************************************************************
## year by species matrices interactions!
## ************************************************************
spec$Int <- paste(spec$GenusSpecies,
                  spec$PlantGenusSpecies)

comms.ints <- lapply(sites, calcYearBeta,
                species.type="Int",
                spec=spec)
comm.int <- makePretty(comms.ints, sites, spec)

save(comm.int, file=file.path(save.dir.comm, 'int-abund.Rdata'))

## ************************************************************
## occurrence nulls
## ************************************************************
load(file=file.path(save.dir.comm, 'int-abund.Rdata'))
nulls.int <- lapply(comm.int$comm, function(x){
replicate(nnull, commsimulator(x, method= 'quasiswap'),
          simplify= FALSE)
})
save(nulls.int, file=file.path(save.dir.nulls, 'int-occ.Rdata'))

## ************************************************************
## alpha div nulls
## ************************************************************
load(file=file.path(save.dir.comm, 'int-abund.Rdata'))
nulls.int <- lapply(comm.int$comm, vaznull.2, N=nnull)
save(nulls.int, file=file.path(save.dir.nulls, 'int-alpha.Rdata'))

## ************************************************************
## individuals nulls
## ************************************************************
load(file=file.path(save.dir.comm, 'int-abund.Rdata'))
nulls.int <- lapply(comm.int$comm, swap.web, N=nnull)
save(nulls.int, file=file.path(save.dir.nulls, 'int-indiv.Rdata'))
