## setwd('~/Dropbox/hedgerow_assembly/analysis/variability')
setwd('variability')
source('src/initialize.R')
source('src/commPrep.R')
save.dir.comm <- "saved/communities"
save.dir.nulls <- "saved/nulls"
nnull <- 99

src.dir2 <- 'analysis/beta-div'
source(file.path(src.dir2, 'src/misc.R'))
source(file.path(src.dir2, 'src/vaznull2.R'))
sites <- unique(spec$Site)
spec$Int <- paste(spec$GenusSpecies,
                  spec$PlantGenusSpecies)

args <- commandArgs(trailingOnly=TRUE)
type <- args[1]

## pols, int or plants
## type <- "plants"
if(type=="pols"){
  species.type="GenusSpecies"
}
if(type=="int"){
  species.type="Int"
}
if(type=="plants"){
  species.type="PlantGenusSpecies"
}

## ************************************************************
## year by species matrices pollinators!
## ************************************************************
comms <- lapply(sites, calcYearBeta,
                species.type=species.type,
                spec=spec)
comm <- makePretty(comms, sites, spec)

save(comm, file=file.path(save.dir.comm,
             sprintf('%s-abund.Rdata', type)))

## ************************************************************
## occurrence nulls
## ************************************************************
load(file=file.path(save.dir.comm,
       sprintf('%s-abund.Rdata', type)))
nulls <- lapply(comm$comm, function(x){
  replicate(nnull, commsimulator(x, method= 'quasiswap'),
            simplify= FALSE)
})
save(nulls, file=file.path(save.dir.nulls,
              sprintf('%s-occ.Rdata', type)))

## ************************************************************
## alpha div nulls
## ************************************************************
load(file=file.path(save.dir.comm,
       sprintf('%s-abund.Rdata', type)))
nulls <- lapply(comm$comm, vaznull.2, N=nnull)
save(nulls, file=file.path(save.dir.nulls,
              sprintf('%s-alpha.Rdata', type)))

## ************************************************************
## individuals nulls
## ************************************************************
load(file=file.path(save.dir.comm,
       sprintf('%s-abund.Rdata', type)))
nulls <- lapply(comm$comm, swap.web, N=nnull)
save(nulls, file=file.path(save.dir.nulls,
              sprintf('%s-indiv.Rdata', type)))
