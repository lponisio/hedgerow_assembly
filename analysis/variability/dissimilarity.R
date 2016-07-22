rm(list=ls())
setwd('~/Dropbox/hedgerow_assembly/analysis/variability')
source('src/initialize.R')
library(bipartite)
source('src/calcDis.R')
traits <- read.csv("../../data/traits.csv")
## ************************************************************
## the variability fo interactions of plants and pols
## ************************************************************
## for each species at each site, 1) create a interaction partner by
## year matrix, 2) calculate the dissimilarity of interaction partners
## between years 3) extract dispersion values for each species, site

load('../../data/networks/expanded_networks.Rdata')

sites <- sapply(strsplit(names(nets), "[.]"), function(x) x[1])

specs.agg <- aggregate(k  ~ GenusSpecies, data=specs, mean)
plants <- getDis(sites, 1, nets, specs.agg, traits, spec)
pols <- getDis(sites, 2, nets, specs.agg, traits, spec)

## dipersion of interaction parteners by core/periphery (continuus
## metric) (categoric metric didn'et make sense bec so few had values
## > 1)
mod.pols <- lmer(Dist~k + (1|Site) + (1|GenusSpecies),
                 data=pols)

## only two species with k >1
hist(unique(data.frame(pols$GenusSpecies, pols$k))$pols.k)

mod.plants <- lmer(Dist~k + (1|Site) + (1|GenusSpecies),
                   data=plants)

## only one species with k >1
hist(unique(data.frame(plants$GenusSpecies, plants$k))$plants.k)

## by specialization
mod.pols.d <- lmer(Dist~d + (1|Site) + (1|GenusSpecies),
                   data=pols)

mod.plants.d <- lmer(Dist~d + (1|Site) + (1|GenusSpecies),
                     data=plants)

## pollinator date occurrence
mod.pols.occ <- lmer(Dist~occ.date + (1|Site) + (1|GenusSpecies),
                     data=pols)

summary(mod.pols)
summary(mod.pols.d)
summary(mod.pols.occ)

summary(mod.plants)
summary(mod.plants.d)


save(pols, plants, mod.pols,
     mod.plants, mod.plants.d,
     mod.pols.d, file="saved/dissMods.Rdata")
