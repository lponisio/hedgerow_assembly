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
plants <- getDis(sites, 1, nets, specs.agg, traits)
pols <- getDis(sites, 2, nets, specs.agg, traits)

mod.pols <- lmer(Dist~k + (1|Site) + (1|GenusSpecies),
                 data=pols)
mod.pols.cat <- lmer(Dist~core + (1|Site) + (1|GenusSpecies),
                     data=pols)
mod.pols.d <- lmer(Dist~d + (1|Site) + (1|GenusSpecies),
                   data=pols)
mod.pols.occ <- lmer(Dist~occ.date + (1|Site) + (1|GenusSpecies),
                     data=pols)
summary(mod.pols)
summary(mod.pols.cat)
summary(mod.pols.d)
summary(mod.pols.occ)

mod.plants <- lmer(Dist~k + (1|Site) + (1|GenusSpecies),
                   data=plants)
mod.plants.cat <- lmer(Dist~core + (1|Site) + (1|GenusSpecies),
                       data=plants)
mod.plants.d <- lmer(Dist~d + (1|Site) + (1|GenusSpecies),
                     data=plants)
summary(mod.plants)
summary(mod.plants.cat)
summary(mod.plants.d)
