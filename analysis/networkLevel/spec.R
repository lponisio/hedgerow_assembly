rm(list=ls())
setwd('~/Dropbox/hedgerow_assembly/analysis/networkLevel')
source('src/initialize.R')

## *******************************************************************
## comparison of pollinator generalization
## *******************************************************************

load('../../data/degree/specializations.Rdata')
sites <- unique(specializations$Site)
metrics <- colnames(specializations)
metrics <- metrics[!metrics %in% c("Site",
                                   "assem",
                                   "GenusSpecies",
                                   "speciesType")]

diff.spec <- lapply(sites, function(x){
  getDiff(specializations[specializations$Site == x,],
          metrics)})

diff.spec <- do.call(rbind, diff.spec)

res.lms.pol <- lapply(metrics, specMods,
                diff.dats = diff.spec[diff.spec$speciesType ==
                                      "pollinator",],
                type= "pol")

res.lms.plant <- lapply(metrics, specMods,
                diff.dats = diff.spec[diff.spec$speciesType ==
                                      "plant",],
                type="plant")
names(res.lms.pol) <- names(res.lms.plant) <- metrics

res.lms.pol$normalised.degree
res.lms.pol$proportional.generality

res.lms.plant$normalised.degree
res.lms.plant$proportional.generality

save(res.lms.pol, res.lms.plant, file="saved/mods/spec.Rdata")
save(diff.spec, metrics, file="saved/spec.Rdata")


## *******************************************************************
## change in the number of visits a plant receives by plant
## generalization
## *******************************************************************
load('../../data/degree/PlantvisitDiffs.Rdata')


## generalists
plot(y=diff.gens$diffVisits, x=diff.gens$d)

plot(y=diff.gens$diffVisits, x=diff.gens$proportional.generality)

summary(lmer(diffVisits ~ d +
       (1|Site) + (1|PlantGenusSpecies),
     data=diff.gens))

summary(lmer(diffVisits ~ proportional.generality +
       (1|Site) + (1|PlantGenusSpecies),
     data=diff.gens))



## specialists
plot(y=diff.spec$diffVisits, x=diff.spec$d)

plot(y=diff.spec$diffVisits, x=diff.spec$proportional.generality)

summary(lmer(diffVisits ~ d +
       (1|Site) + (1|PlantGenusSpecies),
     data=diff.spec))

summary(lmer(diffVisits ~ proportional.generality +
       (1|Site) + (1|PlantGenusSpecies),
     data=diff.spec))
