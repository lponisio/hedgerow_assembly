# rm(list=ls())
setwd('~/Dropbox/hedgerow_assembly/analysis/speciesLevel')
source('src/initialize.R')
load('../../data/networks/all_networks_years.Rdata')

## **********************************************************
## species importance
## **********************************************************
## linear models
load(file=file.path(save.path, 'specs.Rdata'))

## SiteStatus or ypr
xvar <- "ypr"

## anything outputted by specieslevel
ys <- c("proportional.generality", "d", "degree", "betweenness",
        "closeness.log")

formulas <-lapply(ys, function(x) {
  as.formula(paste(x, "~",
                   paste(xvar,
                         "(1|Site)",
                          "(1|GenusSpecies)",
                         sep="+")))
})

mod.pols <- lapply(formulas, function(x){
  lmer(x,
       data=specs[specs$speciesType == "pollinator",])
})

mod.plants <- lapply(formulas, function(x){
  lmer(x,
       data=specs[specs$speciesType == "plant",])
})

names(mod.pols) <- names(mod.plants) <- ys


lapply(mod.plants, summary)
lapply(mod.pols, summary)

print(summary(mod.pols$closeness.log))
print("**************************************************")
print(summary(mod.plants$closeness.log))

save(mod.pols, mod.plants, ys, file=file.path(save.path,
            sprintf('mods/specs_%s.Rdata', xvar)))
