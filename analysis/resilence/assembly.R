rm(list=ls())
setwd('~/Dropbox/hedgerow_assembly/analysis/resilence')
source('src/initialize.R')

## either "abund" or "degree"
extinction.method <- "degree"

## create pp matrix for each site, year
nets <- break.net(spec)

## **********************************************************
## robustness
## **********************************************************
## simulation plant extinction

res <- simExtinction(nets, extinction.method, spec)

save(res, file=file.path(save.path,
            sprintf('resilience_%s.Rdata', extinction.method)))

## no change in robustness by site status
mod.status <- lmer(Robustness ~ SiteStatus
             + (1|Site) + (1|Year),
             data=res)
summary(mod.status)
save(mod.status, file=file.path(save.path,
            sprintf('mods/resilience_status_%s.Rdata', extinction.method)))


## no effect of ypr on robustness
mod.ypr <- lmer(Robustness ~ ypr
             + (1|Site) + (1|Year),
             data=res[!is.na(res$ypr),])
summary(mod.ypr)
save(mod.ypr, file=file.path(save.path,
            sprintf('mods/resilience_ypr_%s.Rdata', extinction.method)))

## **********************************************************
## species importance
## **********************************************************
specs <- calcSpec(nets, spec, spec.metric = "d", 0.3)
save(specs, file=file.path(save.path, 'specs.Rdata'))

## linear models
load(file=file.path(save.path, 'specs.Rdata'))
## SiteStatus or ypr
xvar <- "ypr"

## anything outputted by specieslevel
ys <- c("proportional.generality", "d", "degree", "betweenness",
        "closeness")

## formulas <-lapply(ys, function(x) {
##   as.formula(paste(x, "~",
##                    paste(paste(xvar, "specialization", sep="*"), 
##                          "(1|Site)",
##                           "(1|GenusSpecies)",
##                          sep="+")))
## })



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

save(mod.pols, mod.plants, ys, file=file.path(save.path,
            sprintf('mods/specs_%s.Rdata', xvar)))
