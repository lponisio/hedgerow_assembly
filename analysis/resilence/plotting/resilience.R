rm(list=ls())
setwd('~/Dropbox/hedgerow_assembly/analysis/resilence')
source('plotting/src/predictIntervals.R')
source('plotting/src/CIplotting.R')
source('plotting/src/plotPanels.R')
source('src/initialize.R')

## ************************************************************
## robustness
## ************************************************************
## either "abund" or "degree"
extinction.method <- "abund"

load(file=file.path(save.path,
            sprintf('resilience_%s.Rdata', extinction.method)))

dd <- expand.grid(ypr=seq(from= min(res$ypr, na.rm=TRUE),
                          to= max(res$ypr, na.rm=TRUE),
                          length=10))


load(file=file.path(save.path,
            sprintf('mods/resilience_ypr_%s.Rdata',
                    extinction.method)))

dd.ypr <- cbind(dd, Robustness=0)

ypr.pi <- predict.int(mod= mod.ypr,
                        dd=dd.ypr,
                        y="Robustness",
                        family="gaussian")

plot.predict.ypr(new.dd=ypr.pi,
                 ylabel="Robustness",
                 dats=res,
                 y1="Robustness",
                 extinction.method=extinction.method)



## ************************************************************
## specialization
## ************************************************************
load(file=file.path(save.path, 'specs.Rdata'))
load(file=file.path(save.path, "mods/specs_ypr.Rdata"))

ys <- c("proportional.generality", "d", "degree")
ylabs <- c("Proportional Generality", "Specialization (d')", "Degree")

dd <- expand.grid(ypr=seq(from= min(specs$ypr, na.rm=TRUE),
                          to= max(specs$ypr, na.rm=TRUE),
                          length=10))
                 
pp <- c("plants", "pols")
mods <- list(mod.pols, mod.plants)
names(mods) <- pp

for(j in pp){
  for(i in 1:length(ys)){
    dd.ypr <- cbind(dd, 0)
    colnames(dd.ypr) <- c("ypr", ys[i])

    ypr.pi <- predict.int(mod= mods[[j]][[i]],
                          dd=dd.ypr,
                          y=ys[i],
                          family="gaussian")

    plot.predict.ypr(new.dd=ypr.pi,
                     ylabel=ylabs[i],
                     dats=specs,
                     y1=ys[i],
                     extinction.method=j,
                     agg.col="GenusSpecies")
  }
}


## ************************************************************
## plotting
## ************************************************************
plot.panels()

