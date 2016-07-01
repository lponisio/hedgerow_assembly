rm(list=ls())
setwd('~/Dropbox/hedgerow_assembly/analysis/networkLevel')
source('plotting/src/predictIntervals.R')
source('plotting/src/CIplotting.R')
source('plotting/src/plotPanels.R')
source('src/initialize.R')

## ************************************************************
## robustness
## ************************************************************
## either "abund" or "degree"
extinction.method <- "degree"

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

