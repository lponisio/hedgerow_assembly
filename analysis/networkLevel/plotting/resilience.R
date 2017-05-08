rm(list=ls())
## setwd('~/Dropbox/hedgerow_assembly/analysis/networkLevel')
setwd('analysis/networkLevel')
source('plotting/src/predictIntervals.R')
source('plotting/src/CIplotting.R')
source('plotting/src/plotPanels_resilence.R')
source('src/initialize.R')

## ************************************************************
## robustness to species extinction
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

## ************************************************************
## robustness to perturbation
## ************************************************************

load(file=file.path(save.path, 'mods/AlgCon.Rdata'))

dd.ypr.alg <- expand.grid(ypr=seq(from= min(all.alg.Con.status$ypr,
                                    na.rm=TRUE),
                    to= max(all.alg.Con.status$ypr,
                      na.rm=TRUE),
                    length=10),
                  AlgCon=0)

ypr.pi.alg <- predict.int(mod= alg.con.mod.ypr,
                        dd=dd.ypr.alg,
                        y="AlgCon",
                        family="gaussian")

plot.panels()


