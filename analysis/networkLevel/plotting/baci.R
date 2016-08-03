library(RColorBrewer)
library(lme4)
library(lmerTest)
source('plotting/src/predictIntervals.R')
source('plotting/src/CIplotting.R')
source('plotting/src/plotPanels.R')
source('src/misc.R')
load(file='saved/corMets.Rdata')

dd <- expand.grid(ypr=seq(from= min(cor.dats$ypr, na.rm=TRUE),
                          to= max(cor.dats$ypr, na.rm=TRUE),
                          length=10))

## ************************************************************
## nodf
## ************************************************************
load(file='saved/mods/baci_nodf.Rdata')

dd.nodf <- cbind(dd, zNODF=0)

nodf.pi <- predict.int(mod= baci.nodf.mod,
                        dd=dd.nodf,
                        y="zNODF",
                       family="gaussian")

## ************************************************************
## modularity
## ************************************************************
load(file='saved/mods/baci_mod.Rdata')

dd.mod <- cbind(dd, zmodD=0)

mod.pi <- predict.int(mod= baci.mod.mod,
                        dd=dd.mod,
                        y="zmodD",
                      family="gaussian")

## ************************************************************
## specialization
## ************************************************************
load(file='saved/mods/baci_h2.Rdata')

dd.h2 <- cbind(dd, zH2=0)

h2.pi <- predict.int(mod= baci.h2.mod,
                        dd=dd.h2,
                        y="zH2",
                        family="gaussian")

## ************************************************************
## niche overlap - Plants
## ************************************************************
load(file='saved/mods/baci_no_plant.Rdata')

dd.nop <- cbind(dd, niche.overlap.plants =0)

nop.pi <- predict.int(mod= baci.no.plant.mod,
                        dd=dd.nop,
                        y="niche.overlap.plants",
                        family="gaussian")
                        
## ************************************************************
## niche overlap - Polinators
## ************************************************************
load(file='saved/mods/baci_no_pol.Rdata')

dd.nopol <- cbind(dd, niche.overlap.pol =0)
nopol.pi <- predict.int(mod= baci.no.pol.mod,
                        dd=dd.nopol,
                        y="niche.overlap.pol",
                        family="gaussian")                        
## ************************************************************
## all panels
## ************************************************************
plot.panels()
