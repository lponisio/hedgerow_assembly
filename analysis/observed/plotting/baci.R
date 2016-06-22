library(RColorBrewer)
source('plotting/src/predictIntervals.R')
source('plotting/src/CIplotting.R')
source('plotting/src/plotPanels.R')
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
## all panels
## ************************************************************
plot.panels()

