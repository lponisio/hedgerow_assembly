library(RColorBrewer)
library(lme4)
library(lmerTest)
## setwd('~/Dropbox/hedgerow_assembly/analysis/networkLevel')
setwd('networkLevel')
source('plotting/src/predictIntervals.R')
source('plotting/src/CIplotting.R')
source('plotting/src/plotPanels.R')
source('src/misc.R')
load(file='saved/corMets.Rdata')
load(file='saved/mods/baci_mods.Rdata')

dd <- expand.grid(ypr=seq(from= min(cor.dats$ypr, na.rm=TRUE),
                          to= max(cor.dats$ypr, na.rm=TRUE),
                          length=10))

## ************************************************************
## nodf
## ************************************************************
dd.nodf <- cbind(dd, zNODF=0)
nodf.pi <- predict.int(mod= mods$zNODF,
                        dd=dd.nodf,
                        y="zNODF",
                       family="gaussian")

## ************************************************************
## modularity
## ************************************************************
dd.mod <- cbind(dd, zmod.met.D=0)
mod.pi <- predict.int(mod= mods$zmod.met.D,
                        dd=dd.mod,
                        y="zmod.met.D",
                      family="gaussian")

## ************************************************************
## specialization
## ************************************************************
dd.h2 <- cbind(dd, zH2=0)
h2.pi <- predict.int(mod= mods$zH2,
                        dd=dd.h2,
                        y="zH2",
                     family="gaussian")

## ************************************************************
## connectance
## ************************************************************
dd.conn <- cbind(dd, connectance =0)
conn.pi <- predict.int(mod= mods$connectance,
                        dd=dd.conn,
                        y="connectance",
                        family="gaussian")

## ************************************************************
## plant richness
## ************************************************************
dd.rich.ll <- cbind(dd, number.of.species.LL =0)
rich.ll.pi <- predict.int(mod= mods.poi$number.of.species.LL,
                        dd=dd.rich.ll,
                        y="number.of.species.LL",
                        family="poisson")



## ************************************************************
## pol richness
## ************************************************************
dd.rich.hl <- cbind(dd, number.of.species.HL =0)
rich.hl.pi <- predict.int(mod= mods.poi$number.of.species.HL,
                        dd=dd.rich.hl,
                        y="number.of.species.HL",
                        family="poisson")

## ************************************************************
## all panels
## ************************************************************
plot.panels()
