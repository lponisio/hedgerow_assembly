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

dd.mod <- cbind(dd, zmod.met.D=0)

mod.pi <- predict.int(mod= baci.mod.mod,
                        dd=dd.mod,
                        y="zmod.met.D",
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

dd.nop <- cbind(dd, niche.overlap.LL =0)

nop.pi <- predict.int(mod= baci.no.plant.mod,
                        dd=dd.nop,
                        y="niche.overlap.plants",
                        family="gaussian")
                        
## ************************************************************
## niche overlap - Polinators
## ************************************************************
load(file='saved/mods/baci_no_pol.Rdata')

dd.nopol <- cbind(dd, niche.overlap.HL =0)
nopol.pi <- predict.int(mod= baci.no.pol.mod,
                        dd=dd.nopol,
                        y="niche.overlap.pol",
                        family="gaussian")

## ************************************************************
## connectance
## ************************************************************
load(file='saved/mods/baci_conn.Rdata')

dd.conn <- cbind(dd, connectance =0)
conn.pi <- predict.int(mod= baci.conn.mod,
                        dd=dd.conn,
                        y="connectance",
                        family="gaussian")

## ************************************************************
## plant richness
## ************************************************************
load(file='saved/mods/baci_rich_ll.Rdata')

dd.rich.ll <- cbind(dd, number.of.species.LL =0)
rich.ll.pi <- predict.int(mod= baci.rich.ll.mod,
                        dd=dd.rich.ll,
                        y="number.of.species.LL",
                        family="poisson")



## ************************************************************
## pol richness
## ************************************************************
load(file='saved/mods/baci_rich_hl.Rdata')

dd.rich.hl <- cbind(dd, number.of.species.HL =0)
rich.hl.pi <- predict.int(mod= baci.rich.hl.mod,
                        dd=dd.rich.hl,
                        y="number.of.species.HL",
                        family="poisson")



## ************************************************************
## all panels
## ************************************************************
plot.panels()
