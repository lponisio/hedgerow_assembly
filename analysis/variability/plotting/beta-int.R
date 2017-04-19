rm(list=ls())
library(RColorBrewer)
setwd('~/Dropbox/hedgerow_assembly/analysis/variability')
source("plotting/src/predictIntervals.R")
source("plotting/src/CIplotting.R")
source('../networkLevel/src/misc.R', chdir = TRUE)
load('saved/dissMods.Rdata')

plants$SiteStatus <- "all"
pols$SiteStatus <- "all"

## ************************************************************
## interaction dispersion, closeness, plants
## ************************************************************

dd.plants.closeness <- expand.grid(closeness=seq(
                                     from= min(plants$closeness),
                                     to= max(plants$closeness),
                                     length=10),
                                   SiteStatus= c("all"),
                                   Dist= 0)

plants.closeness.pi <- predict.int(mod= mod.plants.close,
                                   dd=dd.plants.closeness,
                                   y="Dist",
                                   family="gaussian")


## ************************************************************
## interaction dispersion, dprime, plants
## ************************************************************

dd.plants.d <- expand.grid(d=seq(
                             from= min(plants$d),
                             to= max(plants$d),
                             length=10),
                           SiteStatus= c("all"),
                           ## SiteStatus= c("control", "maturing", "mature"),
                           Dist= 0)

plants.d.pi <- predict.int(mod= mod.plants.d,
                           dd=dd.plants.d,
                           y="Dist",
                           family="gaussian")

## plot.predict.div(new.dd=plants.d.pi,
##                  ylabel="Interaction dispersion",
##                  dats=plants,
##                  xs="d",
##                  y1="Dist",
##                  xlabel="Plant specialization",
##                  legend.loc="topright",
##                  height=5,
##                  width=5,
##                  x.adj=0.5,
##                  scaled=FALSE,
##                  treatments="all",
##                  col.lines="black",
##                  f.path = 'figures/intdisper')


## ************************************************************
## interaction dispersion, k, pols
## ************************************************************

dd.pols.k <- expand.grid(k=seq(
                           from= min(pols$k),
                           to= max(pols$k),
                           length=10),
                         SiteStatus= c("all"),
                         ## SiteStatus= c("control", "maturing", "mature"),
                         Dist= 0)

pols.k.pi <- predict.int(mod= mod.pols,
                         dd=dd.pols.k,
                         y="Dist",
                         family="gaussian")

## plot.predict.div(new.dd=pols.k.pi,
##                  ylabel="Interaction dispersion",
##                  dats=pols,
##                  xs="k",
##                  y1="Dist",
##                  xlabel="Pollinator network position (k)",
##                  legend.loc="topright",
##                  height=5,
##                  width=5,
##                  x.adj=0.5,
##                  scaled=FALSE,
##                  treatments="all",
##                  col.lines="black",
##                  f.path = 'figures/intdisper')



## ************************************************************
## interaction dispersion, dprime, pols
## ************************************************************

dd.pols.d <- expand.grid(d=seq(
                           from= min(pols$d),
                           to= max(pols$d),
                           length=10),
                         SiteStatus= c("all"),
                         ## SiteStatus= c("control", "maturing", "mature"),
                         Dist= 0)

pols.d.pi <- predict.int(mod= mod.pols.d,
                         dd=dd.pols.d,
                         y="Dist",
                         family="gaussian")

## plot.predict.div(new.dd=pols.d.pi,
##                  ylabel="Interaction dispersion",
##                  dats=pols,
##                  xs="d",
##                  y1="Dist",
##                  xlabel="Pollinator specialization",
##                  legend.loc="topright",
##                  height=5,
##                  width=5,
##                  x.adj=0.5,
##                  scaled=FALSE,
##                  treatments="all",
##                  col.lines="black",
##                  f.path = 'figures/intdisper')




## ************************************************************
## interaction dispersion, closeness, pols
## ************************************************************

dd.pols.closeness <- expand.grid(closeness=seq(
                                   from= min(pols$closeness),
                                   to= max(pols$closeness),
                                   length=10),
                                 SiteStatus= c("all"),
                                 ## SiteStatus= c("control", "maturing", "mature"),
                                 Dist= 0)

pols.closeness.pi <- predict.int(mod= mod.pols.close,
                                 dd=dd.pols.closeness,
                                 y="Dist",
                                 family="gaussian")

## plot.predict.div(new.dd=pols.closeness.pi,
##                  ylabel="Interaction dispersion",
##                  dats=pols,
##                  xs="closeness",
##                  y1="Dist",
##                  xlabel="Pollinator closeness",
##                  legend.loc="topright",
##                  height=5,
##                  width=5,
##                  x.adj=0.5,
##                  scaled=FALSE,
##                  treatments="all",
##                  col.lines="black",
##                  f.path = 'figures/intdisper')

