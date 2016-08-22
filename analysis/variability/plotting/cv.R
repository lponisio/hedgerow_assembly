rm(list=ls())
library(RColorBrewer)
setwd('~/Dropbox/hedgerow_assembly/analysis/variability')
source("plotting/src/predictIntervals.R")
source("plotting/src/CIplotting.R")
source("plotting/src/plotPanels.R")
source('../networkLevel/src/misc.R', chdir = TRUE)
load('saved/contMods.Rdata')

## ************************************************************
## persistence k
## ************************************************************
occ.k.cv$data$SiteStatus <- "all"
dd.occ <- expand.grid(traits=seq(
                           from= min(occ.k.cv$data$traits, na.rm=TRUE),
                           to= max(occ.k.cv$data$traits, na.rm=TRUE),
                           length=10),
                      SiteStatus="all",
                         cv= 0)

occ.pi.k <- predict.int(mod= occ.k.cv$lm.nss,
                        dd=dd.occ,
                        y="cv",
                        family="gaussian")

plot.predict.div(new.dd=occ.pi.k,
                 ylabel="Network position variability (k)",
                 dats=occ.k.cv$data, 
                 xs="traits",
                 y1="cv",
                 xlabel="Pollinator persistence",
                 legend.loc="bottomright",
                 height=5,
                 width=5,
                 x.adj=0.5,
                 treatments="all",
                 col.lines="black",
                 f.path='figures/cv')

## ************************************************************
## persistence closeness
## ************************************************************
occ.closeness.cv$data$SiteStatus <- "all"
occ.pi.close <- predict.int(mod= occ.closeness.cv$lm.nss,
                        dd=dd.occ,
                        y="cv",
                        family="gaussian")

plot.predict.div(new.dd=occ.pi.close,
                 ylabel="Network closeness variability",
                 dats=occ.closeness.cv$data,
                 xs="traits",
                 y1="cv",
                 xlabel="Pollinator persistence",
                 legend.loc="bottomright",
                 height=5,
                 width=5,
                 x.adj=0.5,
                 treatments="all",
                 col.lines="black",
                 f.path='figures/cv')

## ************************************************************
## degree k
## ************************************************************
degree.k.cv$data$SiteStatus <- "all"
dd.degree <- expand.grid(traits=seq(
                           from= min(degree.k.cv$data$traits, na.rm=TRUE),
                           to= max(degree.k.cv$data$traits, na.rm=TRUE),
                           length=10),
                         ## SiteStatus= c("control", "maturing",
                      ## "mature"),
                      SiteStatus="all",
                         cv= 0)

degree.pi <- predict.int(mod= degree.k.cv$lm.nss,
                        dd=dd.degree,
                        y="cv",
                        family="gaussian")

plot.predict.div(new.dd=degree.pi,
                 ylabel="Network position variability (k)",
                 dats=degree.k.cv$data, 
                 xs="traits",
                 y1="cv",
                 xlabel="Pollinator degree",
                 legend.loc="bottomright",
                 height=5,
                 width=5,
                 x.adj=0.5,
                 treatments="all",
                 col.lines="black",
                 f.path='figures/cv')

## ************************************************************
## degree closeness
## ************************************************************
degree.closeness.cv$data$SiteStatus <- "all"
degree.pi <- predict.int(mod= degree.closeness.cv$lm.nss,
                        dd=dd.degree,
                        y="cv",
                        family="gaussian")

## plot.predict.div(new.dd=degree.pi,
##                  ylabel="Network closeness variability",
##                  dats=degree.closeness.cv$data,
##                  xs="traits",
##                  y1="cv",
##                  xlabel="Pollinator degree",
##                  legend.loc="bottomright",
##                  height=5,
##                  width=5,
##                  x.adj=0.5,
##                  treatments="all",
##                  col.lines="black",
##                  f.path='figures/cv',
##                  dec=0)

plot.panels()




## ************************************************************
## dprime network position - k
## ************************************************************

dd.dprime <- expand.grid(traits=seq(
                           from= min(dprime$data$traits),
                           to= max(dprime$data$traits),
                           length=10),
                         SiteStatus= c("control", "maturing", "mature"),
                         cv= 0)

dprime.pi <- predict.int(mod= dprime.k.cv$lm,
                        dd=dd.dprime,
                        y="cv",
                        family="gaussian")

plot.predict.div(new.dd=dprime.pi,
                 ylabel="Network position variability",
                 dats=dprime.k.cv$data,
                 xs="traits",
                 y1="cv",
                 xlabel="Specialization",
                 legend.loc="topright",
                 height=5,
                 width=5,
                 x.adj=0.5,
                 f.path='figures/cv')

## ************************************************************
## dprime network position - closeness
## ************************************************************

dprime.pi.cl <- predict.int(mod= dprime.closeness.cv$lm,
                        dd=dd.dprime,
                        y="cv",
                        family="gaussian")

plot.predict.div(new.dd=dprime.pi.cl,
                 ylabel="Closeness variability",
                 dats=dprime.closeness.cv$data,
                 xs="traits",
                 y1="cv",
                 xlabel="Specialization",
                 legend.loc="topright",
                 height=5,
                 width=5,
                 x.adj=0.5,
                 f.path='figures/cv')


## ************************************************************
## itd abundance
## ************************************************************

dd.itd <- expand.grid(traits=seq(
                           from= min(itd$data$traits, na.rm=TRUE),
                           to= max(itd$data$traits, na.rm=TRUE),
                           length=10),
                         SiteStatus= c("control", "maturing", "mature"),
                         cv= 0)

itd.pi <- predict.int(mod= itd$lm,
                        dd=dd.itd,
                        y="cv",
                        family="gaussian")

plot.predict.div(new.dd=itd.pi,
                 ylabel="Coefficient of variation",
                 dats=itd$data,
                 xs="traits",
                 y1="cv",
                 xlabel="Body size",
                 legend.loc="bottomright",
                 height=5,
                 width=5,
                 x.adj=0.5,
                 f.path='figures/cv')




## ************************************************************
## dprime abundance
## ************************************************************

dprime.pi <- predict.int(mod= dprime$lm,
                        dd=dd.dprime,
                        y="cv",
                        family="gaussian")

plot.predict.div(new.dd=dprime.pi,
                 ylabel="Coefficient of variation",
                 dats=dprime$data,
                 xs="traits",
                 y1="cv",
                 xlabel="Specialization",
                 legend.loc="bottomright",
                 height=5,
                 width=5,
                 x.adj=0.5,
                 f.path='figures/cv')
