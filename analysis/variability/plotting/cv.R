rm(list=ls())
library(RColorBrewer)
setwd('~/Dropbox/hedgerow_assembly/analysis/variability')
source("plotting/src/predictIntervals.R")
source("plotting/src/CIplotting.R")
source('../networkLevel/src/misc.R', chdir = TRUE)
load('saved/contMods.Rdata')

## ************************************************************
## dprime abundance
## ************************************************************

dd.dprime <- expand.grid(traits=seq(
                           from= min(dprime$data$traits),
                           to= max(dprime$data$traits),
                           length=10),
                         SiteStatus= c("control", "maturing", "mature"),
                         cv= 0)

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
                 x.adj=0.5)


## ************************************************************
## dprime network position - k
## ************************************************************

dprime.pi <- predict.int(mod= dprime.k.sd$lm,
                        dd=dd.dprime,
                        y="cv",
                        family="gaussian")

plot.predict.div(new.dd=dprime.pi,
                 ylabel="Network position variability",
                 dats=dprime.k.sd$data,
                 xs="traits",
                 y1="cv",
                 xlabel="Specialization",
                 legend.loc="topright",
                 height=5,
                 width=5,
                 x.adj=0.5)



## ************************************************************
## dprime network position - closeness
## ************************************************************

dprime.pi.cl <- predict.int(mod= dprime.closeness.sd$lm,
                        dd=dd.dprime,
                        y="cv",
                        family="gaussian")

plot.predict.div(new.dd=dprime.pi.cl,
                 ylabel="Closeness variability",
                 dats=dprime.closeness.sd$data,
                 xs="traits",
                 y1="cv",
                 xlabel="Specialization",
                 legend.loc="topright",
                 height=5,
                 width=5,
                 x.adj=0.5)


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
                 x.adj=0.5)

## ************************************************************
## persistence network position
## ************************************************************
dd.occ <- expand.grid(traits=seq(
                           from= min(occ.k.sd$data$traits, na.rm=TRUE),
                           to= max(occ.k.sd$data$traits, na.rm=TRUE),
                           length=10),
                         SiteStatus= c("control", "maturing", "mature"),
                         cv= 0)

occ.pi <- predict.int(mod= occ.k.sd$lm,
                        dd=dd.occ,
                        y="cv",
                        family="gaussian")

plot.predict.div(new.dd=occ.pi,
                 ylabel="Network position variability",
                 dats=occ.k.sd$data,
                 xs="traits",
                 y1="cv",
                 xlabel="Persistence",
                 legend.loc="topleft",
                 height=5,
                 width=5,
                 x.adj=0.5)

