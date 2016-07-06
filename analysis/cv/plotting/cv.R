rm(list=ls())
setwd('~/Dropbox/hedgerow_assembly/analysis/cv')
source("plotting/src/predictIntervals.R")
source("plotting/src/CIplotting.R")
source('../networkLevel/src/misc.R', chdir = TRUE)
load('saved/contMods.Rdata')

## ************************************************************
## dprime
## ************************************************************

dd.dprime <- expand.grid(traits=seq(
                           from= min(dprime$data$traits),
                           to= max(dprime$data$traits),
                           length=10),
                         status= c("control", "maturing", "mature"),
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
## itd
## ************************************************************

dd.itd <- expand.grid(traits=seq(
                           from= min(itd$data$traits, na.rm=TRUE),
                           to= max(itd$data$traits, na.rm=TRUE),
                           length=10),
                         status= c("control", "maturing", "mature"),
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
