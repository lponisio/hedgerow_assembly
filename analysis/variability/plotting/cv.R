rm(list=ls())
library(RColorBrewer)
setwd('~/Dropbox/hedgerow_assembly/analysis/variability')
source("plotting/src/predictIntervals.R")
source("plotting/src/CIplotting.R")
source("plotting/src/plotPanels.R")
source('../networkLevel/src/misc.R', chdir = TRUE)
load('saved/contMods.Rdata')

## ************************************************************
## persistence closeness
## ************************************************************
dd.occ.pol <- expand.grid(traits=seq(
                            from= min(occ.closeness.cv$data$traits,
                              na.rm=TRUE),
                            to= max(occ.closeness.cv$data$traits,
                              na.rm=TRUE),
                            length=20),
                          SiteStatus="all",
                          cv= 0)

## pols
occ.closeness.cv$data$SiteStatus <- "all"
occ.pi.close <- predict.int(mod= occ.closeness.cv$lm.nss,
                            dd=dd.occ.pol,
                            y="cv",
                            family="gaussian")

## plants
dd.occ.plants <- expand.grid(traits=seq(
                               from= min(plants.occ.closeness.cv$data$traits,
                                 na.rm=TRUE),
                               to= max(plants.occ.closeness.cv$data$traits,
                                 na.rm=TRUE),
                               length=20),
                             SiteStatus="all",
                             cv= 0)


plants.occ.closeness.cv$data$SiteStatus <- "all"
plants.occ.pi.close <- predict.int(mod= plants.occ.closeness.cv$lm.nss,
                                   dd=dd.occ.plants,
                                   y="cv",
                                   family="gaussian")


## ************************************************************
## degree closeness
## ************************************************************
dd.degree.pol <- expand.grid(traits=seq(
                               from=
                               min(degree.closeness.cv$data$traits,
                                   na.rm=TRUE),
                               to=
                               max(degree.closeness.cv$data$traits,
                                   na.rm=TRUE),
                               length=10),
                             SiteStatus="all",
                             cv= 0)
## pols

degree.closeness.cv$data$SiteStatus <- "all"

degree.pi <- predict.int(mod= degree.closeness.cv$lm.nss,
                         dd=dd.degree.pol,
                         y="cv",
                         family="gaussian")

## plants
dd.degree.plants <- expand.grid(traits=seq(
                           from=
                           min(plants.degree.closeness.cv$data$traits,
                               na.rm=TRUE),
                           to=
                           max(plants.degree.closeness.cv$data$traits,
                               na.rm=TRUE),
                           length=10),
                         SiteStatus="all",
                         cv= 0)

plants.degree.closeness.cv$data$SiteStatus <- "all"

plants.degree.pi <- predict.int(mod= plants.degree.closeness.cv$lm.nss,
                                dd=dd.degree.plants,
                                y="cv",
                                family="gaussian")



plot.panels()




## ## ************************************************************
## ## dprime network position - k
## ## ************************************************************

## dd.dprime <- expand.grid(traits=seq(
##                            from= min(dprime$data$traits),
##                            to= max(dprime$data$traits),
##                            length=10),
##                          SiteStatus= c("control", "maturing", "mature"),
##                          cv= 0)

## dprime.pi <- predict.int(mod= dprime.k.cv$lm,
##                          dd=dd.dprime,
##                          y="cv",
##                          family="gaussian")

## plot.predict.div(new.dd=dprime.pi,
##                  ylabel="Network position variability",
##                  dats=dprime.k.cv$data,
##                  xs="traits",
##                  y1="cv",
##                  xlabel="Specialization",
##                  legend.loc="topright",
##                  height=5,
##                  width=5,
##                  x.adj=0.5,
##                  f.path='figures/cv')

## ## ************************************************************
## ## dprime network position - closeness
## ## ************************************************************

## dprime.pi.cl <- predict.int(mod= dprime.closeness.cv$lm,
##                             dd=dd.dprime,
##                             y="cv",
##                             family="gaussian")

## plot.predict.div(new.dd=dprime.pi.cl,
##                  ylabel="Closeness variability",
##                  dats=dprime.closeness.cv$data,
##                  xs="traits",
##                  y1="cv",
##                  xlabel="Specialization",
##                  legend.loc="topright",
##                  height=5,
##                  width=5,
##                  x.adj=0.5,
##                  f.path='figures/cv')


## ## ************************************************************
## ## itd abundance
## ## ************************************************************

## dd.itd <- expand.grid(traits=seq(
##                         from= min(itd$data$traits, na.rm=TRUE),
##                         to= max(itd$data$traits, na.rm=TRUE),
##                         length=10),
##                       SiteStatus= c("control", "maturing", "mature"),
##                       cv= 0)

## itd.pi <- predict.int(mod= itd$lm,
##                       dd=dd.itd,
##                       y="cv",
##                       family="gaussian")

## plot.predict.div(new.dd=itd.pi,
##                  ylabel="Coefficient of variation",
##                  dats=itd$data,
##                  xs="traits",
##                  y1="cv",
##                  xlabel="Body size",
##                  legend.loc="bottomright",
##                  height=5,
##                  width=5,
##                  x.adj=0.5,
##                  f.path='figures/cv')




## ## ************************************************************
## ## dprime abundance
## ## ************************************************************

## dprime.pi <- predict.int(mod= dprime$lm,
##                          dd=dd.dprime,
##                          y="cv",
##                          family="gaussian")

## plot.predict.div(new.dd=dprime.pi,
##                  ylabel="Coefficient of variation",
##                  dats=dprime$data,
##                  xs="traits",
##                  y1="cv",
##                  xlabel="Specialization",
##                  legend.loc="bottomright",
##                  height=5,
##                  width=5,
##                  x.adj=0.5,
##                  f.path='figures/cv')
