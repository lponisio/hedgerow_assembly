rm(list=ls())
setwd('~/Dropbox/hedgerow_assembly/analysis/cv')
source('src/initialize.R')

## ************************************************************
## coefficient of variation of abundance through time
## ************************************************************
byYear <- aggregate(list(Abund=spec$GenusSpecies),
                    list(GenusSpecies= spec$GenusSpecies,
                         Date=spec$Date,
                         SiteStatus= spec$SiteStatus,
                         Site=spec$Site), length)

dprime <- cv.trait(spec, byYear, trait="d",
                   method= "time", time.col="Date",
                   abund.col="Abund")

itd <- cv.trait(spec, byYear, trait="ITD",
                method= "time",  time.col="Date",
                abund.col="Abund")

## ************************************************************
## coefficient of variation of degree thingy through time
## ************************************************************

dprime.k <- cv.trait(spec, specs, trait="d",
                     method= "time", time.col="assem",
                     abund.col="k",
                     cv.function=corCv, zero2na=TRUE)


dprime.k.sd <- cv.trait(spec, specs, trait="d",
                     method= "time", time.col="assem",
                     abund.col="k",
                     cv.function=sd,  zero2na=TRUE, na.rm=TRUE)

dprime.closeness <- cv.trait(spec, specs, trait="d",
                     method= "time", time.col="assem",
                     abund.col="weighted.closeness",
                     cv.function=corCv, zero2na=TRUE)

dprime.closeness.sd <- cv.trait(spec, specs, trait="d",
                     method= "time", time.col="assem",
                     abund.col="weighted.closeness",
                     cv.function=sd)

save(itd, dprime, dprime.k.sd, dprime.closeness.sd,
     file="saved/contMods.Rdata")


## ************************************************************
## coefficient of variation through space
## ************************************************************

## byYr <- aggregate(list(Abund=spec$GenusSpecies),
##                   list(GenSp= spec$GenusSpecies,
##                        status= spec$SiteStatus,
##                        date= spec$Site,
##                        site= spec$Year), length)

## dprime.sp <- cv.trait(spec, byYr, trait="d", xlabel= "Specialization",
##                       method= "space")
## itd.sp <- cv.trait(spec, byYr, trait="ITD", xlabel= "Body size",
##                    method= "space")
## lecty.sp <-  cv.trait(spec, byYr, trait="Lecty", cont=FALSE,
##                       method= "space")
## excavate.sp <-  cv.trait(spec, byYr, trait="Excavate", cont=FALSE,
##                          method=" space")
## nest.sp <-  cv.trait(spec, byYr, trait="NestLoc", cont=FALSE,
##                      method= "space")
## soc.sp <- cv.trait(spec, byYr, trait="Sociality", cont=FALSE,
##                   method= "space")
