rm(list=ls())
library(lme4)
library(vegan)
library(fields)
library(raster)
library(lmerTest)
setwd('~/Dropbox/hedgerow_assembly/analysis/cv')
source('src/misc.R')
source('src/cvCalc.R')
load('../../data/networks/allSpecimens.Rdata')

## ************************************************************
## coefficient of variation through time
## ************************************************************
byDate <- aggregate(list(Abund=spec$GenusSpecies),
                    list(GenSp= spec$GenusSpecies,
                         date=spec$Date,
                         status= spec$SiteStatus,
                         site=spec$Site), length)

dprime <- cv.trait(spec, byDate, trait="d",
                   method= "time")

itd <- cv.trait(spec, byDate, trait="ITD",
                method= "time")

save(itd, dprime, file="saved/contMods.Rdata")


## discrete traits

## lecty <-  cv.trait(spec, byDate, trait="Lecty", cont=FALSE,
##                    method= "time")
## excavate <-  cv.trait(spec, byDate, trait="Excavate", cont=FALSE,
##                       method=" time")
## nest <-  cv.trait(spec, byDate, trait="NestLoc", cont=FALSE,
##                   method= "time")
## soc <- cv.trait(spec, byDate, trait="Sociality", cont=FALSE,
##                   method= "time")


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
