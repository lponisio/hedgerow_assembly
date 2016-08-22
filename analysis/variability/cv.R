rm(list=ls())
setwd('~/Dropbox/hedgerow_assembly/analysis/variability')
source('src/initialize.R')

## ************************************************************
## coefficient of variation of degree thingy through time
## ************************************************************
## ************************************************************
## occurrence
## ************************************************************
## pollinators and k
## not sig
occ.k.cv <- cv.trait(spec,
                     specs[specs$speciesType =="pollinator",],
                     trait="occ.date",
                     method= "time", time.col="assem",
                     abund.col="k",
                     cv.function=cv,
                     zero2na=TRUE,
                     standard.cv=FALSE,
                     na.rm=TRUE)
summary(occ.k.cv$lm.nss)

## plants and k
## medium sig
plants.occ.k.cv <- cv.trait(spec,
                            specs[specs$speciesType =="plant",],
                            trait="occ.plant.date",
                            method= "time", time.col="assem",
                            abund.col="k",
                            cv.function=cv,
                            zero2na=TRUE,
                            standard.cv=FALSE,
                            na.rm=TRUE,
                            species.type="PlantGenusSpecies")
summary(plants.occ.k.cv$lm.nss)

## pollinators and closeness
## medium sig
occ.closeness.cv <- cv.trait(spec,
                             specs[specs$speciesType =="pollinator",],
                             trait="occ.date",
                             method= "time", time.col="assem",
                             abund.col="weighted.closeness",
                             cv.function=cv,
                             zero2na=TRUE,
                             standard.cv=FALSE,
                             na.rm=TRUE)
summary(occ.closeness.cv$lm.nss)

## plants and closeness
## not sig
plants.occ.closeness.cv <- cv.trait(spec,
                                    specs[specs$speciesType =="plant",],
                                    trait="occ.plant.date",
                                    method= "time", time.col="assem",
                                    abund.col="weighted.closeness",
                                    cv.function=cv,
                                    zero2na=TRUE,
                                    standard.cv=FALSE,
                                    na.rm=TRUE,
                                    species.type="PlantGenusSpecies")
summary(plants.occ.closeness.cv$lm.nss)

## ************************************************************
## degree
## ************************************************************
## pollinators and k
## not sig
degree.k.cv <- cv.trait(spec,
                        specs[specs$speciesType =="pollinator",],
                        trait="degree",
                        method= "time", time.col="assem",
                        abund.col="k",
                        cv.function=cv,
                        zero2na=TRUE,
                        standard.cv=FALSE,
                        na.rm=TRUE)
summary(degree.k.cv$lm.nss)

## plants and k
## not sig
plants.degree.k.cv <- cv.trait(spec,
                               specs[specs$speciesType =="plant",],
                               trait="plant.degree",
                               method= "time", time.col="assem",
                               abund.col="k",
                               cv.function=cv,
                               zero2na=TRUE,
                               standard.cv=FALSE,
                               na.rm=TRUE,
                               species.type="PlantGenusSpecies")
summary(plants.degree.k.cv$lm.nss)

## pollinators and closeness
## sig
degree.closeness.cv <- cv.trait(spec,
                                specs[specs$speciesType =="pollinator",],
                                trait="degree",
                                method= "time", time.col="assem",
                                abund.col="weighted.closeness",
                                cv.function=cv,
                                zero2na=TRUE,
                                standard.cv=FALSE,
                                na.rm=TRUE)
summary(degree.closeness.cv$lm.nss)

## plants and closeness
## not sig!
plants.degree.closeness.cv <- cv.trait(spec,
                                       specs[specs$speciesType =="plant",],
                                       trait="plant.degree",
                                       method= "time", time.col="assem",
                                       abund.col="weighted.closeness",
                                       cv.function=cv,
                                       zero2na=TRUE,
                                       standard.cv=FALSE,
                                       na.rm=TRUE,
                                       species.type="PlantGenusSpecies")
summary(plants.degree.closeness.cv$lm.nss)


## check correlation of degree and occ 
occ.k.sd$data$degree <-
  degree$data$traits.ns[match(occ.k.sd$data$GenusSpecies,
                              degree$data$GenusSpecies)]

plot(occ.k.sd$data$traits.ns ~ occ.k.sd$data$degree)

cor.test(occ.k.sd$data$traits.ns, occ.k.sd$data$degree)


## ************************************************************
## dprime
## ************************************************************
## not sig
dprime.k.cv <- cv.trait(spec,
                        specs[specs$speciesType =="pollinator",],
                        trait="d",
                        method= "time", time.col="assem",
                        abund.col="k",
                        cv.function=cv,
                        zero2na=TRUE,
                        standard.cv=TRUE,
                        na.rm=TRUE)
summary(dprime.k.cv$lm.nss)

## not sig
dprime.closeness.cv <- cv.trait(spec,
                                specs[specs$speciesType =="pollinator",],
                                trait="d",
                                method= "time", time.col="assem",
                                abund.col="weighted.closeness",
                                cv.function=cv,
                                zero2na=TRUE,
                                standard.cv=FALSE,
                                na.rm=TRUE)
summary(dprime.closeness.cv$lm.nss)


## ************************************************************
## coefficient of variation of abundance through time
## ************************************************************
byYear <- aggregate(list(Abund=spec$GenusSpecies),
                    list(GenusSpecies= spec$GenusSpecies,
                         Date=spec$Date,
                         SiteStatus= spec$SiteStatus,
                         Site=spec$Site), length)

dprime <- cv.trait(spec, byYear, trait="d",
                   method= "time",
                   time.col="Date",
                   abund.col="Abund")

degree <- cv.trait(spec, byYear, trait="degree",
                   method= "time",
                   time.col="Date",
                   abund.col="Abund")

itd <- cv.trait(spec, byYear, trait="ITD",
                method= "time",
                time.col="Date",
                abund.col="Abund")

## ************************************************************
## save
save(itd, dprime, degree,
     dprime.k.cv, dprime.closeness.cv,
     occ.k.cv, occ.closeness.cv,
     plants.occ.k.cv, plants.occ.closeness.cv,
     degree.k.cv, degree.closeness.cv,
     plants.degree.k.cv, plants.degree.closeness.cv,
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
