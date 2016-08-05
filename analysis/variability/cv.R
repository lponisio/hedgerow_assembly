rm(list=ls())
setwd('~/Dropbox/hedgerow_assembly/analysis/variability')
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
## coefficient of variation of degree thingy through time
## ************************************************************

## ************************************************************
## dprime
## not sig
dprime.k.sd <- cv.trait(spec,
                        specs[specs$speciesType =="pollinator",],
                        trait="d",
                        method= "time", time.col="assem",
                        abund.col="k",
                        cv.function=sd,
                        zero2na=TRUE, standard.cv=FALSE,
                        na.rm=TRUE)
summary(dprime.k.sd$lm.nss)

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
## occurrence
## sig
occ.k.sd <- cv.trait(spec,
                     specs[specs$speciesType =="pollinator",],
                     trait="occ.date",
                     method= "time", time.col="assem",
                     abund.col="k",
                     cv.function=sd,
                     zero2na=TRUE,
                     standard.cv=FALSE,
                     na.rm=TRUE)
summary(occ.k.sd$lm.nss)

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


## check correlation of dprime and occ
occ.k.sd$data$spec <-
  dprime$data$traits.ns[match(occ.k.sd$data$GenusSpecies,
                              dprime$data$GenusSpecies)]

plot(occ.k.sd$data$traits.ns ~ occ.k.sd$data$spec)

cor.test(occ.k.sd$data$traits.ns, occ.k.sd$data$spec)

## ************************************************************
## degree
## sig
degree.k.sd <- cv.trait(spec,
                        specs[specs$speciesType =="pollinator",],
                        trait="degree",
                        method= "time",
                        time.col="assem",
                        abund.col="k",
                        cv.function=sd,
                        zero2na=TRUE,
                        standard.cv=FALSE,
                        na.rm=TRUE)
summary(degree.k.sd$lm.nss)

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


## check correlation of degree and occ 
occ.k.sd$data$degree <-
  degree$data$traits.ns[match(occ.k.sd$data$GenusSpecies,
                              degree$data$GenusSpecies)]

plot(occ.k.sd$data$traits.ns ~ occ.k.sd$data$degree)

cor.test(occ.k.sd$data$traits.ns, occ.k.sd$data$degree)



## ************************************************************
## save
save(itd, dprime, degree,
     dprime.k.sd, dprime.k.cv, dprime.closeness.cv,
     occ.k.sd, occ.k.cv, occ.closeness.cv,
     degree.k.sd, degree.k.cv, degree.closeness.cv,
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
