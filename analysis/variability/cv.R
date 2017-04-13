rm(list=ls())
setwd('~/Dropbox/hedgerow_assembly/analysis/variability')
source('src/initialize.R')

## ************************************************************
## occurrence
## ************************************************************
## pollinators and closeness
## sig
occ.closeness.cv <- cv.trait(spec,
                             specs[specs$speciesType =="pollinator",],
                             trait="occ.date",
                             method= "time", time.col="assem",
                             abund.col="weighted.closeness",
                             cv.function=cv,
                             zero2na=TRUE,
                             standard.cv=TRUE,
                             na.rm=TRUE)

summary(occ.closeness.cv$lm.nss)
occ.closeness.cv$sum.quant


## plants and closeness
## not sig
plants.occ.closeness.cv <- cv.trait(spec,
                                    specs[specs$speciesType =="plant",],
                                    trait="occ.plant.date",
                                    method= "time", time.col="assem",
                                    abund.col="weighted.closeness",
                                    cv.function=cv,
                                    zero2na=TRUE,
                                    standard.cv=TRUE,
                                    na.rm=TRUE,
                                    species.type="PlantGenusSpecies")
summary(plants.occ.closeness.cv$lm.nss)
plants.occ.closeness.cv$sum.quant

## ************************************************************
## degree
## ************************************************************
## pollinators and closeness
## sig
degree.closeness.cv <- cv.trait(spec,
                                specs[specs$speciesType =="pollinator",],
                                trait="degree",
                                method= "time", time.col="assem",
                                abund.col="weighted.closeness",
                                cv.function=cv,
                                zero2na=TRUE,
                                standard.cv=TRUE,
                                na.rm=TRUE)
summary(degree.closeness.cv$lm.nss)
degree.closeness.cv$sum.quant

## plants and closeness
## not sig!
plants.degree.closeness.cv <- cv.trait(spec,
                                       specs[specs$speciesType =="plant",],
                                       trait="plant.degree",
                                       method= "time", time.col="assem",
                                       abund.col="weighted.closeness",
                                       cv.function=cv,
                                       zero2na=TRUE,
                                       standard.cv=TRUE,
                                       na.rm=TRUE,
                                       species.type="PlantGenusSpecies")
summary(plants.degree.closeness.cv$lm.nss)
plants.degree.closeness.cv$sum.quant

## check correlation of degree and occ
## pollinators
layout(matrix(1:2, nrow=1))
check.pol <- unique(cbind(spec$degree,
                          spec$occ.date))
plot(check.pol)

cor.test(check.pol[,1], check.pol[,2])


check.plant <- unique(cbind(spec$plant.degree,
                            spec$occ.plant.date))
plot(check.plant)

cor.test(check.plant[,1], check.plant[,2])


## ************************************************************
## save
save(occ.closeness.cv,
    plants.occ.closeness.cv,
   degree.closeness.cv, plants.degree.closeness.cv,
     file="saved/contMods.Rdata")
