rm(list=ls())
setwd('~/Dropbox/hedgerow_assembly/analysis/variability')
source('src/initialize.R')

## ************************************************************
## occurrence
## ************************************************************
## pollinators and closeness

pol.cv <- cv.trait(spec,
                   specs[specs$speciesType =="pollinator",],
                   trait1="occ.date",
                   trait2="degree",
                   method= "time", time.col="assem",
                   abund.col="weighted.closeness",
                   cv.function=cv,
                   zero2na=TRUE,
                   standard.cv=TRUE,
                   na.rm=TRUE)

pol.mod <- lme4::lmer(formula.cv, data=pol.cv$lm.data)
summary(pol.mod)
vif.mer(pol.mod)

## the reviewers wanted a quantile regression, can only include on
## random effect at a time
pol.lm.cv.quant <- lqmm(fixed=cv ~ occ.date + degree, random=~1,
                        group= GenusSpecies,
                        data=pol.cv$lm.data,
                        control=list(LP_max_iter=10^4))
pol.sum.boot.quant <- summary.boot.lqmm(boot(pol.lm.cv.quant, R=100))


## plants and closeness
plants.cv <- cv.trait(spec,
                      specs[specs$speciesType =="plant",],
                      trait1="occ.plant.date",
                      trait2="plant.degree",
                      method= "time", time.col="assem",
                      abund.col="weighted.closeness",
                      cv.function=cv,
                      zero2na=TRUE,
                      standard.cv=TRUE,
                      na.rm=TRUE,
                      species.type="PlantGenusSpecies")

plants.mod <- lmer(formula.plant.cv, data=plants.cv$lm.data)

summary(plants.mod)
vif.mer(plants.mod)

## the reviewers wanted a quantile regression, can only include on
## random effect at a time
plants.lm.cv.quant <- lqmm(fixed=cv ~ occ.plant.date + plant.degree, random=~1,
                           group= GenusSpecies,
                           data=plants.cv$lm.data,
                           control=list(LP_max_iter=10^3))
plants.sum.boot.quant <- summary.boot.lqmm(boot(plants.lm.cv.quant, R=100))


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

