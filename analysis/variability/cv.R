rm(list=ls())
setwd('~/Dropbox/hedgerow_assembly/analysis/variability')
source('src/initialize.R')

## ************************************************************
## pollinators
## ************************************************************

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

pol.mod <- lmer(formula.cv, data=pol.cv$lm.data)
summary(pol.mod)
vif.mer(pol.mod)

## variance inflation factors > 2, so not great (Zurr et al. 2010)

## models with each explanatory variable
pol.mod.degree <- lmer(cv ~ degree + (1|Site) +
                           (1|GenusSpecies),
                       data=pol.cv$lm.data)
summary(pol.mod.degree)

pol.mod.occ <- lmer(cv ~ occ.date + (1|Site) +
                        (1|GenusSpecies),
                    data=pol.cv$lm.data)
summary(pol.mod.occ)

## ************************************************************
## the reviewers wanted a quantile regression, can only include on
## random effect at a time

pol.lm.cv.quant <- lqmm(fixed=cv ~ occ.date + degree, random=~1,
                        group= GenusSpecies,
                        data=pol.cv$lm.data,
                        control=list(LP_max_iter=10^4))
pol.sum.boot.quant <- summary.boot.lqmm(boot(pol.lm.cv.quant,
                                             R=100))

## ************************************************************
## plant
## ************************************************************
plant.cv <- cv.trait(spec,
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

plant.mod <- lmer(formula.plant.cv, data=plant.cv$lm.data)
summary(plant.mod)
vif.mer(plant.mod)
## variance inflation factors > 2, so not great (Zurr et al. 2010)

plant.mod.degree <- lmer(cv ~ plant.degree + (1|Site) +
                             (1|GenusSpecies),
                         data=plant.cv$lm.data)

plant.mod.occ <- lmer(cv ~ occ.plant.date + (1|Site) +
                          (1|GenusSpecies),
                      data=plant.cv$lm.data)


## ************************************************************
## the reviewers wanted a quantile regression, can only include on
## random effect at a time

plant.lm.cv.quant <- lqmm(fixed=cv ~ occ.plant.date + plant.degree,
                          random=~1,
                          group= GenusSpecies,
                          data=plant.cv$lm.data,
                          control=list(LP_max_iter=10^3))
plant.sum.boot.quant <- summary.boot.lqmm(boot(plant.lm.cv.quant,
                                               R=100))

## ************************************************************
## correlation betwen degree and persistence

cor.test(pol.cv$lm.data$degree, pol.cv$lm.data$occ.date)
cor.test(plant.cv$lm.data$plant.degree,
         plant.cv$lm.data$occ.plant.date)

## ************************************************************
## residual plots

resid.plot <- function(){
    layout(matrix(1:6, nrow=2, byrow=TRUE))
    plot(fitted(pol.mod), residuals(pol.mod),
         xlab = "Fitted Values", ylab = "Residuals",
         main="Pollinators")
    abline(h=0, lty=2)
    lines(smooth.spline(fitted(pol.mod),
                        residuals(pol.mod)))

    qqnorm(residuals(pol.mod), abline(0,1))

    plot(density(residuals(pol.mod)))

    plot(fitted(plant.mod),
         residuals(plant.mod),
         xlab = "Fitted Values", ylab = "Residuals",
         main="Plants")
    abline(h=0, lty=2)
    lines(smooth.spline(fitted(plant.mod),
                        residuals(plant.mod)))
    qqnorm(residuals(plant.mod), abline(0,1))

    plot(density(residuals(plant.mod)))
}


resid.plot()

## save
save(pol.cv, pol.mod, plant.cv, plant.mod,
     file="saved/contMods.Rdata")

