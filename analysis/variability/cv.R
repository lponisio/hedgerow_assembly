## setwd('~/Dropbox/hedgerow_assembly/analysis/variability')
setwd('variability')

source('src/initialize.R')

## @knitr external_cv
## ************************************************************
## pollinators
## ************************************************************

pol.cv <- calcCvTrait(spec, ## specimen data
                   specs[specs$speciesType =="pollinator",], ## network metrics
                   trait1="occ.date", ## first trait of interest, persistence
                   trait2="r.degree", ## second trait of interest
                   ## rarefied degree
                   time.col="assem", ## name of the time column
                   abund.col="weighted.closeness", ## network metric
                   ## of interest
                   cv.function=cv, ## function to use for calculating cv
                   zero2na=TRUE, ## don't convert NAs to zeros in cv calc
                   standard.cv=TRUE, ## log the cv
                   na.rm=TRUE,
                   species.type="GenusSpecies")

## @knitr external_cv_lm

pol.mod <- lmer(formula.cv, data=pol.cv$lm.data)
print(summary(pol.mod))
vif.mer(pol.mod)

## variance inflation factors < 2, so okay!!! (Zurr et al. 2010)

## @knitr external_cv_quantreg

## ************************************************************
## the reviewers wanted a quantile regression, can only include on
## random effect at a time

## occ.date is significant
pol.lm.cv.quant <- lqmm(fixed=cv ~ occ.date + r.degree, random=~1,
                        group= GenusSpecies,
                        data=pol.cv$lm.data,
                        control=list(LP_max_iter=10^4))
pol.sum.boot.quant <- summary.boot.lqmm(boot(pol.lm.cv.quant,
                                             R=100))

## @knitr external_cv_plants
## ************************************************************
## plants
## ************************************************************

plant.cv <- calcCvTrait(spec,
                     specs[specs$speciesType =="plant",],
                     trait1="occ.plant.date",
                     trait2="plant.r.degree",
                     time.col="assem",
                     abund.col="weighted.closeness",
                     cv.function=corCv,
                     zero2na=TRUE,
                     standard.cv=TRUE,
                     na.rm=TRUE,
                     species.type="PlantGenusSpecies")

plant.mod <- lmer(formula.plant.cv, data=plant.cv$lm.data)
print(summary(plant.mod))
vif.mer(plant.mod)

## variance inflation factors <1, so okay (Zurr et al. 2010)

## ************************************************************
## the reviewers wanted a quantile regression, can only include on
## random effect at a time

## @knitr external_cv_plants_quantreg

plant.lm.cv.quant <- lqmm(fixed=cv ~ occ.plant.date + plant.r.degree,
                          random=~1,
                          group= GenusSpecies,
                          data=plant.cv$lm.data,
                          control=list(LP_max_iter=10^3))
plant.sum.boot.quant <- summary.boot.lqmm(boot(plant.lm.cv.quant,
                                               R=100))

## ************************************************************
## correlation betwen degree and persistence
## @knitr external_cv_cor

cor.test(pol.cv$lm.data$r.degree, pol.cv$lm.data$occ.date)
cor.test(plant.cv$lm.data$plant.r.degree,
         plant.cv$lm.data$occ.plant.date)

## ************************************************************
## @knitr external_cv_resids
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

## passable
resid.plot()

## save
save(pol.cv, pol.mod, plant.cv, plant.mod,
     file="saved/contMods.Rdata")

