## rm(list=ls())
## setwd('~/Dropbox/hedgerow_assembly/analysis/networkLevel')
setwd('networkLevel')
source('src/initialize.R')
load('../../data/networks/baci_networks_years.Rdata')

## number of null communities
N <- 999

## ************************************************************
## calculate metrics and zscores ## beware this takes a while!
## ************************************************************
## mets <- lapply(nets, calcNetworkMetrics, N)
## cor.dats <- prepDat(mets,  spec)
## save(cor.dats, file='saved/corMets.Rdata')

## ************************************************************
## effect of years post restoration
## ************************************************************

load(file='saved/corMets.Rdata')
baci <- cor.dats[!is.na(cor.dats$ypr),]
baci$cYear <- as.numeric(baci$Year)-min(as.numeric(baci$Year))


## @knitr external_mets_lm
## ************************************************************
## linear mixed models
ys <- c("zNODF", "zmod.met.D", "zH2", "connectance")

formulas <-lapply(ys, function(x) {
    as.formula(paste(x, "~",
                     paste("scale(ypr)",
                           "(1|Site)",
                           "(1|Year)",
                           sep="+")))
})

formulas.spac <-lapply(ys, function(x) {
    as.formula(paste(x, "~",
                     paste("scale(ypr)",
                           sep="+")))
})

mods <- lapply(formulas, function(x){
    lmer(x,
         data=baci)
})

names(mods) <- ys
## results
lapply(mods, summary)
## @knitr external_mets_lm_end
aic.mods <- sapply(mods, AIC)

## ************************************************************
## generalized linear mixed models
poi.ys <- c("number.of.species.HL", "number.of.species.LL")

formulas.poi <-lapply(poi.ys, function(x) {
    as.formula(paste(x, "~",
                     paste("scale(ypr)",
                           "(1|Site)",
                           "(1|Year)",
                           sep="+")))
})

mods.poi <- lapply(formulas.poi, function(x){
    glmer(x,
          data=baci, family="poisson")
})
names(mods.poi) <- poi.ys

## results
lapply(mods.poi, summary)
## @knitr external_mets_lm_end


## @knitr external_mets_lm_autocorr
## ************************************************************
## standard discrete-time autocorrelation first order model
## takes a while to optimize

formulas.spac <-lapply(ys, function(x) {
    as.formula(paste(x, "~",
                     paste("scale(ypr)",
                           sep="+")))
})

mods.spatial <- lapply(formulas.spac, function(x){
    lme(x,
        random = ~ 1 + cYear | Site,
        correlation=corAR1(form=~cYear),
        data=baci,
        method="REML",
        control=list(maxIter=10^6, niterEM=10^6))
})

names(mods.spatial) <- ys
## results
lapply(mods.spatial, summary)
aic.sp.mods <- sapply(mods.spatial, AIC)

aic.mods-aic.sp.mods
## spatial autocorrelation models have a larger AIC value
## @knitr external_mets_lm_autocorr



## @knitr external_mets_lm_spline
## ************************************************************
## splines

## for gaussian models
mods.splines <- lapply(ys, calcSpine)
names(mods.splines) <- ys

## for poisson models
mods.splines.log <- lapply(poi.ys, calcSpine, calc.log=TRUE)
names(mods.splines.log) <- poi.ys

## @knitr external_mets_lm_spline_end


layout(matrix(1:length(ys), nrow=2))
lapply(mods.splines, plot)


layout(matrix(1:length(poi.ys), nrow=1))
lapply(mods.splines.log, plot)


save(mods, mods.poi, file="saved/mods/baci_mods.Rdata")
