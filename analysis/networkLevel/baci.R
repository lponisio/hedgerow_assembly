rm(list=ls())
setwd('~/Dropbox/hedgerow_assembly/analysis/networkLevel')
source('src/initialize.R')
load('../../data/networks/baci_networks_years.Rdata')
N <- 999

## ************************************************************
## calculate metrics and zscores
## ************************************************************
## mets <- lapply(nets, network.metrics, N)

## cor.dats <- prep.dat(mets,  spec)

## cor.dats$tot.rich <- cor.dats$number.of.species.LL +
##     cor.dats$number.of.species.HL

## save(cor.dats, file='saved/corMets.Rdata')

## ************************************************************
## effect of years post restoration
## ************************************************************
load(file='saved/corMets.Rdata')
baci <- cor.dats[!is.na(cor.dats$ypr),]

## ************************************************************
## linear models

ys <- c("zNODF", "zmod.met.D", "zH2", "niche.overlap.LL",
        "niche.overlap.HL", "connectance")

formulas <-lapply(ys, function(x) {
    as.formula(paste(x, "~",
                     paste("scale(ypr)",
                           "(1|Site)",
                           "(1|Year)",
                           sep="+")))
})

mods <- lapply(formulas, function(x){
    lmer(x,
         data=baci)
})

names(mods) <- ys
## results
lapply(mods, summary)

## ************************************************************
## splines

mods.splines <- lapply(ys, calcSpine)
names(mods.splines) <- ys

layout(matrix(1:length(ys), nrow=2))
lapply(mods.splines, plot)

## ************************************************************
## generalized linear models
poi.ys <- c("number.of.species.HL", "number.of.species.LL", "tot.rich")

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

## ************************************************************
## splines

mods.splines.log <- lapply(poi.ys, calcSpine, calc.log=TRUE)
names(mods.splines.log) <- poi.ys

layout(matrix(1:length(poi.ys), nrow=1))
lapply(mods.splines.log, plot)

save(mods, mods.poi, file="saved/mods/baci_mods.Rdata")
