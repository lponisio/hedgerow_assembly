rm(list=ls())
setwd('~/Dropbox/hedgerow_assembly/analysis/variability')
binary <- FALSE
alpha <- TRUE

library(lme4)
library(vegan)
library(lmerTest)
source('src/misc.R')
source('src/beta.R')
## source('src/plotDists.R')

source('src/initialize_beta.R')

## ************************************************************
## beta diversity as variation between years,
## centroid for each site
## ************************************************************
dis.pols <- mapply(function(a, b, c, d)
                  beta.status(comm= a,
                              status= b,
                              dis.method,
                              nulls=c,
                              occ=binary,
                              years=d,
                              sub="pols"),
                  a=comm.pols$comm,
                  b=comm.pols$status,
                  c= nulls.pols,
                  d= unique(comm.pols$sites),
                  SIMPLIFY=FALSE)

dats.pols <- data.frame(site=unlist(comm.pols$sites),
                       status=unlist(comm.pols$status),
                       year=comm.pols$years,
                       dist=unlist(sapply(dis.pols, function(x)
                         x$distances)))

## run model, plot
mod.pols <- summary(lmer(dist ~ status +  (1|site) + (1|year),
                        data=dats.pols))

plot.beta.div(dis.method =dis.method, dists= dats.pols$dist,
              status= dats.pols$status, type= "space",
              sub= "pols", occ=occ)

plot.coeffs(dis.method =dis.method, mod=mod.pols,
              status= dats.pols$status, type= "space",
              sub= "pols", occ=occ)

dis.pols.pcoa <- mapply(function(a, b)
                 adonis(formula= a ~ b,
                              method=dis.method),
                       a=comm.pols$comm[3:7],
                       b=comm.pols$status[3:7],
                  SIMPLIFY=FALSE)
source('src/plotPcoa.R')

