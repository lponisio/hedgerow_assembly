rm(list=ls())
setwd('~/Dropbox/hedgerow_assembly/analysis/variability')
binary <- TRUE
alpha <- TRUE
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
                  d= comm.pols$comm,
                  SIMPLIFY=FALSE)

dats.pols <- data.frame(site=comm.pols$sites,
                       status=unlist(comm.pols$status),
                       year=unlist(comm.pols$years),
                       dist=unlist(sapply(dis.pols, function(x)
                         x$distances)))

# invid nulls all not sig, alpha nulls mature marginally sig less,
# occurrence nulls all not sig

## run model, plot
mod.pols <- lmer(dist ~ status +  (1|site),
                        data=dats.pols)
summary(mod.pols)

plot.beta.div(dis.method =dis.method, dists= dats.pols$dist,
              status= dats.pols$status, type= "time",
              sub= "pols", occ=occ)

plot.coeffs(dis.method =dis.method, mod=summary(mod.pols),
              status= dats.pols$status, type= "space",
              sub= "pols", occ=occ)

