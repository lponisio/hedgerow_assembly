rm(list=ls())
setwd('~/Dropbox/hedgerow_assembly/analysis/variability')
binary <- FALSE
alpha <- TRUE
## int or pols
type <- "int"
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
rownames(dats.pols) <- NULL

save(dats.pols, file= file.path('saved/speciesTurnover',
                  sprintf('%s.pdf', paste(dis.method, alpha, occ, type,
                                          sep='_'))))

## invid nulls all not sig, alpha nulls mature marginally sig less,
## occurrence nulls all not sig

## run model, plot

load(file= file.path('saved/speciesTurnover', sprintf('%s.pdf',
       paste(dis.method, alpha, occ, type, sep='_'))))


mod.pols <- lmer(dist ~ status +  (1|site),
                 data=dats.pols)
summary(mod.pols)

plot.beta.div(dis.method =dis.method, dists= dats.pols$dist,
              status= dats.pols$status, type= "time",
              sub= type, occ=occ, ylabel=ylabel)

plot.coeffs(dis.method =dis.method, mod=summary(mod.pols),
            status= dats.pols$status, type= "space",
            sub= type, occ=occ)


## plotDistPanels()
