## rm(list=ls())
## setwd('~/Dropbox/hedgerow_assembly/analysis/variability')
setwd('variability')

args <- commandArgs(trailingOnly=TRUE)

binary <- FALSE
alpha <- TRUE
## ints or pols

if(length(args) == 0){
    type <- "pols"
} else{
    type <- args[1]
}
source('src/initialize_beta.R')

## ************************************************************
## beta diversity as variation between years,
## centroid for each site
## ************************************************************

## @knitr external_beta_div
dis <- mapply(function(a, b, c, d)
    calcBetaStatus(comm= a, ## observed communities
                   status= b, ## vector of site types
                   dis.method, ## dissimilarity metric
                   nulls=c, ## null communities
                   occ=binary, ## binary or abundance weighted?
                   years=d, ## calculate beta div within?
                   sub=type,
                   zscore=FALSE), ## use Chase method not zscores
    a=comm$comm,
    b=comm$status,
    c= nulls,
    d= comm$comm,
    SIMPLIFY=FALSE)
## @knitr external_beta_div_end

dats <- data.frame(site=comm$sites,
                   status=unlist(comm$status),
                   year=unlist(comm$years),
                   dist=unlist(sapply(dis, function(x)
                       x$distances)))
rownames(dats) <- NULL


save(dats, file= file.path('saved/speciesTurnover',
                           sprintf('%s.pdf', paste(dis.method, alpha, occ, type,
                                                   sep='_'))))

## ************************************************************
## invid nulls all not sig, alpha nulls mature marginally sig less,
## occurrence nulls all not sig

## run model, plot

load(file= file.path('saved/speciesTurnover', sprintf('%s.pdf',
                                                      paste(dis.method, alpha, occ, type, sep='_'))))

## @knitr external_beta_div_lm
mod <- lmer(dist ~ status +  (1|site),
            data=dats)
print(summary(mod))
## @knitr external_beta_div_lm_end


## plots each type individually.
## plot.beta.div(dis.method =dis.method, dists= dats$dist,
##               status= dats$status, type= "time",
##               sub= type, occ=occ, ylabel=ylabel)

## plot.coeffs(dis.method =dis.method, mod=summary(mod),
##             status= dats$status, type= "space",
##             sub= type, occ=occ)

## plots all of the panels that appear in the ms
if(type == "plants") plotDistPanels()

## ************************************************************
##  temporal autocorrelation

dats$cYear <- as.numeric(dats$year) - min(as.numeric(dats$year))

## @knitr external_beta_div_lm_autocorr
mod.spatial <- lme(dist ~ status,
                   random = ~ 1 + cYear | site,
                   correlation=corAR1(form=~cYear),
                   data=dats,
                   control=list(maxIter=10^5, niterEM=10^5))

summary(mod.spatial)

AIC(mod) -AIC(mod.spatial)

## spatial auto-correlation has a larger AIC
