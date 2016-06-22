rm(list=ls())
setwd('~/Dropbox/hedgerow_assembly/analysis/resilence')
source('src/initialize.R')
library(lmerTest)

source("src/prepNets.R")

## create pp matrix for each site, year
nets <- break.net(spec)

## **********************************************************
## robustness
## **********************************************************
ext <- lapply(nets, second.extinct,
              participant="lower",
              method="abund")

rob <- sapply(ext, robustness)

sites <- sapply(strsplit(names(rob), "[.]"), function(x) x[1])
years <- sapply(strsplit(names(rob), "[.]"), function(x) x[2])

dats <- data.frame(Site= sites,
                   Year=years,
                   Robustness=rob)
rownames(dats) <- NULL

dats$SiteStatus <- spec$SiteStatus[match(paste(dats$Site, dats$Year),
                           paste(spec$Site, spec$Year))]

dats$ypr <- spec$ypr[match(paste(dats$Site, dats$Year),
                           paste(spec$Site, spec$Year))]

## no change in robustness by site status
boxplot(dats$Robustness~dats$SiteStatus)
summary(lmer(Robustness~SiteStatus
             + (1|Site) + (1|Year),
             data=dats))

## no effect of ypr on robustness
boxplot(dats$Robustness~dats$ypr)
summary(lmer(Robustness~ypr
             + (1|Site) + (1|Year),
             data=dats[!is.na(dats$ypr),]))

## **********************************************************
## species importance
## **********************************************************

species.lev <- lapply(nets, function(x){
  sl <- specieslevel(x)
  sl$'higher level'$tot.int <- colSums(x)
  sl$'lower level'$tot.int <- rowSums(x)
  return(sl)
})

specs  <-  mapply(function(a, b)
  getSpec(species.lev = a,
          names.net = b,
          seps="[.]"),
  a = species.lev,
  b = names(nets),
  SIMPLIFY = FALSE)

specs <- do.call(rbind, specs)
rownames(specs) <- NULL

specs$ypr <- spec$ypr[match(paste(specs$Site,
                                            specs$assem),
                           paste(spec$Site, spec$Year))]

plot(specs$d[specs$speciesType == "pollinator"] ~
     specs$ypr[specs$speciesType == "pollinator"])

pols <- unique(specs$GenusSpecies[specs$speciesType == "pollinator"])
cols <- rainbow(length(pols))
names(cols) <- pols

spec.metric <- "proportional.generality"

specs$overall.spec <- traits[,spec.metric][match(specs$GenusSpecies,
                                       traits$GenusSpecies)]

specs$specialization <- "spec"
specs$specialization[specs$overall.spec > 0.2] <- "gen"

summary(lmer(proportional.generality ~ ypr
             + (1|Site) + (1|GenusSpecies),
             data= specs[specs$speciesType == "pollinator",]))

summary(lmer(d ~ ypr +
             (1|Site) + (1|GenusSpecies),
             data= specs[specs$speciesType == "plant",]))
