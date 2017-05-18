## rm(list=ls())
## setwd('~/Dropbox/hedgerow_assembly/analysis/variability')
setwd('variability')
source('src/initialize.R')
source('src/phyloIntBeta.R')
source('src/plotDists.R')

## ************************************************************
## prepare link community in terminal
## ************************************************************
edges.com <- aggregate(list(abund=spec$GenusSpecies),
                       list(GenusSpecies=spec$GenusSpecies,
                            PlantGenusSpecies=spec$PlantGenusSpecies),
                       length)

## lc <- getLinkCommunities(edges.com,
##                          hcmethod = "average",
##                          bipartite=TRUE)
## save(lc, file="saved/lc.Rdata")

## ************************************************************
## turnover of phylo interactions through time
## ************************************************************
load(file="saved/lc.Rdata")
spec$Int <- paste(spec$GenusSpecies,
                  spec$PlantGenusSpecies)

## @knitr external_beta_link
phylo.int <- calcCommDis(spec, "Int", lc, abund.w=TRUE)
## @knitr external_beta_link_end

save(phylo.int, file="saved/phyloInt.Rdata")

## linear model of phylo int by years between samples mature and
## maturing hedgerows have more intearction turnover between years

load(file="saved/phyloInt.Rdata")

## @knitr external_beta_link_lm
phylo.int.mod <- lmer(PhyloInt ~ SiteStatus +
                          (1|Site),
                      data=phylo.int$phylo.int)
print(summary(phylo.int.mod))
## @knitr external_beta_link_end

## ## for comparing maturing and mature
## ## mature and maturing are not sig different

phylo.int$phylo.int$SiteStatus <-
    factor(phylo.int$phylo.int$SiteStatus,
           levels=c("mature", "control", "maturing"))
phylo.int.mod <- lmer(PhyloInt ~ SiteStatus +
                          (1|Site),
                      data=phylo.int$phylo.int)
print(summary(phylo.int.mod))


plot.node.box(ylabel="Weighted interaction turnover",
              dats=phylo.int$phylo.int,
              y1="PhyloInt")

