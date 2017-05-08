## rm(list=ls())
## setwd('~/Dropbox/hedgerow_assembly/analysis/variability')
setwd('analysis/variability')
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
phylo.int <- calcCommDis(spec, "Int", lc, abund.w=TRUE)
save(phylo.int, file="saved/phyloInt.Rdata")

## linear model of phylo int by years between samples mature and
## maturing hedgerows have more intearction turnover between years

load(file="saved/phyloInt.Rdata")
phylo.int.mod <- lmer(PhyloInt ~ SiteStatus +
                      (1|Site),
                      data=phylo.int$phylo.int)
print(summary(phylo.int.mod))

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

## ## ************************************************************
## interaction lifespans
## ************************************************************

int.lives <- lapply(phylo.int$comm, function(x) sort(colSums(x)))

plot(NA, ylim=c(0,1), xlim=c(0,100))
lapply(int.lives, function(x) points(density(x), type='l'))


## ## ************************************************************
## dengrogram based on shared plant interactions
## ************************************************************

## pols <- cbind(spec$GenusSpecies, spec$PlantGenusSpecies)
## ## unique rows (interactions)
## pols <- as.data.frame(unique(pols, MARGIN=1))
## colnames(pols) <- c('GenusSpecies', 'PlantGenusSpecies')

## clust <- calcDendDis(pols, c('GenusSpecies', 'PlantGenusSpecies'))


## phylo.int.mod <- lmer(PhyloInt ~ Dist*SiteStatus +
##                       (1|Site),
##                       data=phylo.int$phylo.int)

## dd.phylo <- expand.grid(Dist=seq(
##                           from= min(phylo.int$phylo.int$Dist),
##                           to= max(phylo.int$phylo.int$Dist),
##                           length=10),
##                         SiteStatus=c("control", "maturing", "mature"),
##                         PhyloInt = 0)

## phylo.pi <- predict.int(mod= phylo.int.mod,
##                         dd=dd.phylo,
##                         y="PhyloInt")

## plot.predict(new.dd=phylo.pi,
##                  ylabel="Node turnover",
##                  dats=phylo.int$phylo.int,
##                  y1="PhyloInt",
##                  legend.loc="bottomright")
