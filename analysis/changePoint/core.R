rm(list=ls())
setwd('~/Dropbox/hedgerow_assembly/analysis/changePoint')
source('plotting/src/initialize.R')
method <- "jaccard"
site.cores <- makeCores(site.nodes, l.path, names.v, num.v)

## create community matrics
comm.mats <- lapply(site.cores, makeComm)

## run perm anovas
## plants
plant.perm <- adonis(comm.mats$plant$comm ~ comm.mats$plant$status,
                     method=method)
plant.beta <- betadisper(vegdist(comm.mats$plant$comm,
                                 method=method),
                         group=comm.mats$plant$status)
permutest(plant.beta, pairwise = TRUE, permutations = 99)

## pollinators
pol.perm <- adonis(comm.mats$pollinator$comm ~ comm.mats$pollinator$status,
                   method=method)
pol.beta <- betadisper(vegdist(comm.mats$pollinator$comm,
                               method=method),
                         group=comm.mats$pollinator$status)
permutest(pol.beta, pairwise = TRUE, permutations = 99)

## plotting
plot.beta.div(list(comm.mats$pollinator$comm, comm.mats$plant$comm),
              list(comm.mats$pollinator$status, comm.mats$plant$status),
              fig.path,
              method= method)

## who are the core species? Species that are never in the core?
## pollinators

stayed.pol <-
  site.cores$pollinator$GenusSpecies[site.cores$pollinator$comm ==
                                     "stayed"]
count.stayed.pol <- table(stayed.pol)
count.stayed.pol <- count.stayed.pol[count.stayed.pol != 0]
always.stayed.pol <- count.stayed.pol[count.stayed.pol == 5]

hist(count.stayed.pol, breaks=0:5)

left.pol <-
  site.cores$pollinator$GenusSpecies[site.cores$pollinator$comm ==
                                     "left"]
count.left.pol <- table(left.pol)
count.left.pol <- count.left.pol[count.left.pol != 0]
always.left.pol <- count.left.pol[count.left.pol == 5]

## plants
stayed.plant <-
  site.cores$plant$GenusSpecies[site.cores$plant$comm ==
                                "stayed"]
count.stayed.plant <- table(stayed.plant)
count.stayed.plant <- count.stayed.plant[count.stayed.plant != 0]
always.stayed.plant <- count.stayed.plant[count.stayed.plant == 5]

hist(count.stayed.plant, breaks=0:5)

left.plant <-
  site.cores$plant$GenusSpecies[site.cores$plant$comm ==
                                "left"]
count.left.plant <- table(left.plant)
count.left.plant <- count.left.plant[count.left.plant != 0]
always.left.plant <- count.left.plant[count.left.plant == 5]


## traits
load('~/Dropbox/hedgerow_network/analysis/functional_traits/saved/traitsbee.Rdata')

bee.comm <- comm.mats$pollinator$comm[, colnames(comm.mats$pollinator$comm)
                                          %in% rownames(traits)]
traits <- traits[rownames(traits) %in% colnames(bee.comm),]
traits <- traits[colnames(bee.comm),]

bee.fdiv <- dbFD(traits, bee.comm)

sites.bee <- sapply(strsplit(rownames(bee.comm), "[.]"),
                    function(x) x[1])
statuses.bee <- sapply(strsplit(rownames(bee.comm), "[.]"),
                    function(x) x[2])

fdis.bee.mod <- lmer(bee.fdiv$FDis ~ statuses.bee + (1|sites.bee))
summary(fdis.bee.mod)

fdiv.bee.mod <- lmer(bee.fdiv$FDiv ~ statuses.bee + (1|sites.bee))
summary(fdiv.bee.mod)
