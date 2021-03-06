rm(list=ls())
setwd('changePoint')
source('plotting/src/initialize.R')
method <- "jaccard"

## **********************************************************
## catagorize species as core or periferal
## **********************************************************
site.cores <- makeCores(site.nodes, l.path, names.v, num.v)

## create community matrics
comm.mats.all <- lapply(site.cores$all, makeComm)
comm.mats.yr <- lapply(site.cores$by.year, makeComm, by.year=TRUE)

## **********************************************************
## perm anovas and dispersion of species that always stayed of left at
## least once
## **********************************************************
## plants
plant.perm.all <- adonis(comm.mats.all$plant$comm ~
                     comm.mats.all$plant$status,
                     method=method)
plant.perm.all

plant.beta.all <- betadisper(vegdist(comm.mats.all$plant$comm,
                                 method=method),
                         group=comm.mats.all$plant$status)
permutest(plant.beta.all, pairwise = TRUE,
          permutations = 999)

## pollinators
pol.perm.all <- adonis(comm.mats.all$pollinator$comm ~
                   comm.mats.all$pollinator$status,
                   method=method)
pol.perm.all
pol.beta.all <- betadisper(vegdist(comm.mats.all$pollinator$comm,
                               method=method),
                       group=comm.mats.all$pollinator$status)
permutest(pol.beta.all, pairwise = TRUE, permutations = 999)

## **********************************************************
## perm anovas and dispersion by year
## **********************************************************
## plants
plant.perm <- adonis(comm.mats.yr$plant$comm ~
                     comm.mats.yr$plant$status,
                     method=method,
                     strata=comm.mats.yr$plant$site)
plant.beta <- betadisper(vegdist(comm.mats.yr$plant$comm,
                                 method=method),
                         group=comm.mats.yr$plant$status)
permutest(plant.beta, pairwise = TRUE, permutations = 999,
          strata=comm.mats.yr$plant$site)

## pollinators
pol.perm <- adonis(comm.mats.yr$pollinator$comm ~
                   comm.mats.yr$pollinator$status,
                   method=method,
                   strata=comm.mats.yr$pollinator$site)
pol.beta <- betadisper(vegdist(comm.mats.yr$pollinator$comm,
                               method=method),
                       group=comm.mats.yr$pollinator$status)
permutest(pol.beta, pairwise = TRUE, permutations = 999,
          strata=comm.mats.yr$pollinator$site)

## plotting
plot.beta.div(list(comm.mats.yr$pollinator$comm,
                   comm.mats.yr$plant$comm),
              list(comm.mats.yr$pollinator$status,
                   comm.mats.yr$plant$status),
              fig.path,
              method= method)

## **********************************************************
## who are the core species? Species that are never in the core?
## pollinators
## **********************************************************
stayed.pol <-
  site.cores$all$pollinator$GenusSpecies[site.cores$all$pollinator$comm ==
                                         "stayed"]
count.stayed.pol <- table(stayed.pol)
count.stayed.pol <- count.stayed.pol[count.stayed.pol != 0]
always.stayed.pol <- count.stayed.pol[count.stayed.pol == 5]

left.pol <-
  site.cores$all$pollinator$GenusSpecies[site.cores$all$pollinator$comm ==
                                         "left"]
count.left.pol <- table(left.pol)
count.left.pol <- count.left.pol[count.left.pol != 0]
always.left.pol <- count.left.pol[count.left.pol == 5]

## plants
stayed.plant <-
  site.cores$all$plant$GenusSpecies[site.cores$all$plant$comm ==
                                    "stayed"]
count.stayed.plant <- table(stayed.plant)
count.stayed.plant <- count.stayed.plant[count.stayed.plant != 0]
always.stayed.plant <- count.stayed.plant[count.stayed.plant == 5]

hist(count.stayed.plant, breaks=0:5)

left.plant <-
  site.cores$all$plant$GenusSpecies[site.cores$all$plant$comm ==
                                    "left"]
count.left.plant <- table(left.plant)
count.left.plant <- count.left.plant[count.left.plant != 0]
always.left.plant <- count.left.plant[count.left.plant == 5]

## **********************************************************
## trait diversity between core and peripheral
## **********************************************************
## traits
load.path <- "../data/saved/"
load(file.path(load.path, 'traitsbee.Rdata'))

bee.comm <- comm.mats.yr$pollinator$comm[,
                                         colnames(comm.mats.yr$pollinator$comm)
                                         %in% rownames(traits)]

traits <- traits[rownames(traits) %in% colnames(bee.comm),]
traits <- traits[colnames(bee.comm),]

bee.fdiv <- dbFD(traits, bee.comm)

sites.bee <- sapply(strsplit(rownames(bee.comm), "[.]"),
                    function(x) x[1])
statuses.bee <- sapply(strsplit(rownames(bee.comm), "[.]"),
                       function(x) x[2])
statuses.bee <- sub("[0-9]", "", statuses.bee)

fdis.bee.mod <- lmer(bee.fdiv$FDis ~ statuses.bee + (1|sites.bee))
summary(fdis.bee.mod)

fdiv.bee.mod <- lmer(bee.fdiv$FDiv ~ statuses.bee + (1|sites.bee))
summary(fdiv.bee.mod)

feve.bee.mod <- lmer(bee.fdiv$FEve ~ statuses.bee + (1|sites.bee))
summary(feve.bee.mod)
