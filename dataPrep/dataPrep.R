rm(list=ls())
setwd('~/Dropbox/hedgerow_assembly/dataPrep')
load('~/Dropbox/hedgerow/data_sets/traditional/specimens-complete.RData')
load("/Users/lcp/Dropbox/hedgerow_network/data/veg/veg.Rdata")
source('src/misc.R')
source('src/prepNets.R')
library(bipartite)

trait.dir <- '~/Dropbox/hedgerow/data_sets/traditional/functional_traits'
bee.trait <-
  read.csv(file.path(trait.dir, 'bee.csv'),
           row.names=1)
syr.trait <-
  read.csv(file.path(trait.dir, 'syr.csv'),
           row.names=1)

spec <- dd

## subset to net specimens
spec <- spec[spec$NetPan == 'net',]

## subset to just bees and syrphids
spec <- spec[spec$Family == 'Syrphidae' |
             spec$BeeNonbee == 'bee',]

## create species column
spec$PlantGenusSpecies <-  fix.white.space(paste(spec$PlantGenus,
                                                 spec$PlantSpecies))

## drop pollinators and plants without identifications
spec <-  spec[spec$Species != '',]
spec <-  spec[spec$PlantGenusSpecies != '',]
spec$SiteStatus[spec$SiteStatus == "restored"] <- "maturing"

## create a giant network to calculate specialization
agg.spec <- aggregate(list(abund=spec$GenusSpecies),
                      list(GenusSpecies=spec$GenusSpecies,
                           PlantGenusSpecies=spec$PlantGenusSpecies),
                      length)

nets.all <- samp2site.spp(agg.spec$PlantGenusSpecies,
                      agg.spec$GenusSpecies, agg.spec$abund, FUN=sum)

all.specializations <- specieslevel(nets.all,
                                    index=c("proportional generality",
                                    "d"))

traits <- data.frame(GenusSpecies= unlist(sapply(all.specializations,
                       rownames)),
                     do.call(rbind, all.specializations))
rownames(traits) <- NULL

write.csv(traits, file="../data/traits.csv", row.names=FALSE)


## drop forb and natural sites 
to.drop.status <- c("forb", "natural")
spec <- spec[!spec$SiteStatus %in% to.drop.status,]

## add various traits
## specialization
spec$d <- traits$d[match(spec$GenusSpecies, traits$GenusSpecies)]

## occurence
load('~/Dropbox/hedgerow/data_sets/matrices/net/bee.syr.RData')
occ <- apply(mat, c(3,1), function(x){
  sum(x > 0, na.rm=TRUE)/sum(x >= 0, na.rm=TRUE)
})


findOcc <- function(x){
  out <- try(occ[x["GenusSpecies"], paste(x["Site"], x["SiteStatusBACI"],
                                sep=":")], silent=TRUE)
  if(inherits(out, "try-error")) out <- NA
  return(out)
}

spec$occ.date <- apply(spec, 1, findOcc)

## bee functional traits

spec$Lecty <-
  bee.trait$Lecty[match(spec$GenusSpecies, rownames(bee.trait))]
spec$NestLoc <-
  bee.trait$NestLoc[match(spec$GenusSpecies, rownames(bee.trait))]
spec$Excavate <-
  bee.trait$Excavate[match(spec$GenusSpecies, rownames(bee.trait))]
spec$ITD <-
  bee.trait$MeanITD[match(spec$GenusSpecies, rownames(bee.trait))]
spec$Sociality <-
  bee.trait$Sociality[match(spec$GenusSpecies, rownames(bee.trait))]

## syrphid functional traits

spec$LarvalDiet <- syr.trait$LarvalDiet[match(spec$GenusSpecies,
                                              rownames(syr.trait))]
spec$AdultDiet <- syr.trait$AdultDiet[match(spec$GenusSpecies,
                                            rownames(syr.trait))]
spec$WingLength <- syr.trait$WingLength[match(spec$GenusSpecies,
                                              rownames(syr.trait))]

save(spec, file='../data/networks/allSpecimens.Rdata')

site.years <- aggregate(Year~ Site, data=spec,
                        function(x) length(unique(x)))

sites.to.keep <- site.years$Site[site.years$Year >= 5]

## *******************************************************************
## create networks
## all sites with > 5 years
nets <- breakNet(spec[spec$Site %in% sites.to.keep,], 'Site', 'Year')

## save networks for each site, timeframe
f.path <- '../data/networks'
save(nets, file=file.path(f.path, 'all_networks_years.Rdata'))

## *******************************************************************
## keep only BACI sites
BACI.site <- c('Barger', 'Butler', 'Hrdy', 'MullerB', 'Sperandio')
spec <-  spec[spec$Site %in% BACI.site,]
veg <-  veg[veg$Site %in% BACI.site,]

veg$ypr <- spec$ypr[match(paste(veg$Year, veg$Site),
                          paste(spec$Year, spec$Site))]

save(spec, file='../data/networks/specimens.Rdata')

## create a "year" columns grouping years 1-3 post restoration and 4-9
spec$assem <- 'early'
spec$assem[spec$ypr > 3] <- 'late'

veg$assem <- 'early'
veg$assem[veg$ypr > 3] <- 'late'

veg.sum <- aggregate(list(abundance=veg$NumQuads),
                     list(site=veg$Site,
                          status=veg$assem,
                          species=veg$PlantGenusSpecies),
                     sum)

## *******************************************************************
## create networks

## by early/late assembly 
networks <- breakNet(spec, 'Site', 'assem')

## save networks for each site, timeframe
saveDats(networks, names(networks), f.path)
save(networks, file=file.path(f.path, 'networks_stages.Rdata'))

## for each year
nets <- breakNet(spec, 'Site', 'Year')
save(nets, file=file.path(f.path, 'baci_networks_years.Rdata'))

## *******************************************************************
## specialization of each species at each site

species.lev <- lapply(networks, function(x){
  sl <- specieslevel(x)
  sl$'higher level'$tot.int <- colSums(x)
  sl$'lower level'$tot.int <- rowSums(x)
  return(sl)
})

specializations  <-  mapply(function(a, b)
  getSpec(species.lev = a,
          names.net = b),
  a = species.lev,
  b = names(networks),
  SIMPLIFY = FALSE)

specializations <- do.call(rbind, specializations)
rownames(specializations) <- NULL

f.path <- '../data/degree'
save(specializations, file=file.path(f.path, 'specializations.Rdata'))


## *******************************************************************
## change in visits of by the generalized pollinators

## hist(specializations$proportional.generality[specializations$speciesType ==
##                                              "pollinator"],
##      xlab="Generalization")

## plants
diff.gens.plants <- getVisitChange(0, 0.2, "proportional.generality",
                            "PlantGenusSpecies", "GenusSpecies")

diff.spec.plants <- getVisitChange(0.5, 1, "proportional.generality",
                            "PlantGenusSpecies", "GenusSpecies")

diff.all.plants <- getVisitChange(0, 1, "proportional.generality",
                           "PlantGenusSpecies", "GenusSpecies")

## pollinators
diff.gens.pol <- getVisitChange(0, 0.2, "proportional.generality",
                            "GenusSpecies", "PlantGenusSpecies")

diff.spec.pol <- getVisitChange(0.5, 1, "proportional.generality",
                            "GenusSpecies", "PlantGenusSpecies")

diff.all.pol <- getVisitChange(0, 1, "proportional.generality",
                           "GenusSpecies", "PlantGenusSpecies")

f.path <- '../data/degree'
save(diff.gens.plants, diff.spec.plants, diff.all.plants,
     file=file.path(f.path, 'PlantVisitDiffs.Rdata'))
write.csv(diff.all.plants,
          file=file.path(f.path, 'PlantVisitDiffs.csv'), row.names=FALSE)

save(diff.gens.pol, diff.spec.pol, diff.all.pol,
     file=file.path(f.path, 'PolVisitDiffs.Rdata'))
write.csv(diff.all.pol,
          file=file.path(f.path, 'PolVisitDiffs.csv'), row.names=FALSE)


## *******************************************************************
## species lists for each site

plants <- getSpecies(networks, rownames)
pols <- getSpecies(networks, rownames)

f.path <- '../data/species'
write.csv(plants, file.path(f.path, 'plants.csv'),
          row.names=FALSE)
save(plants, file=file.path(f.path, 'plants.Rdata'))

write.csv(pols, file.path(f.path, 'pollinators.csv'),
          row.names=FALSE)
save(pols, file=file.path(f.path, 'pols.Rdata'))

## *******************************************************************
## species added between early and late stages

plant.diffs <- getColExt(plants)
pol.diffs <- getColExt(pols)

f.path <- '../data/speciesChange'
write.csv(plant.diffs, file.path(f.path, 'plants.csv'),
          row.names=FALSE)

write.csv(pol.diffs, file.path(f.path, 'pollinators.csv'),
          row.names=FALSE)

## *******************************************************************
## total plants at a site
f.path <- '../data/species'
plants <- plants[,-2]
plants <- unique(plants)
write.csv(plants, file.path(f.path, 'plants_all.csv'),
          row.names=FALSE)

## *******************************************************************
## pollinator and plant degrees by years post restoration
spec$all <- 'all'

yr.networks <- breakNet(spec, 'all', 'assem')

d.pol <- lapply(yr.networks, colSums)
d.plant <- lapply(yr.networks, rowSums)

f.path <- '../data/degree'
saveDats(d.pol, paste(names(d.pol), 'pollinators', sep="_"), f.path)
saveDats(d.plant, paste(names(d.plant), 'plants', sep="_"), f.path)

by.year <- data.frame(t(sapply(c(d.plant, d.pol), calcStats)))
by.year$group <- rep(c('plants','pollinators'), each=2)
by.year$assembly <- rep(c('early', 'late'), 2)

write.csv(by.year, file.path(f.path, 'stats_by_yr.csv'))

## *******************************************************************
## pollinator and plant degrees across all years and sites

all.networks <- breakNet(spec, 'all', 'all')

d.pol.all <- lapply(all.networks, colSums)
d.plant.all <- lapply(all.networks, rowSums)

saveDats(d.pol.all, 'across_yrs_pollinators', f.path)
saveDats(d.plant.all, 'across_yrs_plants', f.path)

all.dats <- data.frame(t(sapply(c(d.plant.all, d.pol.all), calcStats)))
all.dats$group <- c('plants','pollinators')

write.csv(all.dats, file.path(f.path, 'stats_across_yrs.csv'))

## *******************************************************************
## characteristics of plant colonists
plant.col <- plant.diffs[plant.diffs$class == "colonist",]

plant.col$earlyAbund <- veg.sum$abundance[veg.sum$status ==
                                          "early"][match(
                                           plant.col$species,
                                           veg.sum$species)]
plant.col$lateAbund <- veg.sum$abundance[veg.sum$status ==
                                         "late"][match(
                                          plant.col$species,
                                          veg.sum$species)]

plant.col$degree <- d.plant$all_late[match(plant.col$species,
                                           names(d.plant$all_late))]

## number of species that interact
con.plant <- lapply(yr.networks, getCon, 1)

plant.col$partners <- con.plant$all_late[match(plant.col$species,
                                               names(con.plant$all_late))]


f.path <- '../data/speciesChange'
write.csv(plant.col, file.path(f.path, 'plants_char.csv'),
          row.names=FALSE)

## *******************************************************************
## stable network structure
networks.by.year <- breakNet(spec, 'Site', 'Year')

plant.species <- sapply(networks.by.year, nrow)
plant.species <- data.frame(richness=plant.species,
                            sites =
                              sapply(strsplit(names(plant.species),
                                              "_"),
                                     function(x) x[1]),
                            years =
                              sapply(strsplit(names(plant.species),
                                              "_"),
                                     function(x) x[2]))


write.csv(plant.species, file.path(f.path, 'plants_change.csv'),
          row.names=FALSE)
