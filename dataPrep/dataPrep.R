## rm(list=ls())
setwd('~/Dropbox/hedgerow_assembly/dataPrep')
load('~/Dropbox/hedgerow/data_sets/traditional/specimens-complete.RData')
load("/Users/lcp/Dropbox/hedgerow_network/data/veg/veg.Rdata")
source('src/misc.R')
source('src/prepNets.R')
source('src/specialization.R')
library(bipartite)
library(fossil)

trait.dir <- '~/Dropbox/hedgerow/data_sets/traditional/functional_traits'
bee.trait <-
  read.csv(file.path(trait.dir, 'bee.csv'),
           row.names=1)
syr.trait <-
  read.csv(file.path(trait.dir, 'syr.csv'),
           row.names=1)

both.trait <-
  read.csv(file.path(trait.dir, 'bee.syr.csv'))

spec <- dd

spec <-  spec[spec$Year != "2015",]

## subset to net specimens
spec <- spec[spec$NetPan == 'net',]
## create species column
spec$PlantGenusSpecies <-  fix.white.space(paste(spec$PlantGenus,
                                                 spec$PlantSpecies))

## drop pollinators and plants without identifications
spec <-  spec[spec$PlantGenusSpecies != '',]

to.drop.status <- c("forb", "natural")
spec <- spec[!spec$SiteStatus %in% to.drop.status,]


tot.spec <- nrow(spec)
print("proportion syrphids")
nrow(spec[spec$Family == 'Syrphidae',])/tot.spec
print("proprotion bee")
nrow(spec[spec$BeeNonbee == 'bee',])/tot.spec

## subset to just bees and syrphids
spec <- spec[spec$Family == 'Syrphidae' |
             spec$BeeNonbee == 'bee',]

spec <-  spec[spec$Species != '',]
spec$SiteStatus[spec$SiteStatus == "restored"] <- "maturing"


## total specimens
print("total specimens")
nrow(spec)

## total species
print("total species")
length(unique(spec$GenusSpecies))

## sampling dates
print("sampling dates")
length(unique(paste(spec$Site, spec$Date)))

## families and genera
print("total families")
length(unique(spec$Family))
print("total genera")
length(unique(spec$Genus))

## interactions
print("total interactions")
length(unique(paste(spec$GenusSpecies, spec$PlantGenusSpecies)))

## *************************************************
## create a giant network to calculate specialization
## *************************************************
agg.spec <- aggregate(list(abund=spec$GenusSpecies),
                      list(GenusSpecies=spec$GenusSpecies,
                           PlantGenusSpecies=spec$PlantGenusSpecies),
                      length)

nets.all <- samp2site.spp(agg.spec$PlantGenusSpecies,
                          agg.spec$GenusSpecies, agg.spec$abund, FUN=sum)

all.specializations <- specieslevel(nets.all,
                                    index=c("proportional generality",
                                      "degree",
                                      "d"))
## calculate rarified plant.pol degree
rare.plants.degree <- apply(nets.all, 1, chao1)
rare.pols.degree <- apply(nets.all, 2, chao1)

traits <- data.frame(GenusSpecies= unlist(sapply(all.specializations,
                       rownames)),
                     do.call(rbind, all.specializations))

traits$r.degree <-  rare.pols.degree[match(traits$GenusSpecies,
                                           names(rare.pols.degree))]

traits$r.degree[is.na(traits$r.degree)] <-
    rare.plants.degree[match(traits$GenusSpecies[is.na(traits$r.degree)],
                             names(rare.plants.degree))]

rownames(traits) <- NULL

## *************************************************
## thermal traits
## *************************************************
spec$AverageTemp <- apply(spec, 1, function(x){
  mean(as.numeric(c(x["TempStart"], x["TempEnd"])),
       na.rm=TRUE)
})

temp.tol <- do.call(rbind, tapply(spec$AverageTemp, spec$GenusSpecies,
                                  function(x){
                                    temp.mean <- mean(x, na.rm=TRUE)
                                    temp.range <- range(x, na.rm=TRUE)
                                    max.temp <- temp.range[2]
                                    temp.range <- temp.range[2] - temp.range[1]
                                    return(c(temp.mean=temp.mean,
                                             max.temp=max.temp,
                                             temp.range=temp.range))
                                  }))
temp.tol <- as.data.frame(temp.tol)
temp.tol$GenusSpecies <- rownames(temp.tol)
rownames(temp.tol) <- NULL

traits <- merge(traits, temp.tol, all.x=TRUE)

## *************************************************
## sampling table for manuscript
## *************************************************

site.table <- aggregate(list(Samples=spec$Date),
                        list(Year=spec$Year, Site=spec$Site),
                        function(x) length(unique(x)))

ms.table <- samp2site.spp(site=site.table$Site,
                          spp=site.table$Year,
                          abund=site.table$Samples,
                          FUN=sum)

write.csv(ms.table, file="../data/samples.csv")

ms.table <- cbind(spec$SiteStatus[match(rownames(ms.table),
                                        spec$Site)], ms.table)

colnames(ms.table) <- c("Site type", colnames(ms.table)[-1])

ms.table <- ms.table[order(ms.table[, "Site type"], decreasing=TRUE),]

ms.table[, "Site type"][ms.table[, "Site type"] == "maturing"] <-
  "Assembling HR"

ms.table[, "Site type"][ms.table[, "Site type"] == "mature"] <-
  "Non-assembling HR"

ms.table[, "Site type"][ms.table[, "Site type"] == "control"] <-
  "Non-assembling FM"

write.table(ms.table, file="~/Dropbox/hedgerow_assembly_ms/ms/tables/samples.txt",
            sep=" & ")

## *************************************************
## add various traits
## *************************************************
## specialization
spec$d <- traits$d[match(spec$GenusSpecies, traits$GenusSpecies)]
spec$degree <- traits$degree[match(spec$GenusSpecies,
                                   traits$GenusSpecies)]
spec$r.degree <- traits$r.degree[match(spec$GenusSpecies,
                                   traits$GenusSpecies)]
spec$plant.degree <- traits$degree[match(spec$PlantGenusSpecies,
                                         traits$GenusSpecies)]
spec$plant.r.degree <- traits$r.degree[match(spec$PlantGenusSpecies,
                                         traits$GenusSpecies)]

## occurence
load('~/Dropbox/hedgerow/data_sets/matrices/net/bee.syr.RData')

occ <- apply(mat, c(3,1), calcOccArray)
spec$occ.date <- apply(spec, 1, findOcc)
traits$occ.date <- spec$occ.date[match(traits$GenusSpecies,
                                       spec$GenusSpecies)]


## plant occurrence
## create sample matrix
site.date <- mat[,,1]
site.date[site.date > 0] <- 0
rownames(site.date) <- lapply(strsplit(rownames(site.date),":"),
                              function(x) x[1])
long.site.date <- comm.mat2sample(site.date)
long.site.date <- long.site.date[!is.na(long.site.date$Samp),]

## create site by date matrices with plant presence
plant.mat <- make.by.species(spec, long.site.date, site.date)
pol.mat <- make.by.species(spec, long.site.date, site.date,
                           type="GenusSpecies")

save(plant.mat, pol.mat, file='../data/species/allSamples.Rdata')

occ.plant <- apply(plant.mat, c(3,1), calcOccArray)
## match to dataset!
spec$occ.plant.date <- apply(spec, 1, findOccPlant)

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

traits <- merge(traits, syr.trait[,c(5:7,10,33)], all.x=TRUE)

traits <- merge(traits, bee.trait[,c(1:5,27)], all.x=TRUE)

traits <- merge(traits, both.trait[,c(2:3,7)], all.x=TRUE)

traits$bee.syr <- spec$BeeNonbee[match(traits$GenusSpecies,
                                       spec$GenusSpecies)]

traits$bee.syr[is.na(traits$bee.syr)] <- "plant"

mean.abund.pol <- apply(apply(pol.mat, c(1,3), mean, na.rm=TRUE), 2,
                        mean)
mean.abund.plant <- apply(apply(plant.mat, c(1,3), mean, na.rm=TRUE),
                          2, mean)
mean.abund <- c(mean.abund.pol, mean.abund.plant)


traits$mean.abun.net <- mean.abund[match(traits$GenusSpecies,
                                         names(mean.abund))]


save(spec, file='../data/networks/allSpecimens.Rdata')
write.csv(traits, file="../data/traits.csv", row.names=FALSE)


site.years <- aggregate(Year~ Site, data=spec,
                        function(x) length(unique(x)))

sites.to.keep <- site.years$Site[site.years$Year >= 5]

## *******************************************************************
## create networks
## all sites with > 5 years
spec.for.nets <- spec[spec$Site %in% sites.to.keep,]

nets <- breakNet(spec.dat=spec.for.nets, 'Site', 'Year')


## save networks for each site, timeframe
f.path <- '../data/networks'
save(nets, file=file.path(f.path, 'all_networks_years.Rdata'))

## *******************************************************************
## site-species level metric calculation
## *******************************************************************
specs <- calcSpec(nets, spec, spec.metric = "d", 0.3)
specs$closeness[specs$closeness == 0] <- 1*10^-6
specs$closeness.log <- log(specs$closeness)
specs.save.path <- '../analysis/speciesLevel/saved'
save(specs, file=file.path(specs.save.path, 'specs.Rdata'))


## *******************************************************************
## keep only BACI sites
## *******************************************************************

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
## *******************************************************************

## ## by early/late assembly
## networks <- breakNet(spec, 'Site', 'assem')

## ## save networks for each site, timeframe
## saveDats(networks, names(networks), f.path)
## save(networks, file=file.path(f.path, 'networks_stages.Rdata'))

## for each year
nets <- breakNet(spec, 'Site', 'Year')
save(nets, file=file.path(f.path, 'baci_networks_years.Rdata'))

## *******************************************************************
## specialization of each species at each site
## *******************************************************************

## species.lev <- lapply(networks, function(x){
##   sl <- specieslevel(x)
##   sl$'higher level'$tot.int <- colSums(x)
##   sl$'lower level'$tot.int <- rowSums(x)
##   return(sl)
## })

## specializations  <-  mapply(function(a, b)
##                             getSpec(species.lev = a,
##                                     names.net = b),
##                             a = species.lev,
##                             b = names(networks),
##                             SIMPLIFY = FALSE)

## specializations <- do.call(rbind, specializations)
## rownames(specializations) <- NULL

## f.path <- '../data/degree'
## save(specializations, file=file.path(f.path, 'specializations.Rdata'))


## ## *******************************************************************
## ## change in visits of by the generalized pollinators

## ## hist(specializations$proportional.generality[specializations$speciesType ==
## ##                                              "pollinator"],
## ##      xlab="Generalization")

## ## plants
## diff.gens.plants <- getVisitChange(0, 0.2, "proportional.generality",
##                                    "PlantGenusSpecies", "GenusSpecies")

## diff.spec.plants <- getVisitChange(0.5, 1, "proportional.generality",
##                                    "PlantGenusSpecies", "GenusSpecies")

## diff.all.plants <- getVisitChange(0, 1, "proportional.generality",
##                                   "PlantGenusSpecies", "GenusSpecies")

## ## pollinators
## diff.gens.pol <- getVisitChange(0, 0.2, "proportional.generality",
##                                 "GenusSpecies", "PlantGenusSpecies")

## diff.spec.pol <- getVisitChange(0.5, 1, "proportional.generality",
##                                 "GenusSpecies", "PlantGenusSpecies")

## diff.all.pol <- getVisitChange(0, 1, "proportional.generality",
##                                "GenusSpecies", "PlantGenusSpecies")

## f.path <- '../data/degree'
## save(diff.gens.plants, diff.spec.plants, diff.all.plants,
##      file=file.path(f.path, 'PlantVisitDiffs.Rdata'))
## write.csv(diff.all.plants,
##           file=file.path(f.path, 'PlantVisitDiffs.csv'), row.names=FALSE)

## save(diff.gens.pol, diff.spec.pol, diff.all.pol,
##      file=file.path(f.path, 'PolVisitDiffs.Rdata'))
## write.csv(diff.all.pol,
##           file=file.path(f.path, 'PolVisitDiffs.csv'), row.names=FALSE)


## ## *******************************************************************
## ## species lists for each site

## plants <- getSpecies(networks, rownames)
## pols <- getSpecies(networks, rownames)

## f.path <- '../data/species'
## write.csv(plants, file.path(f.path, 'plants.csv'),
##           row.names=FALSE)
## save(plants, file=file.path(f.path, 'plants.Rdata'))

## write.csv(pols, file.path(f.path, 'pollinators.csv'),
##           row.names=FALSE)
## save(pols, file=file.path(f.path, 'pols.Rdata'))

## ## *******************************************************************
## ## species added between early and late stages

## plant.diffs <- getColExt(plants)
## pol.diffs <- getColExt(pols)

## f.path <- '../data/speciesChange'
## write.csv(plant.diffs, file.path(f.path, 'plants.csv'),
##           row.names=FALSE)

## write.csv(pol.diffs, file.path(f.path, 'pollinators.csv'),
##           row.names=FALSE)

## ## *******************************************************************
## ## total plants at a site
## f.path <- '../data/species'
## plants <- plants[,-2]
## plants <- unique(plants)
## write.csv(plants, file.path(f.path, 'plants_all.csv'),
##           row.names=FALSE)

## ## *******************************************************************
## ## pollinator and plant degrees by years post restoration
## spec$all <- 'all'

## yr.networks <- breakNet(spec, 'all', 'assem')

## d.pol <- lapply(yr.networks, colSums)
## d.plant <- lapply(yr.networks, rowSums)

## f.path <- '../data/degree'
## saveDats(d.pol, paste(names(d.pol), 'pollinators', sep="_"), f.path)
## saveDats(d.plant, paste(names(d.plant), 'plants', sep="_"), f.path)

## by.year <- data.frame(t(sapply(c(d.plant, d.pol), calcStats)))
## by.year$group <- rep(c('plants','pollinators'), each=2)
## by.year$assembly <- rep(c('early', 'late'), 2)

## write.csv(by.year, file.path(f.path, 'stats_by_yr.csv'))

## ## *******************************************************************
## ## pollinator and plant degrees across all years and sites

## all.networks <- breakNet(spec, 'all', 'all')

## d.pol.all <- lapply(all.networks, colSums)
## d.plant.all <- lapply(all.networks, rowSums)

## saveDats(d.pol.all, 'across_yrs_pollinators', f.path)
## saveDats(d.plant.all, 'across_yrs_plants', f.path)

## all.dats <- data.frame(t(sapply(c(d.plant.all, d.pol.all), calcStats)))
## all.dats$group <- c('plants','pollinators')

## write.csv(all.dats, file.path(f.path, 'stats_across_yrs.csv'))

## ## *******************************************************************
## ## characteristics of plant colonists
## plant.col <- plant.diffs[plant.diffs$class == "colonist",]

## plant.col$earlyAbund <- veg.sum$abundance[veg.sum$status ==
##                                           "early"][match(
##                                             plant.col$species,
##                                             veg.sum$species)]
## plant.col$lateAbund <- veg.sum$abundance[veg.sum$status ==
##                                          "late"][match(
##                                            plant.col$species,
##                                            veg.sum$species)]

## plant.col$degree <- d.plant$all_late[match(plant.col$species,
##                                            names(d.plant$all_late))]

## ## number of species that interact
## con.plant <- lapply(yr.networks, getCon, 1)

## plant.col$partners <- con.plant$all_late[match(plant.col$species,
##                                                names(con.plant$all_late))]


## f.path <- '../data/speciesChange'
## write.csv(plant.col, file.path(f.path, 'plants_char.csv'),
##           row.names=FALSE)

## ## *******************************************************************
## ## stable network structure
## networks.by.year <- breakNet(spec, 'Site', 'Year')

## plant.species <- sapply(networks.by.year, nrow)
## plant.species <- data.frame(richness=plant.species,
##                             sites =
##                             sapply(strsplit(names(plant.species),
##                                             "_"),
##                                    function(x) x[1]),
##                             years =
##                             sapply(strsplit(names(plant.species),
##                                             "_"),
##                                    function(x) x[2]))


## write.csv(plant.species, file.path(f.path, 'plants_change.csv'),
##           row.names=FALSE)
