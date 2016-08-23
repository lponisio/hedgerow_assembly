## creates a year by species community
## species.type referes to the column in the dataset with the species
## names to be aggregated over

calcYearBeta <- function(x, species.type, spec){
  this.spec <- spec[spec$Site == x,]
  prep.comm <- aggregate(list(Abund=this.spec[, species.type]),
                         list(GenusSpecies=this.spec[, species.type],
                           Year=this.spec$Year),
                         length)
  comm <- samp2site.spp(prep.comm$Year,
                        prep.comm$GenusSpecies,
                        prep.comm$Abund)
  return(comm)
}
