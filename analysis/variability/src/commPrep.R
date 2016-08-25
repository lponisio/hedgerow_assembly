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

makePretty <- function(comms, sites, sepc){
  names(comms) <- sites
  nyears <- sapply(comms, nrow)
  comms <- comms[nyears >= 3]
  years <- lapply(comms, rownames)

  comm.pp <- list(comm=comms,
                    years=years,
                    sites= rep(names(comms),
                      sapply(comms, nrow)))
  statuses <- spec$SiteStatus[match(comm.pp$sites,
                                    spec$Site)]
  comm.pp$status <- split(statuses,
                            comm.pp$sites)[names(comm.pp$comm)]
  return(comm.pp)
}
