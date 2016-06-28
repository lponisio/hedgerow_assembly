
## calculates the robustness of a network using Memmot et al.'s method
## takes the adjacency matric, whether to drop species by abundance or
## degree and whther to drop the "higer" or "lower" level of the
## network
## returns a data frame with the site, robustness score and YPR

simExtinction <- function(nets,
                          extinction.method,
                          spec,
                          participant="lower"){
  ext <- lapply(nets, second.extinct,
                participant="lower",
                method=extinction.method)

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
  return(dats)
}                        
