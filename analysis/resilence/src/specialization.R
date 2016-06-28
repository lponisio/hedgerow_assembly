## calculates the species roles from a network and returns a dataframe
## with site status and ypr
## takes networks and specimen data  

calcSpec <- function(nets, spec){
  ## applies specieslevel from bipartite to networks
  species.lev <- lapply(nets, function(x){
    sl <- specieslevel(x)
    sl$'higher level'$tot.int <- colSums(x)
    sl$'lower level'$tot.int <- rowSums(x)
    return(sl)
  })

  ## extract the values and make a dataframe
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

  specs$SiteStatus <- spec$SiteStatus[match(paste(specs$Site,
                                                  specs$assem),
                                            paste(spec$Site,
                                                  spec$Year))]
  specs$SiteStatus <- factor(specs$SiteStatus,
                             levels= c("control", "maturing", "mature"))

  
  return(specs)
}
