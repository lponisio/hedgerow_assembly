
getSpecies <- function(out, names.v, num.v){
  out$names <- names.v[match(out$label, num.v)]
  memb <- apply(out[, -c(1, ncol(out))], 2, table)
  ## which are thecommunity names of the "cores"
  core <- sapply(memb, function(x) names(x)[x == max(x)])
  stayed <- list()
  ## what are the species that are in the core
  for(i in 1:length(core)){
    stayed[[i]] <- out[, i+1] == core[i]
  }
  ## which species were always in the core
  ind.stayed <- Reduce("*", stayed)
  ## always in the core
  stayed.core <- out$names[ind.stayed == 1]
  ## left the core at least once
  left.core <- out$names[ind.stayed == 0]
  ## for each year, who was in the core
  in.core <- lapply(stayed, function(x) out$names[x])
  ## for each year, who was in the periphery
  in.periphery <- lapply(stayed, function(x) out$names[!x])
  return(list(across.yrs=list(stayed=stayed.core, left=left.core),
              by.yr=list(core=in.core, periphery=in.periphery)))
}

prepDat <- function(sp.core, spec, by.year=FALSE){
  if(by.year){
    num.sp <- rapply(sp.core, length, how="unlist")
  }else{
    num.sp <- sapply(sp.core, length)
  }
  sp.core <- as.data.frame(unlist(sp.core))
  rownames(sp.core) <- NULL
  colnames(sp.core) <- "GenusSpecies"
  sp.core$SiteStat <- rep(names(num.sp), num.sp)
  sp.core$speciesType <- ifelse(sp.core$GenusSpecies %in% 
                                unique(spec$GenusSpecies),
                                "pollinator", "plant")
  sp.core$Site <- sapply(strsplit(sp.core$SiteStat, "[.]"), 
                         function(x) x[1])
  sp.core$comm <- sapply(strsplit(sp.core$SiteStat, "[.]"),
                         function(x) x[2])
  if(by.year){
    sp.core$ChangePoint <- gsub("[a-z]", "", sp.core$comm)
    sp.core$comm <- sub("[0-9]", "", sp.core$comm)
  }
  sp.core$Abund <- 1
  ## split by species type
  sp.core <- split(sp.core, sp.core$speciesType)
  return(sp.core)
}

makeCores <- function(site.nodes, l.path, names.v, num.v){
  site.cores <- list()
  for(i in 1:length(site.nodes)){
    load(file.path(l.path, site.nodes[i]))
    site.cores[[i]] <- getSpecies(out, names.v, num.v)
  }
  ## create data structure
  names(site.cores) <- unique(sites.trees)
  site.cores.all <- unlist(lapply(site.cores, function(x) x$across.yrs),
                           recursive=FALSE)
  site.cores.all <- prepDat(site.cores.all, spec)
  site.cores.yr <- unlist(lapply(site.cores, function(x) x$by.yr),
                          recursive=FALSE)
  site.cores.yr <- prepDat(site.cores.yr, spec, by.year=TRUE)
  return(list(all=site.cores.all, by.year=site.cores.yr))
}

## create community matrics

makeComm <- function(site.core, by.year=FALSE){
  comm.mat <-  samp2site.spp(site.core$SiteStat,
                             site.core$GenusSpecies,
                             site.core$Abund)
  site <- sapply(strsplit(rownames(comm.mat), "[.]"),
                   function(x) x[1])
  status <- sapply(strsplit(rownames(comm.mat), "[.]"),
                   function(x) x[2])
  if(by.year) status <- sub("[0-9]", "", status)
  return(list(comm=comm.mat, status=status, site=site))
}
