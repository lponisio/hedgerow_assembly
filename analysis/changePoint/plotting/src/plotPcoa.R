
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
  return(list(stayed=stayed.core, left=left.core))
}


makeCores <- function(site.nodes, l.path, names.v, num.v){
  site.cores <- list()
  for(i in 1:length(site.nodes)){
    load(file.path(l.path, site.nodes[i]))
    site.cores[[i]] <- getSpecies(out, names.v, num.v)
  }

  ## create data structure
  names(site.cores) <- unique(sites.trees)
  site.cores <- unlist(site.cores, recursive=FALSE)
  num.sp <- sapply(site.cores, length)

  site.cores <- as.data.frame(unlist(site.cores))
  rownames(site.cores) <- NULL
  colnames(site.cores) <- "GenusSpecies"
  site.cores$SiteStat <- rep(names(num.sp), num.sp)
  site.cores$speciesType <- ifelse(site.cores$GenusSpecies %in% 
                                   unique(spec$GenusSpecies),
                                   "pollinator", "plant")
  site.cores$Site <- sapply(strsplit(site.cores$SiteStat, "[.]"), 
                            function(x) x[1])
  site.cores$comm <- sapply(strsplit(site.cores$SiteStat, "[.]"),
                            function(x) x[2])

  site.cores$Abund <- 1
  ## split by species type
  site.cores <- split(site.cores, site.cores$speciesType)
  
  return(site.cores)
}

## create community matrics

makeComm <- function(site.core){
  comm.mat <-  samp2site.spp(site.core$SiteStat,
                             site.core$GenusSpecies,
                             site.core$Abund)

  status <- sapply(strsplit(rownames(comm.mat), "[.]"),
                   function(x) x[2])
  return(list(comm=comm.mat, status=status))
}

## site by species communityes, pollinators then plants in a list,
## a list of vectors of statuses, fig.path, dissimilarity method

plot.beta.div <- function(comms, statuses, fig.path, method){
  f <- function(){
    pcoa <- function(comm, status, i, method){
      pcoa.comm <-  cmdscale(vegdist(comm, method=method))
      plot(pcoa.comm[status == 'left',], asp=1,
           col=cols[1], pch=16, cex=1.5,
           ylim=range(pcoa.comm[,2]) + c(0,0.2),
           xlim=range(pcoa.comm[,1]),
           xlab='',
           ylab='',
           xaxt='n',
           yaxt='n',
           cex.lab=1.5)
      points(pcoa.comm[status == 'stayed',],
             col=cols[2], pch=16, cex=1.5)
      ordihull(pcoa.comm, status)
      mtext('PCoA2', 2, line=1.5, cex=1.5)
      if(i==2){
        mtext('Plants', 2, line=4, cex=1.5)
        mtext('PCoA1', 1, line=1.2, cex=1.5)
      }
      
      if(i==1){
        legend('topleft',
               legend=c('Not core', 'Core'),
               col=cols, pch=16, bty='n', cex=1.1)
        mtext('Bees', 2, line=4, cex=1.5)
      }
    }
    layout(matrix(1:2, ncol=1, byrow=TRUE))
    par(oma=c(2,7,1,1), mar=c(0.5,0,1,0.5),
        mgp=c(2,1,0), cex.axis=1.5)
    cols <- c("darkolivegreen3",
              "darkgoldenrod1")
    for(i in 1:length(comms)){
      pcoa(comms[[i]], statuses[[i]], i=i,
           method=method)
    }
  }

  pdf.f(f, file= file.path(fig.path, "pcoa/pcoa.pdf"),
        width=4, height=5)
}
