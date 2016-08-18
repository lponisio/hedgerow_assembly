## phylo-beta of interactions
## takes raw specimen data (spec.dat), the , 
calcCommDis <- function(spec.dat, type, lc, abund.w=TRUE){
  ## prep community of sites and interactions
  prep.comm <- aggregate(spec.dat[, type],
                         list(site=spec.dat$Site,
                              year=spec.dat$Year,
                              status=spec.dat$SiteStatus,
                              sp=spec.dat[, type]),
                         length)
  ## split by status then by site
  bystatus <- split(prep.comm, prep.comm$status)
  bysite1 <- lapply(bystatus, function(x){split(x, x$site)})
  bysite <- unlist(bysite1, recursive=FALSE)
  ## create a site by interaction community
  comm <- lapply(bysite, function(y) {
    samp2site.spp(y$year, y$sp, y$x)
  })
  ## drop sites that were only samples one year
  comm <- comm[lapply(comm,nrow) > 1]
  ## calculate the number of years between sampling times
  years <- lapply(comm, rownames)
  dist.time <- as.matrix(dist(unique(spec.dat$Year), diag=TRUE))
  rownames(dist.time) <- colnames(dist.time) <- unique(spec.dat$Year) 
  dists.yrs <- lapply(years, function(x){dist.time[x,x]})
  c.dist <- lapply(dists.yrs, function(x)
                   {as.matrix(x)[lower.tri(x)]})
  ## rename the species as node numbers in preparation for using
  ## comdist using the dendrogram as distances
  
  node.num <- cbind(lc$hclust$labels, paste(lc$edgelist[,1],
                                            lc$edgelist[,2]))
  node.order <- lapply(comm, function(x){
    node.num[,1][match(colnames(x), node.num[,2])]
  })
  assignCols <- function(a,b){
    colnames(a) <- b
    return(a)
  }
  comm.int <- mapply(function(a, b)
                     assignCols(a,b),
                     a=comm,
                     b=node.order)
  dist.tree <- cophenetic(lc$hclust)
  ## calculate the dissimilarity between years based on the dengrogram
  comm.dis <- lapply(comm.int, comdist,
                     dis=dist.tree,
                     abundance.weighted= abund.w)

  ## return a dataframe 
  c.comm <- lapply(comm.dis, function(x)
                   {as.matrix(x)[lower.tri(x)]})
  site <- sapply(strsplit(names(comm.dis), "\\."),
                 function(x) x[2])
  sites <- rep(site, sapply(c.comm, length))
  out <- as.data.frame(cbind(PhyloInt=unlist(c.comm),
                             Dist= unlist(c.dist)))
  out$SiteStatus <- as.factor(sapply(strsplit(rownames(out), "\\."),
                                 function(x) x[1]))
  out$SiteStatus <- factor(out$SiteStatus,
                           levels=c("control", "maturing", "mature")) 
  out$Site <- as.factor(sites)
  rownames(out) <- NULL 
  return(list(phylo.int = out,
              dist.tree=dist.tree,
              comm=comm))
}


calcDendDis <- function(spec.dat, type){
  prep.comm <- aggregate(spec.dat[, type[1]],
                         list(sp=spec.dat[, type[1]],
                              sp2= spec.dat[, type[2]]), length)
  comm <- samp2site.spp(prep.comm$sp, prep.comm$sp2, prep.comm$x)
  comm.dis <- vegdist(comm, "chao", diag= TRUE)
  dengram <- hclust(comm.dis, method="average")

  fden <- function(){
    plot(dengram, cex.lab=0.5, hang=-1)
  }
  path.dir <- 'figures/dgram'
  pdf.f(fden, file= file.path(path.dir, sprintf("%s.pdf",
                paste("dedrogram", type[1], sep="_"))),
        width=25, height=10)
  return(dengram)
}


plot.box <- function(ylabel,
                     dats,
                     y1){
  plot.ci <- function(){
    col.lines <-  brewer.pal(4, "Greys")[c(2,3,4)]
    col.fill <- add.alpha(col.lines, alpha=0.2)
    layout(matrix(1, ncol=1))
    par(oma=c(1, 5, 2, 1),
        mar=c(2, 0, 0.5, 1))
    boxplot(dats[,y1]~ dats$SiteStatus,
            col=col.lines,
             axes = F)
    axis(1, at= 1:3, labels= c("", "", ""), lty= 1, las= .8)
    mtext(c("Non-assembling \n field margin",
              "Assembling \n  hedgerow",
                     "Non-assembling \n hedgerow"),
          side = 1, line= 2, at = 1:3)
    axis(2, pretty(dats[,y1]), labels=pretty(dats[,y1]))
    mtext(ylabel, 2, line=3.5, cex=1.5)
  }
  path <- 'figures' 
  pdf.f(plot.ci, file=file.path(path,
                   sprintf("%s.pdf", paste(
                     gsub(" ", "", ylabel),
                     "box", sep="_"))),
        width=6, height=4)

}
