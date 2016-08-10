## used wine nestedness calculator and null model

calc.metric <- function(dat.web, H2_integer=FALSE) {
  ## calculates modularity
  calc.mod <- function(dat.web){ 
    ## converts a p-a matrix to a graph for modularity computation 
    mut.adj <- function(x) {
      nr <- dim(x)[1]
      nc <- dim(x)[2]
      to.fill <- matrix(0, ncol=nc + nr, nrow=nc + nr)
      to.fill[1:nr,(nr+1):(nc+nr)] <- x
      adj.mat <- graph.adjacency(to.fill, mode= "upper", weighted=TRUE)
      return(adj.mat)
    }
    graph <- mut.adj(dat.web)
    weights <- as.vector(dat.web)
    weights <- weights[weights != 0]

    ## if matrix is binary, modularity calculate is not affected by
    ## the weights
    greedy <- modularity(graph,
                         membership(fastgreedy.community(graph,
                                                         weights=weights)),
                         weights=weights)
    
    random.walk <-  modularity(graph,
                               membership(walktrap.community(graph,
                                                             weights=weights)),
                               weights=weights)
    dendro <-  modularity(graph,
                          membership(edge.betweenness.community(graph,
                                                                weights=
                                                                weights)),
                          weights=weights)
    return(c(greedy, random.walk, dendro))
  }

  dat.web <- as.matrix(empty(dat.web))
  ## matrix of all the same number
  if(min(dat.web) == max(dat.web)){
    return(c(mets=rep(0, 1),
             mod.met=rep(0,3)))
  } else{
    nodf <- nestednodf(dat.web,
                       weighted=TRUE,
                       wbinary=TRUE)$statistic["NODF"]
    mets <-  c(nodf,
               networklevel(dat.web,
                      index=c("H2", "connectance", "weighted connectance"),
                      H2_integer=H2_integer))
  }
  mod.met <- calc.mod(dat.web)
  return(c(mets, mod.met= mod.met))
  
}

prob.null <- function(M) {
  randiag <- function(vdims) {
    mats <- diag(vdims)
    ## what cell number are the interactions
    posdiag <- 1:vdims + ((1:vdims - 1) * vdims)
    desp <- matrix(rep(1:vdims, vdims), vdims, vdims, 
                   byrow = TRUE) - (1:vdims)
    fdesp <- sample(1:vdims)
    move <- desp[cbind(1:vdims, fdesp)]
    moved <- posdiag + move
    mdesp <- matrix(0, vdims, vdims)
    mdesp[moved] <- 1
    return(mdesp)
  }
  M <- as.matrix(M)
  values <- M[which(M > 0)]
  lvalues <- length(values)
  if (identical(dim(M)[1], dim(M)[2])) {
    vdims <- dim(M)[1]
    vdimb <- dim(M)[1]
  }
  if (!identical(dim(M)[1], dim(M)[2])) {
    dims <- which(dim(M) == min(dim(M)))
    vdims <- dim(M)[dims]
    dimb <- which(dim(M) == max(dim(M)))
    vdimb <- dim(M)[dimb]
  }
  MR <- matrix(0, vdims, vdimb)
  lMR <- vdims * vdimb
  sample1 <- sample(vdimb, vdims)
  diag1 <- randiag(vdims)
  MR[, sample1] <- diag1
  sample2 <- (1:vdimb)[-sample1]
  pos <- sample(vdims, length(sample2), replace = TRUE)
  MR[cbind(pos, sample2)] <- 2
  MRoccupied <- which(MR > 0)
  vleft <- lvalues - vdimb
  if (vleft > 0) 
    MRoccupied <- c(MRoccupied, sample((1:lMR)[-MRoccupied], 
                                       vleft))
  MR[MRoccupied] <- sample(values)
  if (dim(MR)[1] != dim(M)[1]) 
    MR <- t(MR)
  return(MR)
}

## function to simulate 1 null, and calculate statistics on it
null.stat <- function(dat.web) {
  sim.web <- prob.null(dat.web)
  return(calc.metric(sim.web))
}

##  function that computes summary statistics on simulated null matrices
##  (nulls simulated from web N times)
network.metrics <- function (dat.web, N) {
  ## calculate pvalues
  pvals <- function(stats, nnull){
    rowSums(stats >= stats[, rep(1, ncol(stats))])/(nnull + 1)
  }
  ## calculate zvalues two different ways
  zvals <-function(stats){
    z.sd <- (stats[,1] -
             apply(stats, 1, mean, na.rm = TRUE))/
               apply(stats, 1, sd, na.rm = TRUE)
    z.sd[is.infinite(z.sd)] <- NA
    return(z.sd)
  }
  ## check that matrix is proper format (no empty row/col and no NAs)
  if(all(is.na(dat.web) == FALSE)) {
    ## drop empty rows and columns
    dat.web <- as.matrix(empty(dat.web))
    ## check to make sure emptied matrix is large enough
    ## to calculate statistics on
    if(is.matrix(dat.web)){
      if(all(dim(dat.web) >= 2)) {
        ## calculate null metrics
        null.stat <- replicate(N, null.stat(dat.web),
                               simplify=TRUE)
        ## calculate metrics from data
        true.stat <- calc.metric(dat.web)
        out.mets <- cbind(true.stat, null.stat)
        ## compute z scores
        zvalues <- zvals(out.mets)
        ## compute p-values
        pvalues <- pvals(out.mets, N)
        out <- c(true.stat, zvalues, pvalues)
        names(out) <- c("NODF", "H2", "connectance", "weighted connectance",
                        "modularityG", "modularityR","modularityD",
                        "zNODF", "zH2", "zcon", "zWcon",
                        "zmodG", "zmodR", "zmodD",
                        "pNODF", "pH2", "pcon", "pWcon",
                        "pmodG", "pmodR", "pmodD")
        return(out)
      }
    }
  }
  return(rep(NA,5*3))
}

prep.dat <- function(cor.stats, spec.dat){
  dats <- do.call(rbind, cor.stats)
  out <- data.frame(dats)
  out$Site <- sapply(strsplit(names(cor.stats), "\\."),
                 function(x) x[1])
  out$Year <-  sapply(strsplit(names(cor.stats), "\\."),
                 function(x) x[2])
  out$SiteStatus <- spec.dat$SiteStatus[match(paste(out$Site, out$Year),
                                              paste(spec.dat$Site,
                                                    spec.dat$Year))]
  out$ypr <- spec.dat$ypr[match(paste(out$Site, out$Year),
                                paste(spec.dat$Site,
                                      spec.dat$Year))]
  rownames(out) <- NULL
  return(out)
}



## ##  function that computes summary statistics on simulated null matrices
## ##  (nulls simulated from web N times)
## cor.metrics <- function (true.stat, null.stat, N) {
##   ## calculate pvalues
##   pvals <- function(stats, nnull){
##     colSums(stats >= stats[rep(1, nrow(stats)),])/(nnull + 1)
##   }
##   ## calculate zvalues two different ways
##   zvals <-function(stats){
##     z.sd <- (stats[1,] -
##              apply(stats, 2, mean, na.rm = TRUE))/
##                apply(stats, 2, sd, na.rm = TRUE)
##     z.sd[is.infinite(z.sd)] <- NA
##     return(z.sd)
##   }
##   out.mets <- rbind(true.stat, null.stat)
##   ## compute z scores
##   zvalues <- zvals(out.mets)
##   ## compute p-values
##   pvalues <- pvals(out.mets, N)
##   out <- c(true.stat, zvalues, pvalues)
##   names(out) <- c("NODF", "H2",
##                       "modularityG", "modularityR","modularityD",
##                       "zNODF", "zH2", "zmodG", "zmodR", "zmodD",
##                       "pNODF", "pH2", "pmodG", "pmodR", "pmodD")
##   return(out)
## }

