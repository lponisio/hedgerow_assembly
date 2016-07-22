rm(list=ls())
setwd('~/Dropbox/hedgerow_assembly/analysis/changePoint')
library(ape)
library(igraph)
source('../networkLevel/src/misc.R')
source('plotting/src/communities.R', chdir = TRUE)

fig.path <- 'plotting/figures'
f.path <- "cptPeel/baci"
load(file=file.path(f.path, "graphs.Rdata"))

temp <- list.files(f.path, pattern="*.gml")
tree.graphs <-  lapply(temp, function(x){
  read.graph(file.path(f.path, x), format="gml")
})

#### plotting

sites <- sapply(strsplit(names(graphs), "_"), function(x) x[1])
sites.trees <- sapply(strsplit(temp, "_"), function(x) x[1])
all.poss.yrs <- 2006:2015

for(i in unique(sites.trees)){
  print(i)
  net <- nets[sites == i]
  yrs <- sapply(strsplit(names(net), "[.]"), function(x) x[2])
  names(net) <- yrs
  arr <- simplify2array(net)
  all.years <-  apply(arr, c(1,2), sum)
  plant.sums <- rowSums(all.years)
  pol.sums <- colSums(all.years)
  this.tree <- tree.graphs[sites.trees == i]
  pdf.f(plotNet,
        file=file.path(fig.path, sprintf("%s_networks.pdf", i)),
        width=18, height=4)
  pdf.f(plotDend,
        file=file.path(fig.path, sprintf("%s_communities.pdf", i)),
        width=10, height=4)
}

names.v <- names(V(graphs[[1]]))
num.v <- 0:length(names.v)

getSpecies <- function(out, names.v, num.v){
  out$names <- names.v[match(out$label, num.v)]
  memb <- apply(out[, -c(1, ncol(out))], 2, table)
  core <- sapply(memb, function(x) names(x)[x == max(x)])
  stayed <- list()
  for(i in 1:length(core)){
    stayed[[i]] <- out[, i+1] == core[i]
  }
  ind.stayed <- Reduce("*", stayed)
  stayed.core <- out$names[ind.stayed ==1]
  left.core <- out$names[ind.stayed==0]
  return(list(stayed.core, left.core))
}

l.path <- 'plotting/saved'
site.nodes <- list.files(l.path, "nodes")
site.cores <- list()

for(i in 1:length(site.nodes)){
  load(file.path(l.path, site.nodes[i]))
  site.cores[[i]] <- getSpecies(out, names.v, num.v)
}


## out <- list()
## for(i in unique(sites)){
##   this.graph <- graphs[sites == i]
##   out[[i]] <- do.call(this.graph, +)
## }

## plotHgr <- function(ihgr){
##   vn <- sub("Actor ", "", V(ihrg)$name)
##   colbar <- rainbow(length(optcom))
##   vc <- ifelse(is.na(V(ihrg)$prob), colbar[V(karate)$comm], "darkblue")
##   V(ihrg)$label <- ifelse(is.na(V(ihrg)$prob), vn, round(V(ihrg)$prob, 2))
##   par(mar=c(0,0.5,0,0))
##   plot(ihrg, vertex.size=10, edge.arrow.size=0.2,
##        vertex.shape="none", vertex.label.color=vc)
## }


## ap <- c(rep("P", length(rownames(nets1))),rep("A", length(colnames(nets1))))
## V(g1)$type <- ap
## colrs <- c(rep("tomato", length(rownames(nets1))),
##            rep("gold", length(colnames(nets1))))
## V(g1)$color <- colrs



## tree.hrg <- lapply(hrgs, hrg_tree)

## dendro.hrg <- lapply(hrgs, hrg.dendrogram)

## ihrgs <- lapply(hrgs, as.igraph)

##  dendPlot(hrgs[[1]])







## ## hrgs
## hrgs <- mclapply(graphs, fit_hrg)
## ## hrg.trees <- mapply(function(a, b)
## ##                     consensus_tree(graph = a,
## ##                                    hrg = b,
## ##                                    num.samples = 10000),
## ##                     a =graphs,
## ##                     b = hrgs,
## ##                     SIMPLIFY = FALSE)

## save(hrgs, file=file.path(f.path, "hrgs.Rdata"))

## ## ## community clusters
## ## clust.gs <- mclapply(graphs, function(g){
## ##   optcom <- cluster_optimal(g)
## ##   V(g)$comm <- membership(optcom)
## ##   return(optcom)
## ## })

## clust.gs <- mclapply(graphs,  cluster_optimal)

## save(clust.gs, file=file.path(f.path, "clusts.Rdata"))                 

## ihrgs <- lapply(hrgs, as.igraph)


## plotDend<- function(){
##   par(mar=c(0,0.5,0,0))
##   layout(matrix(1:length(hrg), nrow=1))
##   for(j in 1:length(hrg)){
##     plot_dendrogram(hrg[[j]])
##   }
## }
