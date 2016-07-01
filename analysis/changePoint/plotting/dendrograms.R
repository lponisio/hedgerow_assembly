rm(list=ls())
setwd('~/Dropbox/hedgerow_assembly/analysis/changePoint')
library(ape)
library(igraph)
library(parallel)
source('../networkLevel/src/misc.R')
options(mc.cores=5)
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

plotNet <- function(){
  par(mar=c(0,0.5,0,0))
  layout(matrix(1:length(net), nrow=1))
  for(j in 1:length(net)){
    g <- net[[j]][plant.sums != 0, pol.sums !=0]
    gs <- graph.incidence(g, weighted=TRUE)
    cols <- c(rep("darkolivegreen", length(rownames(g))),
           rep("gold", length(colnames(g))))
    V(gs)$color <- cols
    importance <-  (c(rowSums(g) +0.1, colSums(g) + 0.1)/sum(g))*20
    v.labs <- names(importance)
    v.labs[importance < 5] = ""
    V(gs)$size <- importance
 
    E(gs)$width <- (E(gs)$weight/sum(E(gs)$weight))*20
    gs$layout <- layout_in_circle
    plot.igraph(gs, vertex.label="")
                ## ## vertex.label.cex=importance/10,
                ## vertex.label=v.labs)
  }
}


plotDend <- function(){
  par(mar=c(0,0.5,0,0))
  layout(matrix(1:length(this.tree), nrow=1))
  for(j in 1:length(this.tree)){
    comm.tree <- this.tree[[j]]
    l <- layout_with_fr(comm.tree)
    ceb <- cluster_edge_betweenness(comm.tree)
    plot(ceb, comm.tree, vertex.label="",
         edge.width=0.4,
         vertex.size=4,
         edge.arrow.mode=0,
         layout=l) 
  }
}


for(i in unique(sites)){
  net <- nets[sites == i]
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
  
  ## hrg <- hrgs[sites == i]
  ## pdf.f(plotDend,
  ##       file=file.path(fig.path, sprintf("%s_dendrograms.pdf", i)),
  ##       width=16, height=6)

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
