rm(list=ls())
setwd('~/Dropbox/hedgerow_assembly/analysis/changePoint')
library(ape)
library(igraph)
library(parallel)
source('../assembly/src/misc.R', chdir = TRUE)
options(mc.cores=5)
fig.path <- 'plotting/figures'
f.path <- "cptPeel/baci"
load(file=file.path(f.path, "graphs.Rdata"))

## hrgs
## hrgs <- mclapply(graphs, fit_hrg)
## hrg.trees <- mapply(function(a, b)
##                     consensus_tree(graph = a,
##                                    hrg = b,
##                                    num.samples = 10000),
##                     a =graphs,
##                     b = hrgs,
##                     SIMPLIFY = FALSE)

## save(hrgs, hrg.trees, file=file.path(f.path, "hrgs.Rdata"))



## ## community clusters
## clust.gs <- mclapply(graphs, function(g){
##   optcom <- cluster_optimal(g)
##   V(g)$comm <- membership(optcom)
##   return(optcom)
## })

clust.gs <- mclapply(graphs,  cluster_optimal)


save(clust.gs, file=file.path(f.path, "clusts.Rdata"))                 

#### plotting
load(file=file.path(f.path, "hrgs.Rdata"))
load(file=file.path(f.path, "clusts.Rdata"))

ihrgs <- lapply(hrgs, as.igraph)
lapply(ihrgs, function(x) x$layout <- layout.reingold.tilford)

sites <- sapply(strsplit(names(graphs), "_"), function(x) x[1])

plotNet <- function(){
  par(mar=c(0,0.5,0,0))
  layout(matrix(1:length(ihrg), nrow=1))
  for(j in 1:length(ihrg)){
    plot(ihrg[[j]], vertex.size=10,
         edge.arrow.size=0.2, vertex.label="")
  }
}

plotDend<- function(){
  par(mar=c(0,0.5,0,0))
  layout(matrix(1:length(hrg), nrow=1))
  for(j in 1:length(hrg)){
    plot_dendrogram(hrg[[j]])
  }
}

for(i in unique(sites)){
  ihrg <- ihrgs[sites == i]
  hrg <- hrgs[sites == i]
  pdf.f(plotNet,
        file=file.path(fig.path, sprintf("%s_networks.pdf", i)),
        width=16, height=2)
  pdf.f(plotDend,
        file=file.path(fig.path, sprintf("%s_dendrograms.pdf", i)),
        width=16, height=6)

}

plotHgr <- function(ihgr){
  vn <- sub("Actor ", "", V(ihrg)$name)
  colbar <- rainbow(length(optcom))
  vc <- ifelse(is.na(V(ihrg)$prob), colbar[V(karate)$comm], "darkblue")
  V(ihrg)$label <- ifelse(is.na(V(ihrg)$prob), vn, round(V(ihrg)$prob, 2))
  par(mar=c(0,0.5,0,0))
  plot(ihrg, vertex.size=10, edge.arrow.size=0.2,
       vertex.shape="none", vertex.label.color=vc)
}




tree.hrg <- lapply(hrgs, hrg_tree)

dendro.hrg <- lapply(hrgs, hrg.dendrogram)

ihrgs <- lapply(hrgs, as.igraph)

 dendPlot(hrgs[[1]])



