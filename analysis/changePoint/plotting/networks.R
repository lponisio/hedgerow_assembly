rm(list=ls())
setwd('analysis/changePoint')
## setwd("~/Dropbox/hedgerow_assembly/analysis/changePoint")
source('plotting/src/initialize.R')

sites <- sapply(strsplit(names(graphs), "_"), function(x) x[1])
sites.trees <- sapply(strsplit(temp, "_"), function(x) x[1])
all.poss.yrs <- 2006:2015

for(i in unique(sites)){
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
    ## pdf.f(plotDend,
    ##       file=file.path(fig.path, sprintf("%s_dend.pdf", i)),
    ##       width=10, height=4)
}
