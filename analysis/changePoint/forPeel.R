rm(list=ls())
setwd('~/Dropbox/hedgerow_assembly')
load('data/networks/all_networks_years.Rdata')
load('data/networks/allSpecimens.Rdata')
f.path <- "analysis/changePoint/cptPeel/baci"

library(igraph)
library(parallel)
library(ape)

## creates matrix of all combinations of plants and pollinators and
## fills it 
expandNets <- function(sub.mat, all.mat){
  all.mat[match(rownames(sub.mat), rownames(all.mat)),
          match(colnames(sub.mat), colnames(all.mat))] <- sub.mat
  return(all.mat)
}

pols <- unique(spec$GenusSpecies)
plants <- unique(spec$PlantGenusSpecies)

all.pp <- matrix(0, nrow=length(plants),
                 ncol=length(pols))
rownames(all.pp) <- plants
colnames(all.pp) <- pols

nets <- rapply(nets, expandNets, all.mat=all.pp,
               how="replace")

graphs <- lapply(nets, graph.incidence, weighted=TRUE, add.names=NA)
names(graphs) <- gsub("[.]", "_", names(graphs))
save(graphs, nets, file=file.path(f.path, "graphs_num.Rdata"))

sites <- sapply(strsplit(names(graphs), "_"), function(x) x[1])

## for(i in unique(sites)){
##   g <- graphs[sites == i]
##   verts <- sort(unique(unlist(sapply(g, function(x){
##     V(x)[igraph::degree(x)]
##   }))))
##   out <- data.frame(virtual=0:(length(verts) -1), real=verts)
##   ## colnames(out) <- c("virtual", "real")
##   write.table(out,  row.names=FALSE, sep="\t",
##             file=file.path(f.path, sprintf("%s.lut", i)))
## }

for(i in 1:length(graphs)){
  write.graph(graphs[[i]], file=file.path(f.path,
                             sprintf("%s.pairs",
                                     names(graphs)[i])))
}


graphs <- lapply(nets, graph.incidence, weighted=TRUE)
names(graphs) <- gsub("[.]", "_", names(graphs))
save(graphs, nets, file=file.path(f.path, "graphs.Rdata"))

lutfile <- cbind(0:(length(V(graphs[[1]]))-1), 0:(length(V(graphs[[1]]))-1))
                 
colnames(lutfile) <- c("virtual", "real") 
write.table(lutfile,  row.names=FALSE, sep="\t",
            file=file.path(f.path, "graph-names.lut"))


