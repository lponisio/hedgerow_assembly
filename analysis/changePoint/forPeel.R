rm(list=ls())
setwd('~/Dropbox/hedgerow_assembly')
load('data/networks/networks_years.Rdata')
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

mut.adj <- function(x) {
  nr <- dim(x)[1]
  nc <- dim(x)[2]
  to.fill <- matrix(0, ncol=nc + nr, nrow=nc + nr)
  to.fill[1:nr,(nr+1):(nc+nr)] <- x
  adj.mat <- graph.adjacency(to.fill, mode= "upper", weighted=TRUE)
  return(adj.mat)
}

graphs <- lapply(nets, mut.adj)
names(graphs) <- gsub("[.]", "_", names(graphs))

lutfile <- cbind(0:(length(V(graphs[[1]]))-1), 0:(length(V(graphs[[1]]))-1))
                 
colnames(lutfile) <- c("virtual", "real") 
write.table(lutfile,  row.names=FALSE, sep="\t",
            file=file.path(f.path, "graph-names.lut"))

for(i in 1:length(graphs)){
  write.graph(graphs[[i]], file=file.path(f.path,
                             sprintf("%s.pairs",
                                     names(graphs)[i])))
}

save(graphs, file=file.path(f.path, "graphs.Rdata"))

