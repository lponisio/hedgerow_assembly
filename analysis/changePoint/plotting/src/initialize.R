library(ape)
library(vegan)
library(FD)
library(igraph)
library(lme4)
library(lmerTest)
source('../networkLevel/src/misc.R')
source('plotting/src/plotPcoa.R')
load('../../data/networks/allSpecimens.Rdata')
source('plotting/src/communities.R')
source('src/calcCore.R')

fig.path <- 'plotting/figures'
f.path <- "cptPeel/baci"
load(file=file.path(f.path, "graphs.Rdata"))

temp <- list.files(f.path, pattern="*.gml")
tree.graphs <-  lapply(temp, function(x){
  read.graph(file.path(f.path, x), format="gml")
})

l.path <- 'plotting/saved'
site.nodes <- list.files(l.path, "nodes")

sites <- sapply(strsplit(names(graphs), "_"), function(x) x[1])
sites.trees <- sapply(strsplit(temp, "_"), function(x) x[1])

names.v <- names(V(graphs[[1]]))
num.v <- 0:length(names.v)
