rm(list=ls())
library(igraph)
library(networkDynamic)
library(ndtv)
source('~/Dropbox/hedgerow_assembly/analysis/changePoint/plotting/videos/movie_functions.R')

load('~/Dropbox/hedgerow_assembly/analysis/changePoint/cptPeel/baci/graphs.Rdata')
file.path <- '~/Dropbox/hedgerow_assembly/analysis/changePoint/plotting/videos/'
sites <- sapply(strsplit(names(graphs), "_"), function(x) x[1])

#this for is not working
for(i in 1:length(unique(sites))){
  netsS<- nets[sites ==unique(sites)[i]]
  name<-unique(sites)[i]
  nets.dyn<-lapply(netsS,as.network.matrix,matrix.type='incidence',
                 bipartite=TRUE, names.eval='sample')
  nets.dyn.list<-networkDynamic(network.list=nets.dyn)
  #creating an attribute to allow for the degree
  nets.dyn.list%n%'slice.par'<-list(start=1,end=10,interval=1,
                                    aggregate.dur=1,rule='latest')
  #creating just one object to extract the species 
  nets1<-netsS[[1]]
  nets1.dyn<-as.network.matrix(nets1,matrix.type='incidence',
                   bipartite=TRUE, names.eval='sample')
  ap<-c(rep("P", dim(nets1)[1]), rep("A", (dim(nets1)[2])))
  nets.dyn.list%v%'ap'<-ap
  nets.dyn.list%v%'species'<-network.vertex.names(nets1.dyn)
  #animating - the tween.frames arg defines the duration of the "between" frames
  render.animation(nets.dyn.list, render.par = list(tween.frames = 25))
  compute.animation(nets.dyn.list,
                  animation.mode='MDSJ',
                  default.dist=1,
                  verbose=FALSE)
  save.movie(nets.dyn.list=nets.dyn.list, name=name, file.path=file.path)
}
