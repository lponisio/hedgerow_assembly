#creating the attributes for plotting dyn nw and defining the colours
nw.for.plot<-function(graph, cola, colp){
  V(graph)$type[V(graph)$type==FALSE]<-"P"
  V(graph)$type[V(graph)$type==TRUE]<-"A"
  V(graph)$color[V(graph)$type=="A"]<-cola
  V(graph)$color[V(graph)$type=="P"]<-colp
  E(graph)$time<-1:length(E(graph))
  return(graph)
}

#deleting the species that do not interact
delete.isolates <- function(graph) {
  g2<-igraph::delete.edges(graph, which(E(graph)$weight <1))
  g2<-igraph::delete.vertices(graph, which(igraph::degree(graph)<1))
  return(g2)
}

save.movie<- function(nets.dyn.list, name, file.path) {
  render.d3movie(nets.dyn.list,
                 filename = file.path(sprintf("%s_video.html", name)),
                 render.par=list(tween.frames=30),
                 label.cex=0.8,label.col='gray',
                 vertex.cex=function(slice){(sna::degree(slice)/3)},
                 # color 
                 vertex.col=ifelse(nets.dyn.list%v%'ap'=='P','darkolivegreen','#E89E0C'),
                 edge.col="black",
                 output.mode = 'HTML'
  )
}