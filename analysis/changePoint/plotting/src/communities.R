require(googleVis)

plotNet <- function(){
  par(mar=c(0,0.5,0,0))
  layout(matrix(1:length(all.poss.yrs), nrow=1))
  for(j in all.poss.yrs){
    if(j %in% yrs){
      g <- net[[as.character(j)]][plant.sums != 0, pol.sums !=0]
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
    } else{
      plot(NA,
           ylim=c(0,1),xlim=c(0,1),
           xaxt="n", yaxt="n", ylab="", xlab="", bty="n")
    }
  }
}


comm.mat2sample <-  function (z) {
  temp <- data.frame(expand.grid(dimnames(z))[1:2],
                     as.vector(as.matrix(z)))
  temp <- temp[sort.list(temp[, 1]), ]
  data.frame(source = temp[, 1], target = temp[, 2], value = temp[, 3])
}



plotDend <- function(){
  par(mar=c(0,0.5,0,0))
  layout(matrix(1:length(this.tree), nrow=1))
  cebs <- lapply(this.tree, cluster_edge_betweenness)
  attribs <- lapply(this.tree, get.vertex.attribute)
  id.memb <- lapply(attribs, function(x){
    out <- cbind(x$id, x$group)
    colnames(out) <- c("id", "group")
    return(out)
  })
  ## node membership history
  out <- list()
  for(i in 1:length(id.memb)){
    if(i == 1) out <- id.memb[[i]]
    else{
      out <- merge(out, id.memb[[i]], by="id", all=TRUE)
      out[, ncol(out)] <- out[, ncol(out)] + i*20
    }    
  }
  ## counts of nodes that moved from one community to another
  links <- list()
  for(i in 2:ncol(out)){
    if(i==ncol(out)) break
    comm <- table(out[, i:(i+1)])
    if(i==2) links <- comm.mat2sample(comm)
    else{
      links <- rbind(links, comm.mat2sample(comm))
    }
  }
  links <- links[!links$value == 0,]
  ## nodes <-
  ## data.frame(name=paste("n", 1:length(unique(c(links$target,
  ##              links$source))), sep=""))
  ## nodes$name <- as.character(nodes$name)

  ## sankeyNetwork(Links = links,
  ##               Nodes = nodes,
  ##               Source = 'source', Target = 'target',
  ##               Value = 'value', NodeID = 'name',
  ##               fontSize = 12,  nodeWidth = 30, units = 'TWh')

  plot(gvisSankey(links, from="source", 
                  to="target", weight="value",
                  options=list(
                    height=250,
                    sankey="{link:{color:{fill:'lightblue'}}}"
                    ))
       )
  for(j in 1:length(this.tree)){
    comm.tree <- this.tree[[j]]
    l <- layout_with_fr(comm.tree)
    ceb <- cebs[[j]]

    plot(ceb, comm.tree, vertex.label="",
         edge.width=0.4,
         vertex.size=4,
         edge.arrow.mode=0,
         layout=l) 
  }
}
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
