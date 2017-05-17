options(warn=-1)

## add transparency to named colors
add.alpha <- function(col, alpha=0.2){
    apply(sapply(col, col2rgb)/255, 2,
          function(x)
              rgb(x[1], x[2], x[3],
                  alpha=alpha))
}

plotNet <- function(){
    col.white <- add.alpha("white", alpha=0)
    par(mar=c(0,0.5,0,0))
    layout(matrix(1:length(all.poss.yrs), nrow=1))
    for(j in all.poss.yrs){
        if(j %in% yrs){
            g <- net[[as.character(j)]][plant.sums != 0, pol.sums !=0]
            gs <- graph.incidence(g, weighted=TRUE)
            cols <- c(rep("darkolivegreen", length(rownames(g))),
                      rep("gold", length(colnames(g))))
            importance <-  (c(rowSums(g) +0.1, colSums(g) + 0.1)/sum(g))*175
            cols[importance == min(importance)] <- col.white
            V(gs)$color <- cols
            v.labs <- names(importance)
            v.labs[importance == min(importance)] = ""
            V(gs)$size <- importance
            v.boarders <- rep("black", length(importance))
            v.boarders[importance == min(importance)] <- col.white
            E(gs)$width <- (E(gs)$weight/sum(E(gs)$weight))*40
            gs$layout <- layout_in_circle
            plot.igraph(gs, vertex.label="", vertex.frame.color=v.boarders,
                        vertex.label.cex=importance/10, vertex.label=v.labs)
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
        membs <- data.frame(label=x$label, group=x$group)
        membs <- membs[!grepl("D", membs$label),]
        membs <- membs[!grepl("root", membs$label),]
        membs$group <- as.numeric(membs$group)
        return(membs)
    })

    ## node membership history
    out <- list()
    for(k in 1:length(id.memb)){
        if(k == 1) out <- id.memb[[k]]
        else{
            out <- merge(out, id.memb[[k]], by="label", all=TRUE)
            out[, ncol(out)] <- out[, ncol(out)] + k*20
        }
    }

    save(out, file=sprintf('plotting/saved/nodes_%s.Rdata', i))
    ## counts of nodes that moved from one community to another

    links <- list()
    for(z in 2:ncol(out)){
        if(z==ncol(out)) break
        comm <- table(out[, z:(z+1)])
        if(z==2) links <- comm.mat2sample(comm)
        else{
            links <- rbind(links, comm.mat2sample(comm))
        }
    }
    links <- links[!links$value == 0,]
    save(links, file=sprintf('plotting/saved/links_%s.Rdata', i))
    colors_link <- c('green', 'blue')
    colors_link_array <- paste0("[", paste0("'", colors_link,"'",
                                            collapse = ','), "]")

    colors_node <- c('yellow', 'lightblue', 'red', 'black', 'brown')
    colors_node_array <- paste0("[", paste0("'", colors_node,"'",
                                            collapse = ','), "]")

    plot(gvisSankey(links, from="source",
                    to="target", weight="value",
                    options=list(
                        height=250,
                        sankey="{link: {colors: { fill: ['#a6cee3', '#1f78b4',    '#b2df8a',     '#33a02c'] },
    colorMode: 'gradient'},
                        node: { width: 4,
                                colors: { fill:  ['#a6cee3', '#1f78b4',    '#b2df8a',     '#33a02c'] },
                                label: { fontName: 'Times-Roman',
                                         fontSize: 1,
                                         color: '#871b47',
                                         bold: true,
                                         italic: true }}}"
    ))
    )
    for(j in 1:length(this.tree)){
        comm.tree <- this.tree[[j]]
        l <- layout_with_fr(comm.tree)
        ceb <- cebs[[j]]
        colors <-  c('#a6cee3', '#b2df8a', '#fb9a99', '#fdbf6f',
                     '#cab2d6', '#ffff99', '#1f78b4', '#33a02c')

        plot(ceb, comm.tree,
             vertex.label="",
             edge.width=0.4,
             vertex.size=4,
             edge.arrow.mode=0,
             layout=l, palette=colors,
             mark.border = "white",
             mark.col="white")
        ## igraph:::dendPlot(ceb, mode="hclust",
        ##                   tip.label="",
        ##                   layout=l, palette=colors)
    }
}
