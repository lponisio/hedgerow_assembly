
## site by species communityes, pollinators then plants in a list,
## a list of vectors of statuses, fig.path, dissimilarity method

plot.beta.div <- function(comms, statuses, fig.path, method,
                          treatments=c("core", "periphery")){
  f <- function(){
    pcoa <- function(comm, status, i, method){
      pcoa.comm <-  cmdscale(vegdist(comm, method=method))
      plot(pcoa.comm[status == treatments[1],], asp=1,
           col=cols[1], pch=16, cex=1.5,
           ylim=range(pcoa.comm[,2]) + c(-0.1,0.1),
           xlim=range(pcoa.comm[,1]),
           xlab='',
           ylab='',
           xaxt='n',
           yaxt='n',
           cex.lab=1.5)
      points(pcoa.comm[status == treatments[2],],
             col=cols[2], pch=16, cex=1.5)
      ordihull(pcoa.comm, status)
      mtext('PCoA2', 2, line=1.5, cex=1.5)
      if(i==2){
        mtext('Plants', 2, line=4, cex=1.5)
        mtext('PCoA1', 1, line=1.2, cex=1.5)
      }
      
      if(i==1){
        legend('topleft',
               legend=c('Core', 'Periphery'),
               col=cols, pch=16, bty='n', cex=1.1)
        mtext('Bees', 2, line=4, cex=1.5)
      }
    }
    layout(matrix(1:2, ncol=1, byrow=TRUE))
    par(oma=c(2,7,1,1), mar=c(0.5,0,1,0.5),
        mgp=c(2,1,0), cex.axis=1.5)
    cols <- c("darksalmon",
              "black")
    for(i in 1:length(comms)){
      pcoa(comms[[i]], statuses[[i]], i=i,
           method=method)
    }
  }

  pdf.f(f, file= file.path(fig.path, "pcoa/pcoa.pdf"),
        width=4, height=5)
}
