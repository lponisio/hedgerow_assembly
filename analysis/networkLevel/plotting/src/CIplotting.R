plot.panel <- function(dats,
                       new.dd,
                       y1,
                       xs,
                       col.lines,
                       col.fill,
                       col.points=col.fill,
                       plotPoints=TRUE,
                       agg.col="Site",...){
  plotting.loop <- function(){
    ys <- aggregate(list(y=dats[,y1]),
                    list(Site=dats[, agg.col],
                         x=dats[,xs]),
                    mean, na.rm=TRUE)
    if(plotPoints){
      points(x=jitter(ys$x, factor=0.25),
             y=ys$y,
             pch=16,
             col=col.points,
             cex=1.2)
    }
    lines(x=new.dd[,xs],
          y=new.dd[,y1],
          col=col.lines,
          lwd=2)
    lines(x=new.dd[,xs],
          y=new.dd$plo,
          col=col.lines,
          lty="dashed")
    lines(x=new.dd[,xs],
          y=new.dd$phi,
          col=col.lines,
          lty="dashed")
    polygon(c(new.dd[,xs],
              rev(new.dd[,xs])),
            c(new.dd$phi,
              rev(new.dd$plo)),
            col=col.fill, border=NA)
  }
  plot(NA,
       xlim=range(new.dd[,xs], na.rm=TRUE),
       ylim=range(c(0, new.dd$phi,  new.dd$plo, quantile(dats[,y1],
         na.rm=TRUE)["75%"],
         na.rm=TRUE)),
       xlab="",
       ylab="",
       xaxt="n",
       las=1,...)
  plotting.loop()
}


plot.predict.ypr <- function(new.dd,
                             ylabel,
                             dats,
                             y1,
                             xs='ypr',
                             plotPoints=TRUE,
                             xlabel= 'Year of assembly',
                             path = 'figures',
                             extinction.method,
                             agg.col="Site"){
  plot.ci <- function(){
    col.lines <-  brewer.pal(4, "Greys")[3]
    col.fill <- add.alpha(col.lines, alpha=0.2)
    layout(matrix(1, ncol=1))
    par(oma=c(4, 5, 2, 1),
        mar=c(0.5, 0, 0.5, 1))
    plot.panel(dats, new.dd, y1, xs,
               col.lines, col.fill, plotPoints,
               agg.col=agg.col)
    axis(1, pretty(dats[,xs]), labels=pretty(dats[,xs]))
    mtext(ylabel, 2, line=3, cex=1)
    mtext(xlabel, 1, line=2.5, cex=1)
  }
  pdf.f(plot.ci, file=file.path(path,
                   sprintf('%s.pdf', paste(
                     gsub('[[:space:]]', '_', xlabel),
                     gsub('[[:space:]]', '_', ylabel),
                     extinction.method,  sep='_'))),
        width=4, height=4)

}
