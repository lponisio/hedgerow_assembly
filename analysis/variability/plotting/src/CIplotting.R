plot.panel <- function(dats,
                       new.dd,
                       y1,
                       y2,
                       xs,
                       treatments,
                       col.lines,
                       col.fill,
                       col.points=col.fill,
                       ylabel="",
                       ag.col="SiteStatus",
                       plot.x=TRUE,
                       scaled=FALSE,
                       plot.y=TRUE,
                       pchs=c(16, 16, 16),
                       dec=1){
  plotting.loop <- function(){
    for(i in 1:length(treatments)){
      print(treatments[i])
      ## subset data to specific treatment
      sub.dd <- new.dd[new.dd[,ag.col]==treatments[i],]
      sub.dats <- dats[dats[,ag.col]==treatments[i],]
      ys <- aggregate(list(y=sub.dats[,y1]),
                      list(sp=sub.dats[, "GenusSpecies"],
                           x=sub.dats[,xs]),
                      mean, na.rm=TRUE)
      ## plots means
      points(x=jitter(ys$x, factor=0.25),
             y=ys$y,
             col=col.points[treatments[i]],
             cex=1,
             pch=pchs[i])
      ## points(x=jitter(sub.dats[,xs], factor=0.25),
      ##      y=sub.dats[,y1],
      ##      col=col.fill[i],
      ##      cex=1,
      ##      pch=pchs[i])
      ## plots CI
      lines(x=sub.dd[,xs],
            y=sub.dd[,y1],
            col=col.lines[treatments[i]],
            lwd=2)
      lines(x=sub.dd[,xs],
            y=sub.dd$plo,
            col=col.lines[treatments[i]],
            lty="dashed")
      lines(x=sub.dd[,xs],
            y=sub.dd$phi,
            col=col.lines[treatments[i]],
            lty="dashed")
      ## add fill from ci.up to ci.lb
      polygon(c(sub.dd[,xs],
                rev(sub.dd[,xs])),
              c(sub.dd$phi,
                rev(sub.dd$plo)),
              col=col.fill[i], border=NA)
    }
  }

  plot(NA,
       xlim=range(new.dd[,xs], na.rm=TRUE),
       ylim=range(c(new.dd$phi,  new.dd$plo, y2,
          quantile(dats[,y1],
         na.rm=TRUE)[c("25%","75%")]),
         na.rm=TRUE),
       xlab="",
       ylab="",
       xaxt="n",
       yaxt="n",
       las=1)
  if(plot.y){
    axis(2,
         pretty(range(c(new.dd$phi,  new.dd$plo, y2), na.rm=TRUE),
                       min.n=3),
         las=1)
    mtext(ylabel, 2, line=3, cex=1.5)
  }
  if(plot.x){
    if(scaled){
      axis(1, pretty(new.dd[,xs], 3),
           labels= round((pretty(new.dd[,xs], 5)*
             sd(dats[, "traits.ns"], na.rm=TRUE)) +
             mean(dats[, "traits.ns"], na.rm=TRUE), dec))
    } else{
      axis(1, pretty(new.dd[,xs], 4))
    }
  }
  plotting.loop()
}



plot.predict.div <- function(new.dd,
                             ylabel,
                             dats,
                             y1,
                             y2=NA,
                             xs="",
                             legend.loc="bottomright",
                             legend.loc.year="topleft",
                             x.adj=-0.5,
                             width=8, height=5,
                             xlabel,
                             scaled=TRUE,
                             treatments = c("control", "maturing",
                               "mature"),
                             col.lines =  brewer.pal(3, "Dark2")[c(2,
                             1, 3)],
                             f.path = 'figures',
                             dec=1){
  plot.ci <- function(){
    col.fill <- add.alpha(col.lines, alpha=0.2)
    names(col.lines) <- names(col.fill) <- treatments
    layout(matrix(1, ncol=1))
    par(oma=c(6, 5, 2, 1),
        mar=c(0.5, 0, 0.5, 1))
    plot.panel(dats, new.dd, y1, y2, xs,
               treatments,
               col.lines,
               col.fill,
               ylabel,
               scaled=scaled,
               dec=dec)
    mtext(xlabel, 1, line=3.5, cex=1.5, adj=x.adj)
    if(length(treatments) > 1){
      legend(legend.loc,
             legend=c("Unrestored", "Maturing", "Mature"),
             col=col.lines,
             pch=16, bty="n", cex=1)
    }
  }
  pdf.f(plot.ci, file=file.path(f.path,
                   sprintf("%s.pdf", paste(
                     gsub(" ", "", ylabel),
                     gsub(" ", "", xlabel),
                     sep="_"))),
        width=width, height=height)

}
