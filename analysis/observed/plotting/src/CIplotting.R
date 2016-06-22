plot.panel <- function(dats,
                       new.dd,
                       y1,
                       xs,
                       col.lines,
                       col.fill,
                       plotPoints=TRUE,...){
  plotting.loop <- function(){
    if("Site" %in% colnames(dats)){
      ys <- aggregate(list(y=dats[,y1]),
                      list(Site=dats$Site,
                           x=dats[,xs]),
                      mean, na.rm=TRUE)
    } else{
      ys <- dats
    }
    if(plotPoints){
      points(x=jitter(ys$x, factor=0.25),
             y=ys$y,
             pch=16,
             col=col.lines,
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
       ylim=range(c(0, new.dd$phi,  new.dd$plo,
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
                             xs="ypr",
                             plotPoints=TRUE,
                             xlabel= "Years Post Restoration"){
  plot.ci <- function(){
    col.lines <-  brewer.pal(3, "Dark2")[2]
    col.fill <- add.alpha(col.lines, alpha=0.2)
    layout(matrix(1, ncol=1))
    par(oma=c(6, 5, 2, 1),
        mar=c(0.5, 0, 0.5, 1))
    plot.panel(dats, new.dd, y1, xs,
               col.lines, col.fill, plotPoints)
    axis(1, pretty(dats[,xs]), labels=pretty(dats[,xs]))
    mtext(ylabel, 2, line=3.5, cex=1.5)
    mtext(xlabel, 1, line=3.5, cex=1.5)
  }
  path <- 'figures'
  pdf.f(plot.ci, file=file.path(path,
                                sprintf("%s.pdf", paste(
                                                    gsub(" ", "", ylabel),
                                                    xlabel, sep="_"))),
        width=4, height=4)

}


plotSpecs <- function(metric, diff.spec, type, path.f){
  f <- function(){
    load(file=sprintf("saved/mods/%s%s.Rdata", type, metric))
    y <- diff.spec[, paste(metric, "y", sep=".")]
    x <- diff.spec[, paste(metric, "x", sep=".")]

    dd.mod <- data.frame(cbind(x=seq(min(x), max(x), length=10),
                               y=seq(min(y), max(y), length=10)))

    mod.pi <- try(predict.int(mod= out.mod,
                              dd=dd.mod,
                              y="y",
                              family="gaussian"), silent=TRUE)
    if(inherits(mod.pi, "try-error")) {
      return(NA)
    }else{
      dats <- data.frame(cbind(x=x, y=y))
      col.lines <- "darkred"
      col.fill <- add.alpha(col.lines, alpha=0.2)
      layout(matrix(1, ncol=1))
      par(oma=c(6, 5, 2, 1),
          mar=c(0.5, 0, 0.5, 1))
      plot.panel(new.dd=mod.pi,
                 dats=dats,
                 y1="y",
                 xs="x",
                 col.fill=col.fill,
                 col.lines=col.lines,
                 plotPoints=TRUE)
      axis(1, pretty(dats$y))
      mtext("Generalization early assembly", 1, line=3.5, cex=1.5)
      mtext("Generalization late assembly", 2, line=3.5, cex=1.5)
      abline(a=0, b=1, lty="dashed")
    }
  }
  pdf.f(f, file=file.path(path.f,
                          sprintf("%s/%s.pdf", type, metric)),
        width=5, height=5)
}
