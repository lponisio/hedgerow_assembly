plot.panels <- function(){
  f <- function(){
    col.lines <-  brewer.pal(4, "Greys")[3]
    col.fill <- add.alpha(col.lines, alpha=0.2)
    layout(matrix(1:2, ncol=1))
    par(oma=c(4, 6, 0.5, 1),
        mar=c(0.5, 0, 2, 1), cex.axis=1.5)
    ## nodf
    plot.panel(new.dd=ypr.pi.pol,
               dats=specs,
               y1="closeness",
               xs="ypr",
               col.fill=col.fill,
               col.lines=col.lines,
               plotPoints=TRUE,
               agg.col="GenusSpecies")

    mtext("Pollinators", 3, line=0.5, cex=1.5)

    plot.panel(new.dd=ypr.pi.plant,
               dats=specs,
               y1="closeness",
               xs="ypr",
               col.fill=col.fill,
               col.lines=col.lines,
               plotPoints=TRUE,
               agg.col="GenusSpecies")
    
    axis(1, pretty(specs$ypr), labels=pretty(specs$ypr))
    mtext("Closeness", 2, line=4.5, cex=1.5, adj=1.75)
    mtext("Plants", 3, line=0.5, cex=1.5)
    mtext("Years of assembly", 1, line=3, cex=1.5)

  }
  path <- 'figures' 
  pdf.f(f, file=file.path(path,
             sprintf("%s.pdf", "closenessPanel")),
        width=4, height=6)

}

