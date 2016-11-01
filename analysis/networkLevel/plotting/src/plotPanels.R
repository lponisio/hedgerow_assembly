plot.panels <- function(){
  f <- function(){
    col.white <- add.alpha("white", alpha=0)
    col.lines <- brewer.pal(4, "Greys")[3]
    col.fill <- add.alpha(col.lines, alpha=0.2)
    layout(matrix(1:6, ncol=3))
    par(oma=c(6, 1, 2, 1),
        mar=c(1, 6, 0.5, 1), cex.axis=1.5)
    ## nodf
    plot.panel(new.dd=nodf.pi,
               dats=cor.dats,
               y1="zNODF",
               xs="ypr",
               col.fill=col.fill,
               col.lines=col.lines,
               plotPoints=TRUE)
    mtext("Nestedness", 2, line=4, cex=1.5)
    ## modularity
    plot.panel(new.dd=mod.pi,
               dats=cor.dats,
               y1="zmod.met.D",
               xs="ypr",
               col.fill=col.white,
               col.lines=col.white,
               col.points=col.fill,
               plotPoints=TRUE)
    mtext("Modularity", 2, line=4, cex=1.5)

    axis(1, pretty(cor.dats$ypr), labels=pretty(cor.dats$ypr))
    mtext("Year of assembly", 1, line=3.5, cex=1.5)
    
    ## specialization
    plot.panel(new.dd=h2.pi,
               dats=cor.dats,
               y1="zH2",
               xs="ypr",
               col.fill=col.fill,
               col.lines=col.lines,
               plotPoints=TRUE)
    mtext("Specialization", 2, line=4, cex=1.5)

    plot.panel(new.dd=conn.pi,
               dats=cor.dats,
               y1="connectance",
               xs="ypr",
               col.fill=col.fill,
               col.lines=col.lines,
               plotPoints=TRUE)
    mtext("Connectance", 2, line=4, cex=1.5)
    axis(1, pretty(cor.dats$ypr), labels=pretty(cor.dats$ypr))
    mtext("Year of assembly", 1, line=3.5, cex=1.5)

    plot.panel(new.dd=rich.hl.pi,
               dats=cor.dats,
               y1="number.of.species.HL",
               xs="ypr",
               col.fill=col.fill,
               col.lines=col.lines,
               plotPoints=TRUE)
    mtext("Pollinator richness", 2, line=4, cex=1.5)

    plot.panel(new.dd=rich.ll.pi,
               dats=cor.dats,
               y1="number.of.species.LL",
               xs="ypr",
               col.fill=col.fill,
               col.lines=col.lines,
               plotPoints=TRUE)
    mtext("Plant richness", 2, line=4, cex=1.5)
    
    axis(1, pretty(cor.dats$ypr), labels=pretty(cor.dats$ypr))
    mtext("Year of assembly", 1, line=3.5, cex=1.5)
  }
  path <- 'figures'
  pdf.f(f, file=file.path(path,
             sprintf("%s.pdf", "baci")),
        width=10, height=6)

}

