plot.panels <- function(){
  f <- function(){
    col.lines <- brewer.pal(3, "Dark2")[2]
    col.fill <- add.alpha(col.lines, alpha=0.2)
    layout(matrix(1:2, ncol=2))
    par(oma=c(6, 7, 2, 1),
        mar=c(0.5, 0, 0.5, 1), cex.axis=1.5)
    ## nodf
    plot.panel(new.dd=nodf.pi,
               dats=cor.dats,
               y1="Robustness",
               xs="ypr",
               col.fill=col.fill,
               col.lines=col.lines,
               plotPoints=TRUE)
    mtext("Robustness", 2, line=5, cex=1.5)
    axis(1, pretty(cor.dats$ypr), labels=pretty(cor.dats$ypr))
    mtext("Years Post Restoration", 1, line=3.5, cex=1.5)

  }
  path <- 'figures' 
  pdf.f(f, file=file.path(path,
             sprintf("%s.pdf", "baci")),
        width=5, height=12)

}

