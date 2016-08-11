plot.panels <- function(){
  f <- function(){
    col.lines <- brewer.pal(4, "Greys")[3]
    col.fill <- add.alpha(col.lines, alpha=0.2)
    layout(matrix(1:2, ncol=1))
    par(oma=c(4, 1, 2, 1),
        mar=c(1, 6, 0.5, 1), cex.axis=1.5)
    ## nodf
    plot.panel(new.dd=ypr.pi,
               dats=res,
               y1="Robustness",
               xs="ypr",
               col.fill=col.fill,
               col.lines=col.lines,
               plotPoints=TRUE)
    mtext("Robustness to \n species extinction", 2, line=4, cex=1.5)
    ## modularity
    plot.panel(new.dd=ypr.pi.alg,
               dats=all.alg.Con.status,
               y1="AlgCon",
               xs="ypr",
               col.fill=col.fill,
               col.lines=col.lines,
               plotPoints=TRUE)
    mtext("Robustness to \n perturbation", 2, line=4, cex=1.5)

    axis(1, pretty(res$ypr), labels=pretty(res$ypr))
    mtext("Year of assembly", 1, line=3.5, cex=1.5)
  }
  path <- 'figures'
  pdf.f(f, file=file.path(path,
             sprintf("%s.pdf", "robustness")),
        width=4, height=5.5)

}

