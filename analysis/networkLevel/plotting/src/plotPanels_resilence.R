plot.panels <- function(){
  f <- function(){
    col.white <- add.alpha("white", alpha=0)
    col.lines <- brewer.pal(4, "Greys")[3]
    col.fill <- add.alpha(col.lines, alpha=0.2)
    cols.points <- add.alpha(brewer.pal(5, "Set1"), alpha=0.6)
    layout(matrix(1:2, ncol=1))
    par(oma=c(4, 1, 2, 1),
        mar=c(1, 6, 0.5, 1), cex.axis=1.5)
    ## species extinction
    plot.panel(new.dd=ypr.pi,
               dats=res,
               y1="Robustness",
               xs="ypr",
               col.fill=col.white,
               col.lines=col.white,
               col.points=cols.points,
               plotPoints=TRUE)
    mtext("Robustness to \n species extinction", 2, line=4, cex=1.5)
    ## sensitive to perturbations
    plot.panel(new.dd=ypr.pi.alg,
               dats=all.alg.Con.status,
               y1="AlgCon",
               xs="ypr",
               col.fill=col.white,
               col.lines=col.white,
               col.points=cols.points,
               plotPoints=TRUE)
    mtext("Sensitivity to \n perturbations", 2, line=4, cex=1.5)

    axis(1, pretty(res$ypr), labels=pretty(res$ypr))
    mtext("Year of assembly", 1, line=3.5, cex=1.5)
  }
  path <- 'figures'
  pdf.f(f, file=file.path(path,
             sprintf("%s.pdf", "robustness")),
        width=5, height=6.5)

}

