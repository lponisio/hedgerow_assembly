plot.panels <- function(){
  f <- function(){
    col.lines <- "black"
    col.fill <- add.alpha(col.lines, alpha=0.2)
    names(col.lines) <- names(col.fill) <- "all"
    layout(matrix(1:4, nrow=2, byrow=TRUE))
    par(oma=c(4, 5, 0.5, 1),
        mar=c(0.5, 0, 0, 1), cex.axis=1.2)
    ## persistence pollinators
    plot.panel(new.dd=occ.pi,
               dats=occ.closeness.cv$data,
               y1="cv",
               y2=60,
               xs="traits",
               col.fill=col.fill,
               col.lines=col.lines,
               plot.y=TRUE,
               treatments="all")
    ## mtext("Pollinator persistence", 1, line=2.5, cex=1.2)
    mtext("Closeness variability", 2, line=3, cex=1.2)
    ## degree pollinators
    plot.panel(new.dd=degree.pi,
               dats=degree.closeness.cv$data,
               y1="cv",
               y2=60,
               xs="traits",
               col.fill=col.fill,
               col.lines=col.lines,
               plot.y=FALSE,
               treatments="all",
               dec=0)
    ## mtext("Pollinator degree", 1, line=2.5, cex=1.2)

        plot.panel(new.dd=occ.pi,
               dats=occ.closeness.cv$data,
               y1="cv",
               y2=60,
               xs="traits",
               col.fill=col.fill,
               col.lines=col.lines,
               plot.y=TRUE,
               treatments="all")
    mtext("Pollinator persistence", 1, line=2.5, cex=1.2)
    mtext("Closeness variability", 2, line=3, cex=1.2)
    ## degree pollinators
    plot.panel(new.dd=degree.pi,
               dats=degree.closeness.cv$data,
               y1="cv",
               y2=60,
               xs="traits",
               col.fill=col.fill,
               col.lines=col.lines,
               plot.y=FALSE,
               treatments="all",
               dec=0)
    mtext("Pollinator degree", 1, line=2.5, cex=1.2)
  }
  path <- 'figures/cv'
  pdf.f(f, file=file.path(path,
             sprintf("%s.pdf", "occ_degree")),
        width=6.2, height=3)
}

