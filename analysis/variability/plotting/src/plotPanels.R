plot.panels <- function(){
  f <- function(){
    col.lines <- "black"
    col.fill <- add.alpha(col.lines, alpha=0.2)
    names(col.lines) <- names(col.fill) <- "all"
    layout(matrix(1:4, nrow=2, byrow=TRUE))
    par(oma=c(2, 5, 0.5, 1),
        mar=c(2.5, 0, 0, 1), cex.axis=1.2)
    ## persistence pollinators
    plot.panel(new.dd=occ.pi.close,
               dats=occ.closeness.cv$data,
               y1="cv",
               y2= range(c(degree.pi$plo, degree.pi$phi)),
               xs="traits",
               col.fill=col.fill,
               col.lines=col.lines,
               plot.y=TRUE,
               plot.x=TRUE,
               treatments="all")
    ## degree pollinators
    plot.panel(new.dd=degree.pi,
               dats=degree.closeness.cv$data,
               y1="cv",
               y2=range(c(occ.pi.close$plo, occ.pi.close$phi)),
               xs="traits",
               col.fill=col.fill,
               col.lines=col.lines,
               plot.y=FALSE,
               plot.x=TRUE,
               treatments="all",
               dec=0)

    ## persistence plants
    plot.panel(new.dd=plants.occ.pi.close,
               dats=plants.occ.closeness.cv$data,
               y1="cv",
               y2= range(c(plants.degree.pi$plo, plants.degree.pi$phi)),
               xs="traits",
               col.fill=col.fill,
               col.lines=col.lines,
               plot.y=TRUE,
               treatments="all")
    
    mtext("Persistence", 1, line=3, cex=1.2)
    mtext("Closeness variability (log)", 2,
          line=3.2, cex=1.2, adj=-0.75)

    ## degree plants
    plot.panel(new.dd=plants.degree.pi,
               dats=plants.degree.closeness.cv$data,
               y2=range(c(plants.occ.pi.close$plo, plants.occ.pi.close$phi)),
               y1="cv",
               xs="traits",
               col.fill=col.fill,
               col.lines=col.lines,
               plot.y=FALSE,
               treatments="all",
               dec=0)
    mtext("Degree", 1, line=3, cex=1.2)
  }
  path <- 'figures/cv'
  pdf.f(f, file=file.path(path,
             sprintf("%s.pdf", "occ_degree")),
        width=5, height=4)
}

