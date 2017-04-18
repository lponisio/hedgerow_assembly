plot.panels <- function(){
  f <- function(){
    col.lines <- "black"
    col.fill <- add.alpha(col.lines, alpha=0.2)
    col.white <- add.alpha("white", alpha=0)
    names(col.lines) <- names(col.fill) <- "all"
    layout(matrix(1:4, nrow=2, byrow=TRUE))
    par(oma=c(2, 6.5, 0.5, 1),
        mar=c(2.5, 0, 0, 1), cex.axis=1.2)

    ## persistence pollinators
    plot.panel(new.dd=pol.occ.pi.close,
               dats=pol.cv$lm.data,
               y1="cv",
               y2= range(c(pol.degree.pi.close$plo, pol.degree.pi.close $phi)),
               xs="occ.date",
               col.fill=col.fill,
               col.lines=col.lines,
               plot.y=TRUE,
               plot.x=TRUE,
               treatments="all")
    mtext("Pollinators", 2, line=5, cex=1.2)

    ## degree pollinators
    plot.panel(new.dd=pol.degree.pi.close,
               dats=pol.cv$lm.data,
               y1="cv",
               y2=range(c(pol.occ.pi.close$plo, pol.occ.pi.close$phi)),
               xs="r.degree",
               col.fill=col.white,
               col.lines=col.white,
               col.points=col.fill,
               plot.y=FALSE,
               plot.x=TRUE,
               treatments="all",
               dec=0)
    ## persistence plants
    plot.panel(new.dd=plant.occ.pi.close,
               dats=plant.cv$lm.data,
               y1="cv",
               y2= range(c(plant.degree.pi.close$plo, plant.degree.pi.close$phi)
    + c(-0.2,0.2)),
               xs="occ.plant.date",
               col.fill=col.white,
               col.lines=col.white,
               col.points=col.fill,
               plot.y=TRUE,
               treatments="all")

    mtext("Plants", 2, line=5, cex=1.2)

    mtext("Persistence", 1, line=3, cex=1.2)
    mtext("Closeness variability (log)", 2,
          line=3.2, cex=1.2, adj=-0.75)


    ## degree plants
    plot.panel(new.dd=plant.degree.pi.close,
               dats=plant.cv$lm.data,
               y2=range(c(plant.occ.pi.close$plo,
    plant.occ.pi.close$phi)) + c(-0.2, 0.1),
               y1="cv",
               xs="plant.r.degree",
               col.fill=col.white,
               col.lines=col.white,
               col.points=col.fill,
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

