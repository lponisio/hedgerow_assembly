



plot.panels <- function(){
  f <- function(){
    col.lines <- c(rev(brewer.pal(4, "RdYlGn"))[c(2,3,4)], "black")
    col.fill <- add.alpha(col.lines, alpha=0.3)
    treatments <- c("LOW", "MOD", "HIGH")

    layout(matrix(1:6, ncol=3))
    par(oma=c(6, 7, 3, 1),
        mar=c(0.5, 0, 0.5, 1), cex.axis=1.5)
    ## ************************************************************
    ## floral across site statuses

    plot(NA, xlim=range(scale(floral.dist.nostat$dist)),
         ylim=c(0, 1),
         xlab="",
         ylab="",
         xaxt="n",
         yaxt="n",
         las=1)
    mtext("Flowering plants", 3, line=1, cex=1.5)

    points(x=jitter(scale(floral.dist.nostat$dist), factor=0.25),
           y=floral.dist.nostat$comm,
           pch=16,
           col=col.fill[4],
           cex=1.5)
    ## plots CI
    lines(x=floral.pi.nostat$dist,
          y=floral.pi.nostat$comm,
          col=col.lines[4],
          lwd=2)
    lines(x=floral.pi.nostat$dist,
          y=floral.pi.nostat$plo,
          col=col.lines[4],
          lty="dashed")
    lines(x=floral.pi.nostat$dist,
          y=floral.pi.nostat$phi,
          col=col.lines[4],
          lty="dashed")
    ## add fill from ci.up to ci.lb
    polygon(c(floral.pi.nostat$dist,
              rev(floral.pi.nostat$dist)),
            c(floral.pi.nostat$phi,
              rev(floral.pi.nostat$plo)),
            col=col.fill[4], border=NA)
    axis(2, pretty(c(-0.3, 1), 4), las=1)
    ## ************************************************************
    ## floral within a site status
    plot.panel(dats=floral.dist,
               new.dd=floral.pi,
               xs="dist",
               y1="comm",
               y2="comm",
               year=NA,
               treatments=treatments,
               col.lines=col.lines,
               col.fill=col.fill,
               ylabel=NA,
               plot.x=TRUE)
    ## mtext("Geographic Distance (km)", 1, line=3.5, cex=1.5, at=2)
    mtext("Corrected Dissimilarity", 2, line=4, at=1.2, cex=1.5)

    ## ************************************************************
    ## pollinators between
    plot(NA, xlim=range(scale(pol.dist.nostat$dist)),
         ylim=c(0, 1),
         xlab="",
         ylab="",
         xaxt="n",
         yaxt="n",
         las=1)
    mtext("Bees", 3, line=1, cex=1.5)

    points(x=jitter(scale(pol.dist.nostat$dist), factor=0.25),
           y=pol.dist.nostat$comm,
           pch=16,
           col=col.fill[4],
           cex=1.5)
    ## plots CI
    lines(x=pol.pi.nostat$dist,
          y=pol.pi.nostat$comm,
          col=col.lines[4],
          lwd=2)
    lines(x=pol.pi.nostat$dist,
          y=pol.pi.nostat$plo,
          col=col.lines[4],
          lty="dashed")
    lines(x=pol.pi.nostat$dist,
          y=pol.pi.nostat$phi,
          col=col.lines[4],
          lty="dashed")
    ## add fill from ci.up to ci.lb
    polygon(c(pol.pi.nostat$dist,
              rev(pol.pi.nostat$dist)),
            c(pol.pi.nostat$phi,
              rev(pol.pi.nostat$plo)),
            col=col.fill[4], border=NA)
    ## ************************************************************
    ## pollinators within
    plot.panel(dats=pol.dist,
               new.dd=pol.pi,
               xs="dist",
               y1="comm",
               y2="comm",
               year=NA,
               treatments=treatments,
               col.lines=col.lines,
               col.fill=col.fill,
               ylabel=NA,
               plot.x=TRUE,
               plot.y=FALSE)
    mtext("Geographic Distance (km)", 1, line=3.5, cex=1.5)
    ## ************************************************************
    ## interactions between
    plot(NA, xlim=range(scale(int.dist.nostat$dist)),
         ylim=c(0, 1),
         xlab="",
         ylab="",
         xaxt="n",
         yaxt="n",
         las=1)
    mtext("Interactions", 3, line=1, cex=1.5)

    points(x=jitter(scale(int.dist.nostat$dist), factor=0.25),
           y=int.dist.nostat$comm,
           pch=16,
           col=col.fill[4],
           cex=1.5)
    ## plots CI
    lines(x=int.pi.nostat$dist,
          y=int.pi.nostat$comm,
          col=col.lines[4],
          lwd=2)
    lines(x=int.pi.nostat$dist,
          y=int.pi.nostat$plo,
          col=col.lines[4],
          lty="dashed")
    lines(x=int.pi.nostat$dist,
          y=int.pi.nostat$phi,
          col=col.lines[4],
          lty="dashed")
    ## add fill from ci.up to ci.lb
    polygon(c(int.pi.nostat$dist,
              rev(int.pi.nostat$dist)),
            c(int.pi.nostat$phi,
              rev(int.pi.nostat$plo)),
            col=col.fill[4], border=NA)
    ## ************************************************************
    ## interactions within
    plot.panel(dats=int.dist,
               new.dd=int.pi,
               xs="dist",
               y1="comm",
               y2="comm",
               year=NA,
               treatments=treatments,
               col.lines=col.lines,
               col.fill=col.fill,
               ylabel=NA,
               plot.x=TRUE,
               plot.y=FALSE)

    legend("bottomright",
           legend=c("Low", "Moderate","High"),
           col=col.lines,
           pch=c(15, 16, 17), cex=1.2)
  }
  path <- 'figures'
  pdf.f(f, file=file.path(path,
             sprintf("%s%s.pdf", "alldecay", dis.method)),
        width=8, height=6)

}


## plot.panels <- function(){
##   f <- function(){
##     col.lines <- c("darkolivegreen3",
##                    "darkgoldenrod1", "brown3", "black")
##     col.fill <- add.alpha(col.lines, alpha=0.2)
##     treatments <- c("LOW", "MOD", "HIGH")

##     layout(matrix(1:2, ncol=1))
##     par(oma=c(6, 7, 2, 1),
##         mar=c(0.5, 0, 0.5, 1), cex.axis=1.5)

##     ## across site statuses

##     plot(NA, xlim=range(scale(floral.dist.nostat$dist)),
##          ylim=range(c(-0.3, 1)),
##          xlab="",
##          ylab="",
##          xaxt="n",
##          yaxt="n",
##          las=1)

##     points(x=jitter(scale(floral.dist.nostat$dist), factor=0.25),
##            y=floral.dist.nostat$comm,
##            pch=16,
##            col=col.fill[4],
##            cex=1.5)
##     ## plots CI
##     lines(x=floral.pi.nostat$dist,
##           y=floral.pi.nostat$comm,
##           col=col.lines[4],
##           lwd=2)
##     lines(x=floral.pi.nostat$dist,
##           y=floral.pi.nostat$plo,
##           col=col.lines[4],
##           lty="dashed")
##     lines(x=floral.pi.nostat$dist,
##           y=floral.pi.nostat$phi,
##           col=col.lines[4],
##           lty="dashed")
##     ## add fill from ci.up to ci.lb
##     polygon(c(floral.pi.nostat$dist,
##               rev(floral.pi.nostat$dist)),
##             c(floral.pi.nostat$phi,
##               rev(floral.pi.nostat$plo)),
##             col=col.fill[4], border=NA)
##     axis(2, pretty(c(-0.3, 1), 4), las=1)

##     plot.panel(dats=floral.dist,
##                new.dd=floral.pi,
##                xs="dist",
##                y1="comm",
##                y2="comm",
##                year=NA,
##                treatments=treatments,
##                col.lines=col.lines,
##                col.fill=col.fill,
##                ylabel=NA,
##                plot.x=TRUE)
##     legend("bottomright",
##            legend=c("Low", "Moderate","High"),
##            col=col.lines,
##            pch=16, bty="n", cex=0.9)
##     mtext("Geographic Distance (km)", 1, line=3.5, cex=1.5)
##     mtext("Corrected Dissimilarity", 2, line=4, at=1.2, cex=1.5)

##   }
##   path <- 'figures'
##   pdf.f(f, file=file.path(path,
##              sprintf("%s.pdf", "poldecay")),
##         width=6, height=8)

## }

## plot.panels <- function(){
##   f <- function(){
##     col.lines <- c("darkolivegreen3",
##                    "darkgoldenrod1", "brown3", "black")
##     col.fill <- add.alpha(col.lines, alpha=0.2)
##     treatments <- c("LOW", "MOD", "HIGH")

##     layout(matrix(1:2, ncol=1))
##     par(oma=c(6, 7, 2, 1),
##         mar=c(0.5, 0, 0.5, 1), cex.axis=1.5)

##     ## across site statuses

##     plot(NA, xlim=range(scale(pol.dist.nostat$dist)),
##          ylim=range(c(-0.3, 1)),
##          xlab="",
##          ylab="",
##          xaxt="n",
##          yaxt="n",
##          las=1)

##     points(x=jitter(scale(pol.dist.nostat$dist), factor=0.25),
##            y=pol.dist.nostat$comm,
##            pch=16,
##            col=col.fill[4],
##            cex=1.5)
##     ## plots CI
##     lines(x=pol.pi.nostat$dist,
##           y=pol.pi.nostat$comm,
##           col=col.lines[4],
##           lwd=2)
##     lines(x=pol.pi.nostat$dist,
##           y=pol.pi.nostat$plo,
##           col=col.lines[4],
##           lty="dashed")
##     lines(x=pol.pi.nostat$dist,
##           y=pol.pi.nostat$phi,
##           col=col.lines[4],
##           lty="dashed")
##     ## add fill from ci.up to ci.lb
##     polygon(c(pol.pi.nostat$dist,
##               rev(pol.pi.nostat$dist)),
##             c(pol.pi.nostat$phi,
##               rev(pol.pi.nostat$plo)),
##             col=col.fill[4], border=NA)
##     axis(2, pretty(c(-0.3, 1), 4), las=1)

##     plot.panel(dats=pol.dist,
##                new.dd=pol.pi,
##                xs="dist",
##                y1="comm",
##                y2="comm",
##                year=NA,
##                treatments=treatments,
##                col.lines=col.lines,
##                col.fill=col.fill,
##                ylabel=NA,
##                plot.x=TRUE)
##     legend("bottomright",
##            legend=c("Low", "Moderate","High"),
##            col=col.lines,
##            pch=16, bty="n", cex=0.9)
##     mtext("Geographic Distance (km)", 1, line=3.5, cex=1.5)
##     mtext("Corrected Dissimilarity", 2, line=4, at=1.2, cex=1.5)

##   }
##   path <- 'figures'
##   pdf.f(f, file=file.path(path,
##              sprintf("%s.pdf", "poldecay")),
##         width=6, height=8)

## }
