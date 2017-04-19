library(RColorBrewer)
library(beeswarm)
## box plot

plot.beta.div <- function(dis.method, ## dissmiliarity method for file name
                          dists, ## beta-div estimate (box data)
                          status, ## vector of site statuses (box catagoties)
                          type, ## time or space for ylab
                          sub, ## bees/syrphids for file name
                          occ,## occ data? TRUE/FALSE
                          fig.path = 'figures/betadisper',
                          ylabel){
  makeBetaBox <- function(){

    par(oma=c(2,6,1,1), mar=c(5,0,2,0.5), mgp=c(2,1,0),
        cex.axis=1.5)
    status <- factor(status,
                     levels=c('control', 'maturing',
                       'mature'))
    cols <- brewer.pal(4, "Greys")[c(2,3,4)]
    cols.fills <- add.alpha(cols, alpha=0.5)
    boxplot(dists ~ status, col=cols.fills,
            xlab='', ylab='',
            names=c("","",""),
            las=1,
            ylim=c(0,1))
    mtext(c("Non-assembling \n field margin",
            "Assembling \n  hedgerow",
            "Non-assembling \n hedgerow"),
          side = 1, line= 2, at = 1:3)
    beeswarm(dists ~ status, col = cols, add = TRUE)

    mtext(ylabel,
          2, line=3.5, cex=1.5)

  }
  pdf.f(makeBetaBox,
        file= file.path(fig.path, sprintf('%s.pdf',
          paste(type,
                occ,
                dis.method,
                sub,
                sep='_'))),
        width=6, height=6)
}

## extracts model coefficients and SE and plots them as points and
## lines

plot.coeffs <- function(dis.method, ## dissmiliarity method for file name
                        mod, ## model object
                        status, ## vector of site statuses
                        type,  ## time or space for ylab
                        sub, ## bees/syrphids for file name
                        occ, ## occ data? TRUE/FALSE
                        fig.path = 'figures/betadisper',
                        add.labels=TRUE){
  ## cols <-brewer.pal(6, 'Dark2')[c(6,2,1)]
  cols <- brewer.pal(4, "Greys")[c(2,3,4)]
  f <- function(){
    par(oma=c(2,6,1,1), mar=c(5,0,2,0.5), mgp=c(2,1,0),
        cex.axis=1.5)
    means <- c(coef(mod)[1,1],
               coef(mod)[1,1] + coef(mod)[3,1],
               coef(mod)[1,1] + coef(mod)[2,1])
    ci.ub <- means +
      c(coef(mod)[1,2], coef(mod)[3,2], coef(mod)[2,2])
    ci.lb <- means -
      c(coef(mod)[1,2], coef(mod)[3,2], coef(mod)[2,2])
    plot(x=1:3, y=means,
         col=cols,
         pch=16,
         ylim=c(0,0.6),
         xlim=c(0.5,3.5),
         xlab='',
         xaxt='n',
         cex=1.5, las=1)
    arrows(x0=1:3,
           y0= ci.lb, y1=ci.ub, angle=90,
           length=0, code=3, col=cols, lwd=1.5)
    mtext("Species temporal turnover", 2, line=3.5 , cex=2)
    if(add.labels){
      axis(1, at= 1:3,
           labels=c("Non-assembling \n field margin",
             "Assembling \n  hedgerow",
             "Non-assembling \n hedgerow"))
    }
  }
  pdf.f(f, file= file.path(fig.path, sprintf('%s.pdf',
             paste('coeff', type,
                   occ,
                   dis.method,
                   sub,
                   sep='_'))),
        width=6, height=6)
}


plot.node.box <- function(ylabel,
                          dats,
                          y1){
  makeNodeBox <- function(){
    par(oma=c(2,6,1,1), mar=c(5,0,2,0.5), mgp=c(2,1,0),
        cex.axis=1.5)
    cols <- brewer.pal(4, "Greys")[c(2,3,4)]
    cols.fills <- add.alpha(cols, alpha=0.5)
    boxplot(dats[,y1]~ dats$SiteStatus,
            col=cols.fills,
            names=c("","",""))
    mtext(c("Non-assembling \n field margin",
            "Assembling \n  hedgerow",
            "Non-assembling \n hedgerow"),
          side = 1, line= 2, at = 1:3)
    beeswarm(dats[,y1]~ dats$SiteStatus, col = cols, add = TRUE)
    mtext(ylabel, 2, line=3.5, cex=1.5)
  }
  path <- 'figures'
  pdf.f(makeNodeBox, file=file.path(path,
                       sprintf("%s.pdf", paste(
                         gsub(" ", "", ylabel),
                         "box", sep="_"))),
        width=6, height=6)

}


plotDistPanels <- function(){
  f3 <- function(){
    layout(matrix(1:4, ncol=2, byrow=TRUE))
    par(oma=c(2.5, 6, 1, 1),
        mar=c(1, 1, 2, 1), cex.axis=1.5)
    cols <- brewer.pal(4, "Greys")[c(2,3,4)]
    cols.fills <- add.alpha(cols, alpha=0.5)
    ## pollinator species turnover
    load(file= file.path('saved/speciesTurnover', sprintf('%s.pdf',
           paste(dis.method, alpha, occ, type="pols", sep='_'))))
    dists <- dats$dist
    status <- dats$status
    status <- factor(status,
                     levels=c('control', 'maturing',
                       'mature'))
    boxplot(dists ~ status, col=cols.fills,
            xlab='', ylab='', ylim=c(0,1),
            names=c("","",""),
            las=1, xaxt="n")
    text(x=1:3, y=0.9, c("a", "a", "b"))
    beeswarm(dists ~ status, col = cols.fills,
             add = TRUE, method="hex", corral="gutter")
    mtext("Pollinators",
          3, line=0.9, cex=1.5)

    ## plant species turnover
    load(file= file.path('saved/speciesTurnover', sprintf('%s.pdf',
           paste(dis.method, alpha, occ, type="plants", sep='_'))))
    dists <- dats$dist
    status <- dats$status
    status <- factor(status,
                     levels=c('control', 'maturing',
                       'mature'))
    boxplot(dists ~ status, col=cols.fills,
            xlab='', ylab='',ylim=c(0,1),
            names=c("","",""),
            las=1, xaxt="n", yaxt="n")
    beeswarm(dists ~ status, col = cols.fills, add = TRUE,
             method="hex", corral="gutter")
    mtext("Plants",
          3, line=0.9, cex=1.5)

    ## interaction turnover
    load(file= file.path('saved/speciesTurnover', sprintf('%s.pdf',
           paste(dis.method, alpha, occ, type="ints", sep='_'))))
    dists <- dats$dist
    status <- dats$status
    status <- factor(status,
                     levels=c('control', 'maturing',
                       'mature'))
    boxplot(dists ~ status, col=cols.fills,
            xlab='', ylab='',ylim=c(0,1),
            names=c("","",""),
            las=1)
    beeswarm(dists ~ status, col = cols.fills, add = TRUE,
             method="hex", corral="gutter")
    mtext("Interactions",
          3, line=0.9, cex=1.5)
    mtext(c("Non-assembling \n field margin",
            "Assembling \n  hedgerow",
            "Non-assembling \n hedgerow"),
          side = 1, line= 2, at = 1:3, cex=0.9)
    mtext("Turnover", 2, line=4, cex=2, adj=1.35)

    ## weighted link turnover
    load(file="saved/phyloInt.Rdata")
    dats <- phylo.int$phylo.int
    y1 <- "PhyloInt"
    boxplot(dats[,y1]~ dats$SiteStatus,
            col=cols.fills,
            names=c("","",""), las=1,
            ylim=c(0,1), yaxt="n")
    mtext(c("Non-assembling \n field margin",
            "Assembling \n  hedgerow",
            "Non-assembling \n hedgerow"),
          side = 1, line= 2, at = 1:3, cex=0.9)
    text(x=1:3, y=0.2, c("a", "b", "b"))
    beeswarm(dats[,y1]~ dats$SiteStatus, col = cols.fills, add = TRUE,
             method="hex", corral="gutter")
    mtext("Weighted interactions", 3, line=0.9, cex=1.5)
  }
  path <- 'figures'
  pdf.f(f3, file=file.path(path,
              sprintf("%s.pdf", "turnover_panels")),
        width=10, height=9)
}

