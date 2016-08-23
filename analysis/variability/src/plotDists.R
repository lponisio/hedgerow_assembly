library(RColorBrewer)
library(beeswarm)
## box plot

f1 <- function(){
  makeBetaBox <- function(dists=dists, ## beta-div estimate (box data)
                          status=status){ ## vector of site statuses
    ## (box catagoties)
    status <- factor(status,
                     levels=c('control', 'maturing',
                       'mature'))
    cols <- brewer.pal(4, "Greys")[c(2,3,4)]
    cols.fills <- add.alpha(cols, alpha=0.5)
    boxplot(dists ~ status, col=cols.fills,
            xlab='', ylab='',
            names=c("","",""),
            las=1)
    if(add.labels){
      mtext(c("Non-assembling \n field margin",
              "Assembling \n  hedgerow",
              "Non-assembling \n hedgerow"),
            side = 1, line= 2, at = 1:3)
    }
    beeswarm(dists ~ status, col = cols, add = TRUE)
    if(type == 'time'){
      mtext("Species temporal turnover",
            2, line=3.5, cex=1.5)
    } else {
      mtext(expression(paste(beta-'diversity', ' (corrected)')),
            2, line=3.3, cex=2)
    }
  }
}

plot.beta.div <- function(dis.method, ## dissmiliarity method for file name
                          dists, ## beta-div estimate (box data)
                          status, ## vector of site statuses (box catagoties)
                          type, ## time or space for ylab
                          sub, ## bees/syrphids for file name
                          occ,## occ data? TRUE/FALSE
                          fig.path = 'figures/betadisper'){
  par(oma=c(2,6,1,1), mar=c(5,0,2,0.5), mgp=c(2,1,0),
      cex.axis=1.5)

  pdf.f(f1,
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


f2 <- function(){
  makeNodeBox <- function(ylabel=ylabel,
                          dats=dats,
                          y1=y1){
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
}

plot.node.box <- function(ylabel,
                          dats,
                          y1){
  par(oma=c(2,6,1,1), mar=c(5,0,2,0.5), mgp=c(2,1,0),
      cex.axis=1.5)
  path <- 'figures'
  pdf.f(f2, file=file.path(path,
              sprintf("%s.pdf", paste(
                gsub(" ", "", ylabel),
                "box", sep="_"))),
        width=6, height=6)
  
}


plotDistPanels <- function(){
  f3 <- function(){
    layout(matrix(1:2, nrow=2))
    par(oma=c(3, 5, 0.5, 1),
        mar=c(1, 0, 0, 1), cex.axis=1.2)
    cols <- brewer.pal(4, "Greys")[c(2,3,4)]
    cols.fills <- add.alpha(cols, alpha=0.5)
    ## species turnover
    dists <- dats.pols$dist
    status <- dats.pols$status
    status <- factor(status,
                     levels=c('control', 'maturing',
                       'mature'))
    boxplot(dists ~ status, col=cols.fills,
            xlab='', ylab='',
            names=c("","",""),
            las=1, xaxt="n")
    beeswarm(dists ~ status, col = cols, add = TRUE)
    mtext("Species turnover",
          2, line=3.5, cex=1.5)
    
    ## node turnover
    load(file="saved/phyloInt.Rdata")
    ylabel <- "Weighted species turnover"
    dats <- phylo.int$phylo.int
    y1 <- "PhyloInt"
    boxplot(dats[,y1]~ dats$SiteStatus,
            col=cols.fills,
            names=c("","",""), las=1)
    mtext(c("Non-assembling \n field margin",
            "Assembling \n  hedgerow",
            "Non-assembling \n hedgerow"),
          side = 1, line= 2, at = 1:3, cex=1.2)
    beeswarm(dats[,y1]~ dats$SiteStatus, col = cols, add = TRUE)
    mtext(ylabel, 2, line=3.5, cex=1.5)
  }
  path <- 'figures'
  pdf.f(f3, file=file.path(path,
              sprintf("%s.pdf", "turnover_panels")),
        width=8, height=11)
}

