library(abind)
library(FD)
library(lme4)
library(RColorBrewer)

fig <- function(dats, bee.syr, sites) {
  g <- function(){
    layout(matrix(1:3, 3, 1, byrow=TRUE),
           widths=c(1,1,1,1), heights=c(1,1,1,1))
    par(oma=c(2,2,0.1,0.1), mar=c(0.5,3,0.3,0.3),
        mgp=c(0,1,0), cex.axis=1)

    panel <- function(metric, axes, lab) {
      ## add points and standard errors
      vals <- list(c=dats[dats$SiteStatus=='control',],
                   h=dats[dats$SiteStatus=='maturing',],
                   m=dats[dats$SiteStatus=='mature',])

      means <- sapply(vals, function (x) mean(x[, metric],
                                              na.rm=TRUE))
      errors <- sapply(vals, function (x) se(x[, metric]))
      ranges <-  range(c(means + errors, means-errors))

      plot(NA, xlim=c(0.15,0.85),
           ylim=ranges,
           xlab='', ylab='',
           xaxt='n',
           yaxt="n",
           las=1)
      axis(2, pretty(ranges, n=4), las=1)
      
      x.coord <- c(1/4, 2/4, 3/4)
      cols <- brewer.pal(3, "Dark2")
      points(x=x.coord, y=means, pch=16, col=cols)
      arrows(x.coord, means-errors, x.coord, means+errors,
             code=0, angle=90, length=0.02, col=cols)

      ## add axes labels
      if(axes[1]==1)
        axis(1, at=x.coord,
             labels=c('Unrestored', 'Maturing', 'Mature'))
    }

    panel(metric='zNODF', axes=c(0,1), 'a')
    panel(metric='zmodG', axes=c(0,1), 'b')
    panel(metric='H2', axes=c(1,1), 'c')
    
    mtext('Nestedness',
          side=2, line=0.25, cex=1, at=0.85, outer=TRUE)
    mtext('Modularity',
          side=2, line=0.25, cex=1, at=0.5, outer=TRUE)
    mtext('Specialization',
          side=2, line=0.25, cex=1, at=0.15, outer=TRUE)
  }
  fn <- sprintf('metrics_%s_%s.pdf', bee.syr, sites)
  file.name <- file.path('figures', fn)
  pdf.f(g, file.name, height=4, width=2.5)
}

fig(dats=cor.dats, bee.syr="all", sites='total')
## ************************************************************
