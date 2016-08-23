library(RColorBrewer)

## box plot 
plot.beta.div <- function(dis.method, ## dissmiliarity method for file name
                          dists, ## beta-div estimate (box data)
                          status, ## vector of site statuses (box catagoties)
                          type, ## time or space for ylab
                          sub, ## bees/syrphids for file name
                          occ,## occ data? TRUE/FALSE
                          fig.path <- 'figures/betadisper'){ 
  cols <- brewer.pal(4, "Greys")[c(2,3,4)]
  f <- function(){
    par(oma=c(2,6,1,1), mar=c(5,0,2,0.5), mgp=c(2,1,0),
        cex.axis=1.5)
    status <- factor(status,
                     levels=c('control', 'maturing',
                       'mature'))
    boxplot(dists ~ status, col=cols,
            xlab='', ylab='',
            names=c('Unrestored', 'Maturing',
              'Mature'),
            las=1)
    if(type == 'time'){
      mtext(expression(paste('corrected ',
          beta-'diversity through time')),
            2, line=3, cex=1.5)
    } else {
      mtext(expression(paste(beta-'diversity', ' (corrected)')),
            2, line=3.3, cex=2)
    }
  }
  pdf.f(f, file= file.path(fig.path, sprintf('%s.pdf',
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
                        fig.path <- 'figures/betadisper'){  
  cols <-brewer.pal(6, 'Dark2')[c(6,2,1)]
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
    mtext(expression(paste(beta-'diversity', ' (corrected)')), 2,
          line=3.5 , cex=2)
    ## mtext('Site status',
    ##       side=1, line=4, cex=1.5)
    ## x.labs <- c('Unrestored', 'Maturing', 'Mature')
    ## text(x=1:3 + 0.4,
    ##      y=par('usr')[3] - 0.1,
    ##      labels=x.labs,
    ##      adj=c(1.1,1.1),
    ##      xpd=TRUE,
    ##      cex=1.4)
    axis(1, at= 1:3,
         labels=c('Unrestored', 'Maturing', 'Mature'))
    
  }
  pdf.f(f, file= file.path(fig.path, sprintf('%s.pdf',
             paste('coeff', type,
                   occ,
                   dis.method,
                   sub,
                   sep='_'))),
        width=6, height=6)
}

