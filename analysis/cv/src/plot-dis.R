## dissimilarity of plants, pol, int

plot.dis <- function(c.pol, lm.pol, c.plant, lm.plant, c.int,
                     lm.int, sitetype = "hedgerow", dis.method,
                     trait="d",
                     method= "time"){

  f.curve <- function(x, intercept, slope){
    y <- intercept + x * slope
    return(y)
  }

  if(sitetype=="cv"){
    print(sitetype)
    plot.3.coeff <- function(){
      cols <- hsv(c(0.1, 0.35, 0.75), alpha=0.2)
      cols.lines <- hsv(c(0.1, 0.35, 0.75), alpha=1)
      pc <- function(c.ob, lm.ob,...){

        plot(c.ob[c.ob[,"status"] == "mature","cv"]~
             c.ob[c.ob[,"status"] == "mature", "traits"],
             pch =16, col=cols[1], ylim=range(c.ob[,"cv"],
                                     na.rm=TRUE),
             xlim=range(c.ob[, "traits"], na.rm=TRUE),
             ylab="Coefficient of variation",...)

        curve(f.curve(x=x, intercept =
                      (coef(summary(lm.ob))[1,1] +
                       coef(summary(lm.ob))[2,1]),
                      slope=(coef(summary(lm.ob))[4,1] +
                             coef(summary(lm.ob))[5,1])),
              lwd=1.5, col=cols.lines[1], add=TRUE)
        
        points(c.ob[c.ob[,"status"] == "maturing","cv"]~
               c.ob[c.ob[,"status"] == "maturing", "traits"],
               pch=16, col=cols[2])
        curve(f.curve(x=x, intercept =
                      (coef(summary(lm.ob))[1,1] +
                       coef(summary(lm.ob))[3,1]),
                      slope=(coef(summary(lm.ob))[4,1] +
                             coef(summary(lm.ob))[6,1])),
              lwd=1.5, col=cols.lines[2], add=TRUE)

        points(c.ob[c.ob[,"status"] == "control","cv"]~
               c.ob[c.ob[,"status"] == "control", "traits"], pch=16,
               col=cols[3])
        curve(f.curve(x=x, intercept =
                      coef(summary(lm.ob))[1,1],
                      slope=coef(summary(lm.ob))[4,1]),
              lwd=1.5, col=cols.lines[3], add=TRUE)
      }

      f <- function(){
        par(mfrow=c(1,1), oma=c(2.5, 2.5, 0.25, 1),
            mar=c(5, 5, 2, 0.5), cex.axis=1, cex.lab=1)
        pc(c.pol, lm.pol, xlab=trait)
        legend("bottomright", legend= c("Mature", "Maturing",
                                "Unrestored"), col=cols.lines,
               pch=16, bty="n", cex=1.5)
      }
      path <- 'figures' 
      pdf.f(f, file= file.path(path, sprintf("%s.pdf",
                 paste(trait, method, sep="_"))),
            width=6, height=6)      
    }
    plot.3.coeff()
  } else if(sitetype=="box"){
    plot.box <- function(){
      pc <- function(dats){
        cols <- hsv(c(0.1, 0.35, 0.75), alpha=1)
        dats <- dats[!is.na(dats[,"traits"]),]
        cats <-    length(unique(dats[,"status"]))*
          length(unique(dats[,"traits"]))
        boxplot(dats[,"cv"] ~ dats[,"status"]*dats[,"traits"],
                col=cols, names=rep("", cats),
                ylim= range(dats[, "cv"], na.rm=TRUE) + c(0, 100),
                ylab=" Coefficient of variation")
        if(cats == 6) {
          mtext(at=c(2, 5), side= 1,
                text= unique(dats[, "traits"]), line=3, cex=2)
        }
        if(cats == 9) {
          mtext(at=c(2, 5, 8), side= 1,
                text= unique(dats[, "traits"]), line=3, cex=2)
        }
        legend("topright", legend= c("Mature", "Maturing",
                             "Unrestored"), col=cols,
               pch=16, bty="n", cex=1.5)
      }
      f <- function(){
        par(mfrow=c(1,1), oma=c(2.5, 2.5, 0.25, 1),
            mar=c(5, 5, 2, 0.5), cex.axis=1.5, cex.lab=1.5)
        pc(c.pol)
      }
      path <- 'figures' 
      pdf.f(f, file= file.path(path, sprintf("%s.pdf",
                 paste(trait, method, sitetype, sep="_"))),
            width=15, height=10)
    }
    plot.box()
  }
}
