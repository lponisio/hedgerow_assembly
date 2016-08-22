## regresses coefficent of variation against traits

corCv <- function(x){
  cv(x)*(1 + (1/4*length(x)))
}

## cv.function can be corCv, cv, or sd
## standard cv devides cv by 100 to avoid too large of numbers for the
## lm
cv.trait <- function(spec.dat,
                     byType,
                     trait,
                     cont=TRUE,
                     method,
                     time.col,
                     abund.col,
                     status.order=c("control", "maturing", "mature"),
                     cv.function=corCv,
                     zero2na =FALSE,
                     standard.cv=TRUE,
                     species.type="GenusSpecies",...){
  byStatus <- split(byType, byType$SiteStatus)
  bySite <- lapply(byStatus, function(x) {split(x, x$Site)})
  bySite <- unlist(bySite, recursive=FALSE)
  prep.cv <- lapply(bySite, function(y) {
    samp2site.spp(y[, time.col], y[,"GenusSpecies"], y[, abund.col])
  })
  if(zero2na){
    prep.cv <- lapply(prep.cv, function(x){
      x[x == 0] <- NA
      return(x)
    })
  }
  coeff.cv <- lapply(prep.cv, function(x){apply(x, 2, cv.function, ...)})
  dats <- data.frame(cv=unlist(coeff.cv))
  if(standard.cv){
    dats$cv <- dats$cv/100
  }
  dats$SiteStatus <-  gsub('\\..*', '', rownames(dats))
  dats$SiteStatus <- factor(dats$SiteStatus, levels=status.order)
  dats$GenusSpecies <- unlist(lapply(coeff.cv, names))
  dats$Site <-  sapply(strsplit(rownames(dats), "\\."),
                       function(x) x[2])
  rownames(dats) <- NULL
  if(cont){
    dats$traits.ns <- spec.dat[,trait][match(dats$GenusSpecies,
                                             spec.dat[, species.type])]
    dats$traits <- scale(dats$traits.ns)
  } else{
    dats$traits <- spec.dat[,trait][match(dats$GenusSpecies,
                                          spec.dat[, species.type])]
  }
  lm.cv <- lmer(cv ~ SiteStatus*traits + (1|Site) + (1|GenusSpecies),
                data=dats[!is.na(dats$cv),])
  lm.cv.nss <- lmer(cv ~ traits + (1|Site) + (1|GenusSpecies),
                data=dats[!is.na(dats$cv),])
  return(list(data=dats, lm=lm.cv, lm.nss=lm.cv.nss))
}


## previous plotting code. Takes xlabel argument
## if(cont){
##   plot.dis(dats, lm.cv, sitetype="cv", trait=trait, method= method)
##   predict.int(lm.cv, xlabel=xlabel, dats=dats$traits, method= method)
## } else{
##   plot.dis(dats, lm.cv, sitetype="box", trait=trait, method= method)
## }
