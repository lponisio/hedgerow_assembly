## regresses coefficent of variation against traits

cv.trait <- function(spec.dat,
                     byType,
                     trait,
                     cont=TRUE,
                     method,
                     status.order=c("control", "maturing", "mature")){
  byStatus <- split(byType, byType$status)
  bySite <- lapply(byStatus, function(x) {split(x, x$site)})
  bySite <- unlist(bySite, recursive=FALSE)
  prep.cv <- lapply(bySite, function(y) {
    samp2site.spp(y[,"date"], y[,"GenSp"], y[,"Abund"])
  })
  coeff.cv <- lapply(prep.cv, function(x){apply(x, 2, cv)})
  dats <- data.frame(cv=unlist(coeff.cv))
  dats$status <-  gsub('\\..*', '', rownames(dats))
  dats$status <- factor(dats$status, levels=status.order)
  dats$sp <- unlist(lapply(coeff.cv, names))
  dats$site <-  sapply(strsplit(rownames(dats), "\\."),
                       function(x) x[2])
  rownames(dats) <- NULL
  if(cont){
    dats$traits.ns <- spec.dat[,trait][match(dats$sp,
                                             spec.dat$GenusSpecies)]
    dats$traits <- scale(dats$traits.ns)
  } else{
    dats$traits <- spec.dat[,trait][match(dats$sp,
                                          spec.dat$GenusSpecies)]
  }
  lm.cv <- lmer(cv ~ status*traits + (1|site) + (1|sp),
                data=dats, REML=FALSE)
  return(list(data=dats, lm=lm.cv))
}


## previous plotting code. Takes xlabel argument
## if(cont){
##   plot.dis(dats, lm.cv, sitetype="cv", trait=trait, method= method)
##   predict.int(lm.cv, xlabel=xlabel, dats=dats$traits, method= method)
## } else{
##   plot.dis(dats, lm.cv, sitetype="box", trait=trait, method= method)
## }
