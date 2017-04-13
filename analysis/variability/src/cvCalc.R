library(lme4)
library(lqmm)

## regresses coefficent of variation against traits

corCv <- function(x){
    cv(x)*(1 + (1/4*length(x)))
}

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
    ## cv.function can be corCv, cv, or sd
    ## standard cv devides cv by 100 to avoid too large of numbers for the
    ## lm
    byStatus <- split(byType, byType$SiteStatus)
    bySite <- lapply(byStatus, function(x) {split(x, x$Site)})
    bySite <- unlist(bySite, recursive=FALSE)
    prep.cv <- lapply(bySite, function(y) {
        samp2site.spp(y[, time.col], y[,"GenusSpecies"], y[, abund.col])
    })
    ## to avoid zeros in calculation of cv
    if(zero2na){
        prep.cv <- lapply(prep.cv, function(x){
            x[x == 0] <- NA
            return(x)
        })
    }
    coeff.cv <- lapply(prep.cv, function(x){apply(x, 2, cv.function, ...)})
    dats <- data.frame(cv=unlist(coeff.cv))
    if(standard.cv){
        dats$cv[!is.na(dats$cv)] <- log(dats$cv[!is.na(dats$cv)])
    }
    dats$SiteStatus <-  gsub('\\..*', '', rownames(dats))
    dats$SiteStatus <- factor(dats$SiteStatus, levels=status.order)
    dats$GenusSpecies <- unlist(lapply(coeff.cv, names))
    dats$Site <-  sapply(strsplit(rownames(dats), "\\."),
                         function(x) x[2])
    rownames(dats) <- NULL
    if(cont){
        dats$traits.ns <- spec.dat[,trait][match(dats$GenusSpecies,
                                                 spec.dat[,
                                                          species.type])]
        ## trying out not scaling... does it break everything?
        dats$traits <- dats$traits.ns
    } else{
        dats$traits <- spec.dat[,trait][match(dats$GenusSpecies,
                                              spec.dat[, species.type])]
    }
    lm.dats <- dats[!is.na(dats$cv) & !is.na(dats$traits),]
    ## reviwers wanted a quantile model
    ## can only include one random effect in the quantile mixed
    ## effects model
    lm.cv.quant <- lqmm(fixed=cv ~ traits, random=~1, group=
                     GenusSpecies,
                        data=lm.dats,
                     control=list(LP_max_iter=10^3))
    sum.boot.quant <- summary.boot.lqmm(boot(lm.cv.quant))
    lm.cv.nss <- lmer(cv ~ traits + (1|Site) + (1|GenusSpecies),
                      data=lm.dats)
    return(list(data=dats,
                lm.nss=lm.cv.nss,
                lm.cv.quant=lm.cv.quant,
                sum.quant=sum.boot.quant))
}
