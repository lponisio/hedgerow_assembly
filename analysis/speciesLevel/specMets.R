rm(list=ls())
setwd('~/Dropbox/hedgerow_assembly/analysis/speciesLevel')
source('src/initialize.R')
load('../../data/networks/all_networks_years.Rdata')

## **********************************************************
## species importance
## **********************************************************
## linear models
load(file=file.path(save.path, 'specs.Rdata'))

## SiteStatus or ypr
xvar <- "ypr"

## anything outputted by specieslevel
ys <- c("proportional.generality", "d", "degree", "betweenness",
        "closeness.log")

## formulas <-lapply(ys, function(x) {
##   as.formula(paste(x, "~",
##                    paste(paste(xvar, "specialization", sep="*"), 
##                          "(1|Site)",
##                           "(1|GenusSpecies)",
##                          sep="+")))
## })

formulas <-lapply(ys, function(x) {
  as.formula(paste(x, "~",
                   paste(xvar, 
                         "(1|Site)",
                          "(1|GenusSpecies)",
                         sep="+")))
})

mod.pols <- lapply(formulas, function(x){
  lmer(x,
       data=specs[specs$speciesType == "pollinator",])
})

mod.plants <- lapply(formulas, function(x){
  lmer(x,
       data=specs[specs$speciesType == "plant",])
})
names(mod.pols) <- names(mod.plants) <- ys


lapply(mod.plants, summary)
lapply(mod.pols, summary)

save(mod.pols, mod.plants, ys, file=file.path(save.path,
            sprintf('mods/specs_%s.Rdata', xvar)))

## **********************************************************
## degree distributions (abundance distributions)
## **********************************************************

baci.sites <- c("Barger", "Butler", "MullerB", "Sperandio", "Hrdy")
specs <- specs[specs$Site %in% baci.sites,]

layout(matrix(1:6, nrow=2))
cols <- rainbow(length(unique(specs$ypr)))
lapply(unique(specs$Site), function(x){
  print(x)
  this.specs <- specs[specs$Site == x, c("degree", "ypr")]
  plot(NA, ylim=c(0,0.8), xlim=c(0,25),
       ylab="Frequency",
       xlab="Abundance",
       main=x)
  for(i in 1:length(unique(this.specs$ypr))){
    this.ypr <- unique(this.specs$ypr)[i]
    print(this.ypr)
    points(density(this.specs$degree[this.specs$ypr == this.ypr]),
           col=cols[i], type="l", lwd=2)
  }
})

plot(NA, ylim=c(0,1), xlim=c(0,1), xaxt="n", yaxt="n", ylab="", xlab="")
legend("center", col=cols, lwd="2",
       legend=sort(unique(specs$ypr)),
       bty="n")
