library(igraph)
library(bipartite)
library(lme4)
library(lmerTest)
library(nlme)
library(sme)

source('../../dataPrep/src/prepNets.R')
source('src/CalcMetrics.R')
source('src/misc.R')
source("src/resilience.R")
source('src/diffs.R')

traits <- read.csv("../../data/traits.csv")
load('../../data/networks/allSpecimens.Rdata')

save.path <- 'saved'



## splines
calcSpine <- function(x, calc.log=FALSE){
    if(!calc.log){
        spline.mod <- sme(baci[,x], tme=baci$ypr, ind=baci$Site)
    } else{
        spline.mod <- sme(log(baci[,x]), tme=baci$ypr, ind=baci$Site)
    }
    print(spline.mod[c("info","smoothingParameters")])
    return(spline.mod)
}
