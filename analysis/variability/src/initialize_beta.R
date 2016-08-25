library(lme4)
library(vegan)
library(lmerTest)
library(beeswarm)
source('src/misc.R')
source('src/beta.R')
source('src/plotDists.R')

if(type == "pols"){
  if(!binary & alpha){
    occ <- "abund"
    dis.method <- "chao"
    load(file='saved/communities/pols-abund.Rdata')
    load(file='saved/nulls/pols-alpha.Rdata')
  }

  if(!binary & !alpha){
    occ <- "indiv"
    dis.method <- "chao"
    load(file='saved/communities/pols-abund.Rdata')
    load(file='saved/nulls/pols-indiv.Rdata')
  }

  if(binary){
    occ <- "occ"
    dis.method <- "jaccard"
    load(file='saved/communities/pols-abund.Rdata')
    load(file='saved/nulls/pols-occ.Rdata')
  }
  ylabel <- "Species turnover"
}


if(type == "int"){
  if(!binary & alpha){
    occ <- "abund"
    dis.method <- "chao"
    load(file='saved/communities/int-abund.Rdata')
    load(file='saved/nulls/int-alpha.Rdata')
  }

  if(!binary & !alpha){
    occ <- "indiv"
    dis.method <- "chao"
    load(file='saved/communities/int-abund.Rdata')
    load(file='saved/nulls/int-indiv.Rdata')
  }

  if(binary){
    occ <- "occ"
    dis.method <- "jaccard"
    load(file='saved/communities/int-abund.Rdata')
    load(file='saved/nulls/int-occ.Rdata')
  }
  comm.pols <- comm.int
  nulls.pols <- nulls.int
  ylabel <- "Interaction turnover"
}
