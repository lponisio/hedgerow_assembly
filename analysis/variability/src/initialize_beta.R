library(lme4)
library(vegan)
library(lmerTest)
library(beeswarm)
library(nlme)
source('src/misc.R')
source('src/beta.R')
source('src/plotDists.R')

if(!binary & alpha){
  occ <- "abund"
  dis.method <- "chao"
  load(file=file.path('saved/communities',
         sprintf('%s-abund.Rdata', type)))
  load(file=file.path('saved/nulls',
         sprintf('%s-alpha.Rdata', type)))
}

if(!binary & !alpha){
  occ <- "indiv"
  dis.method <- "chao"
  load(file=file.path('saved/communities',
         sprintf('%s-abund.Rdata', type)))
  load(file=file.path('saved/nulls',
         sprintf('%s-indiv.Rdata', type)))
}

if(binary){
  occ <- "occ"
  dis.method <- "jaccard"
  load(file=file.path('saved/communities',
         sprintf('%s-abund.Rdata', type)))
  load(file=file.path('saved/nulls',
         sprintf('%s-occ.Rdata', type)))
}

if(type=="pols"){
  ylabel <- "Pollinator species turnover"
}
if(type=="ints"){
  ylabel <- "Interaction turnover"
}
if(type=="plants"){
  ylabel <- "Plant species turnover"
}
