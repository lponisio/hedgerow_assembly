library(igraph)
library(bipartite)
library(lme4)
library(lmerTest)
source('../../dataPrep/src/prepNets.R')
source('src/misc.R')

traits <- read.csv("../../data/traits.csv")
load('../../data/networks/allSpecimens.Rdata')

save.path <- 'saved'
