library(igraph)
library(bipartite)
library(lme4)
source('../../dataPrep/src/prepNets.R')
source('src/CalcMetrics.R')
source('src/misc.R')

traits <- read.csv("../../data/traits.csv")
load('../../data/networks/allSpecimens.Rdata')
