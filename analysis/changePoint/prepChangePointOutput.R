rm(list=ls())
library(reshape)
setwd('~/Dropbox/hedgerow_assembly/analysis/changePoint/cptPeel')
source('../src/extractingOutput.R', chdir = TRUE)

results <- read.table('results_4.txt', sep=' ')
logs <- read.table('LogLs_4.txt', sep=' ')

chpts <- makeChangepointData(results=results, logs=logs)

makeConsensusTable(changing.points=chpts)
