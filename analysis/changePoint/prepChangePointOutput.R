rm(list=ls())
library(reshape)
setwd('~/Dropbox/hedgerow_assembly/analysis/changePoint/cptPeel')
source('../src/extractingOutput.R', chdir = TRUE)

samples <- read.csv("../../../data/samples.csv")
results <- read.table('results_4.txt', sep=' ')
logs <- read.table('LogLs_4.txt', sep=' ')

chpts <- makeChangepointData(results=results, logs=logs, samples=samples,
                             value=0.949, w=4)
## change the value argument to the "p value" to be considered, like 0.949

makeConsensusTable(changing.points=chpts)
