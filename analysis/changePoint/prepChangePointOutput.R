setwd('changePoint/cptPeel')
## setwd('~/Dropbox/hedgerow_assembly/analysis/changePoint/cptPeel')
source('../src/extractingOutput.R')

args <- commandArgs(trailingOnly=TRUE)
print(args)

w <- 4

samples <- read.csv("../../../data/samples.csv")
results <- read.table(sprintf('results_%s.txt', w), sep=' ')
logs <- read.table(sprintf('LogLs_%s.txt', w), sep=' ')

chpts <- makeChangepointData(results=results, logs=logs, samples=samples,
                             value=0.949, w=4, file.name=args[1])
## change the value argument to the "p value" to be considered, like 0.949

makeConsensusTable(changing.points=chpts)
