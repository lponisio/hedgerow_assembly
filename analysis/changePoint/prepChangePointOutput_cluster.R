setwd('cptPeel')

setwd('~/Dropbox/hedgerow_assembly/analysis/changePoint/cptPeel')
source('../src/extractingOutput.R')


w <- 4

samples <- read.csv("../../../data/samples.csv")
results <- list.files(pattern=sprintf('results_%s', w))
logs <- list.files(pattern=sprintf('LogLs_%s', w))

runs.res <- gsub(pattern=sprintf('results_%s', w), replacement="",
                 x=results)
runs.logs <- gsub(pattern=sprintf('LogLs_%s', w), replacement="",
                  x=logs)

## order the log and result files in case runs were broken on the
## cluster and log files were not produced
logs <- logs[match(runs.res, runs.logs)]
results <- results[match(runs.logs, runs.res)]

logs <- logs[!is.na(logs)]
logs <- logs[!is.na(results)]
results <- results[!is.na(logs)]
results <- results[!is.na(results)]

load.results <- lapply(results, read.table, sep=' ')
load.logs <- lapply(results, read.table, sep=' ')

chpts <- mapply(function(a, b, c)
    makeChangepointData(results= a,
           logs= b,
           file.name=c,
           value=0.949,
           samples=samples,
           w=w),
    a=load.results,
    b=load.logs,
    c=results,
    SIMPLIFY=FALSE)
