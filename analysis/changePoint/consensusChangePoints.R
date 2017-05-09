source('src/extractingOutput.R')
save.dir <- "saved/runs"

run.change.points <- list.files(save.dir)

all.runs <- lapply(run.change.points, function(x){
    read.csv(file.path(save.dir, x))
})

removeDup <- function(y){
    y <- y[!duplicated(cbind(y$sites, y$cp)),]
}

all.runs <- lapply(all.runs, removeDup)

all.runs <- do.call(rbind, all.runs)

sum.runs <- aggregate(list(num.runs = all.runs$sites),
                      list(sites=all.runs$sites,
                           cp=all.runs$cp),
             function(x) length(x)/length(run.change.points))

write.csv(sum.runs, file="saved/consensusChangePoints.csv",
          row.names=FALSE)


## drop change points that do not occur in 95% of runs
sum.runs <- sum.runs[sum.runs$num.runs >= 0.9,]

makeConsensusTable(changing.points=sum.runs)
