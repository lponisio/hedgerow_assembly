library(reshape)
library(tidyr)

makeChangepointData <- function(results, logs,
                                save.path="../saved/"){
    ## makes a useful data set with the change points at each site
    ## from the output of the change point detection analysis
    ## split the year and creates an extra table with year 1 and
    ## period

    results <- results[sort(results$V1),]
    logs <- logs[sort(logs$V1),]
    results$V1 <- as.character(results$V1)
    results$V2 <- as.character(results$V2)
    years1  <-  as.numeric(sapply(strsplit(results$V1, '_'),
                                  function(x) x[2]))
    years2  <-  as.numeric(sapply(strsplit(results$V2, '_'),
                                  function(x) x[2]))
    anos <- min(years1):max(years2)
    anos.unique <- unique(c(years1, years2))
    anos <- anos[anos %in% anos.unique]

    years <- as.data.frame(matrix(NA, (length(anos)-1),2))
    for(n in 1:(length(anos)-1)){
        years[n,2] <- paste(anos[n], anos[n+1], sep='-')
    }
    years[,1] <- anos[-length(anos)]
    colnames(years) <- c('first', 'period')

    ## creates other variables for results
    sites  <- sapply(strsplit(results$V1, '_'), function(x) x[1])
    max <- apply(results[,c(3:5)], 1, max)
    numberOfOnes <- apply(results[,c(3:5)], 1,
                          function(x) length(which(x==1)))
    max.index <- rep(NA, dim(results)[1])
    results <- cbind(results,  sites, years1, max, numberOfOnes, max.index)

    ## checking if there is more than 1 value, and if yes, returns the
    ## index with the maxlog
    for(row in 1:dim(results)[1]){
        if(results$numberOfOnes[row] > 1){
            results$max.index[row] <- which.max(logs[row,c(3:5)])
        }
        else {## return the position of the max
            results$max.index[row] <- which(results[row,c(3:5)] ==
                                            results$max[row])[1]
        }
    }
    sigs <- results[results$max > 0.949,]
    cp <- rep(NA, dim(sigs)[1])
    value <- rep(NA, dim(sigs)[1])
    sigs <- cbind(sigs, cp, value)

    ## cp is always the year + next year
    for(i in 1:dim(sigs)[1]){
        y1 <- which(years$first %in% sigs$years1[i] == TRUE)
        sigs$cp[i] <- years$period[y1+(sigs$max.index[i]-1)]
        sigs$value[i] <- sigs[i,c(3:5)][sigs$max.index[i]]
    }
    sigs$value <- unlist(sigs$value)
    changing.points <- sigs[,c(6,11,12)]

    write.csv(changing.points,
              file=file.path(save.path, 'changing_points.csv'),
              row.names=FALSE)
    return(changing.points)
}

###

makeConsensusTable <- function(changing.points,
                               pairs.path='baci',
                               save.path='../saved/'){
    ## makes a table of consesus networks between each of the changing
    ## points
    changing.points$sites <- as.character(changing.points$sites)
    by.site <- split(changing.points, changing.points$sites)

    files <-
        list.files(pairs.path,
                   pattern='pairs')

    files <- files[sapply(strsplit(files, '_'), function(x) x[1])
                   %in% changing.points$sites]

    split.ch.yrs <- lapply(by.site, function(x){
        separate(x, cp, into=c('min','max'))
    })

    clusters <- vector(mode='list', length(split.ch.yrs))

    for(j in 1:length(split.ch.yrs)){
        this.site <- split.ch.yrs[[j]]
        this.site <- this.site[!duplicated(this.site$min),]
        this.site <- this.site[order(this.site$min),]
        clusters[[j]] <- list()
        for(i in 1:nrow(this.site)){
            clusters[[j]][[i]] <- try(2006:this.site$min[i],
                                      silent=TRUE)
            if(inherits(clusters[[j]][[i]], "try-error")) browser()
            if(i >= 2){
                clusters[[j]][[i]] <-
                    clusters[[j]][[i]][!clusters[[j]][[i]] %in%
                                       clusters[[j]][[i - 1]]]
            }
            if(i == nrow(this.site)){
                clusters[[j]][[i + 1]] <- this.site$min[i]:2014
            }
        }
    }
    names(clusters) <- names(split.ch.yrs)
    clusters.pairs <- clusters

    for(i in 1:length(clusters.pairs)){
        for(j in 1:length(clusters.pairs[[i]])){
            clusters.pairs[[i]][[j]] <- paste(names(split.ch.yrs)[i], '_',
                                              clusters.pairs[[i]][[j]],
                                              '.pairs', sep='')
        }
    }

    clusters.pairs <- lapply(clusters.pairs, function(x){
        lapply(x, function(y){
            y <- y[y %in% files]
        })
    })
    pasteComma <- function(...){
        paste(..., sep=', ')
    }

    ## create table for consensus function
    prep.table <- unlist(lapply(clusters.pairs,
                                function(x) do.call(pasteComma, x)))

    ## create table to convert consensus trees to something useable
    last.yr <- rapply(clusters, max, how="unlist")
    names.yr <- substr(names(last.yr), 1, nchar(names(last.yr))-1)

    out <- character(length(last.yr))
    for(i in 1:length(last.yr)){
        out[i] <- sprintf("baci/%s_%s_joint_s1-consensus.tree",
                          names.yr[i], last.yr[i])
    }

    write.table(prep.table, file=file.path(save.path, 'consensus.txt'),
                col.names=FALSE,
                row.names=FALSE)

    write.table(out, file=file.path(save.path, 'lastyr_consensus.txt'),
                col.names=FALSE,
                row.names=FALSE)


}
