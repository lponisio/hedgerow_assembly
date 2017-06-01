library(reshape)
library(tidyr)

makeChangepointData <- function(results, logs, value, samples,
                                save.path="../saved/runs",
                                w=4,
                                file.name="changing_points.csv"){
    ## makes a useful data set with the change points at each site

    ## from the output of the change point detection analysis
    ## split the year and creates an extra table with year 1 and
    ## period

    ##The change point analysis (cp) results gives a table with five
    ##columns: the first two represents the first and last years, and the
    ##other three the probabilities that a cp occured, for each
    ##interval. For example, if the first two collumns are 2006 and 2009,
    ##the three values represent the cp values for intervals 2006-2007,
    ##2007-2008, 2008-2009. A cp occures when there is a value=1.00.  When
    ##you have more than four years, there is more than one line for that
    ##site, and thus, it can have more than one value per period. When that
    ##happens, you have to check in the logs file to see which one has the
    ##greater likelihood.

    ##The code below is to identify the years that were sampled in order to
    ##create the channging points table from a samples table that has the
    ##sites in the rows and the number of times each year was sampled in
    ##the columns, for all years. The output is a list in which each
    ##dimension is a sample site, and there is a dataframe with the year
    ##the sample started and the period that follows.
    print(file.name)
    years <- as.numeric(sapply(strsplit(colnames(samples), 'X'),
                               function(x) x[2]))
    ## create a vector of the periods
    periods <- apply(samples, 1, function(x) which(x!=0))
    names(periods) <- samples[,1]
    periods2 <- list()
    ## for each site checks if there is more than one year sampled.
    for(p in 1: length(periods)){
        anos <- years[periods[[p]]]
        anos <- anos[!is.na(anos)]
        if(length(anos)!=1){
            years2 <- as.data.frame(matrix(NA, (length(anos)-1),2))
            for(n in 1:(length(anos)-1)){
                years2[n,2] <- paste(anos[n], anos[n+1], sep='-')
            }
            years2[,1] <- anos[-length(anos)]
            colnames(years2) <- c('first', 'period')
            periods2[[p]]<-years2
        }else{
            periods2[p]<- anos}
    }

    names(periods2) <- samples[,1]


    results <- results[sort(results$V1),]
    logs <- logs[sort(logs$V1),]
    ## checking to make sure it is the same order
    if("FALSE" %in% (results[,1:2]==logs[,1:2])){print("Stop! Names in files don't match")}

    ## The code below is just to create the periods between years
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
    ## Identifyes the sites' name
    sites  <- sapply(strsplit(results$V1, "_"), function(x) x[1])
    ## Return the maximum values per line

    maxs <- apply(as.matrix(results[,c(3:w + 1)]), 1, max)
    ## Return the number of changing points values per line
    numberOfcps <- apply(as.matrix(results[,c(3:w+1)]), 1,
                         function(x) length(which(x > value)))
    ##creating a dummy column to identify the period, which will be
    ##the year that specific period (line) started, plus the column
    ##number.
    max.index <- rep(NA, dim(results)[1])
    results <- cbind(results,  sites, years1, maxs, numberOfcps, max.index)

    ##Separating only the lines that have values of 1.00 (changing
    ##points idetified) for the results and the logs
    sigs <- results[results$maxs > value,]
    logs.sigs <- logs[results$maxs > value,]
    ## checking if there is more than 1 value, and if yes, returns the
    ## index with the maxlog (log table)

    for(row in 1:dim(sigs)[1]){
        if(sigs$numberOfcps[row] > 1){
            sigs$max.index[row] <- which.max(logs.sigs[row,c(3:w+1)])
        }
        else {## return the position of the max
            sigs$max.index[row] <- which(sigs[row,c(3:w+1)] ==
                                         sigs$maxs[row])[1]
        }
    }
    ##aadding the period
    period <- rep(NA, dim(sigs)[1])
    sigs <- cbind(sigs, period)

    ## cp occurs between two years, indicated by the first year of the
    ## line plus the column where the 1.0 is found
    for(i in 1:dim(sigs)[1]){
        local <- which(sigs$sites[i] == names(periods2))
        y1 <- which(periods2[[local]][,1] == sigs$years1[i])
        sigs$period[i] <- periods2[[local]][y1+(sigs$max.index[i]-1),2]
    }

    changing.points <- sigs[,c("sites", "maxs", "period")]
    colnames(changing.points)<-c("sites", "value", "cp")
    write.csv(changing.points,
              file=file.path(save.path, file.name),
              row.names=FALSE)
    return(changing.points)
}

###

makeConsensusTable <- function(changing.points,
                               pairs.path='cptPeel/baci',
                               save.path='saved/'){
    ## makes a table of consensus networks between each of the changing
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
            clusters[[j]][[i]] <- 2006:this.site$min[i]
            if(i >= 2){
                clusters[[j]][[i]] <-
                    clusters[[j]][[i]][!clusters[[j]][[i]] %in%
                                       clusters[[j]][[i - 1]]]
            }
            if(i == nrow(this.site)){
                clusters[[j]][[i + 1]] <- (as.numeric(this.site$min[i])+1):2014
            }
        }
    }

    names(clusters) <- names(split.ch.yrs)

    for(i in 1:length(clusters)){
        for(j in 1:length(clusters[[i]])){
            clusters[[i]][[j]] <- paste(names(split.ch.yrs)[i], '_',
                                        clusters[[i]][[j]], '.pairs', sep='')
        }
    }

    obs.clusters <- lapply(clusters, function(x){
        lapply(x, function(y){
            y <- y[y %in% files]
        })
    })

    pasteComma <- function(...){
        paste(..., sep=', ')
    }

    prep.table <- rapply(obs.clusters, function(x){
        paste(x, collapse=", ")
    },
    how="unlist")
    names(prep.table) <- NULL

    last.year <- rapply(obs.clusters, function(x){
        x[length(x)]
    }, how="unlist")
    names(last.year) <- NULL
    last.year <- gsub(".pairs", "_joint_s1-consensus.tree", last.year)
    last.year <- paste0("baci/", last.year)

    write.table(last.year, file=file.path(save.path, 'lastyr_consensus.txt'),
                col.names=FALSE, row.names=FALSE)
    write.table(prep.table, file=file.path(save.path, 'consensus.txt'),
                col.names=FALSE, row.names=FALSE)


}
