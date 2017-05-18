
calcBetaStatus <- function(comm,
                        status,
                        dis.method,
                        nulls,
                        sub= "pol",
                        occ= FALSE,
                        years,
                        zscore=FALSE){ ## calculate zscores?
    ## computes dispersion of community matrices, returns output of
    ## vegan function betadisper

    ## create community dissimilarity matrix
    comm.dis <-  as.matrix(vegdist(comm, method= dis.method,
                                   diag= TRUE, binary= occ))
    ## null dissimilarity matrices
    null.dis <- lapply(nulls, function(x) {
        as.matrix(vegdist(x, method= dis.method, diag=TRUE))
    })

    null.dis[[length(nulls) + 1]] <- comm.dis
    arr <- array(unlist(null.dis), c(dim(comm.dis)[1],
                                     dim(comm.dis)[2],
                                     length(nulls) + 1))

    ## standardize dissimilarities
    if(!zscore){
        less.than  <-   apply(arr, 1:2, function(x){
            sum(x[length(null.dis)] > x)
        })
        equal.2  <-   apply(arr, 1:2, function(x){
            sum(x[length(null.dis)] == x)
        })
        cor.dis <- as.dist((less.than + 0.5*equal.2)/
                           length(null.dis), diag= TRUE)
    }else{
        cor.dis  <-  (comm.dis -
                      apply(arr , 1:2 , mean))/
            (apply(arr , 1:2 , sd) + 10^-10)
        cor.dis <-  as.dist(((cor.dis - min(cor.dis))/diff(range(cor.dis))),
                            diag= TRUE)
    }

    ## run model
    mod.beta <- try(betadisper(cor.dis, status,
                               type="centroid"), silent=TRUE)
    if(inherits(mod.beta, "try-error")) browser()
    return(mod.beta)
}

