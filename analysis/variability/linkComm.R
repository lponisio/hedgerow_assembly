function (network, hcmethod = "average", use.all.edges = FALSE,
    edglim = 10^4, directed = FALSE, dirweight = 0.5, bipartite = FALSE,
    dist = NULL, plot = TRUE, check.duplicates = TRUE, removetrivial = TRUE,
    verbose = TRUE)
{
    if (is.character(network) && !is.matrix(network)) {
        if (file.access(network) == -1) {
            stop(cat("\nfile not found: \"", network, "\"\n",
                sep = ""))
        }
        else {
            network <- read.table(file = network, header = FALSE)
        }
    }
    x <- network
    rm(network)
    if (ncol(x) == 3) {
        wt <- as.numeric(as.character(x[, 3]))
        if (length(which(is.na(wt) == TRUE)) > 0) {
            stop("\nedge weights must be numerical values\n")
        }
        x <- cbind(as.character(x[, 1]), as.character(x[, 2]))
    }
    else if (ncol(x) == 2) {
        x <- cbind(as.character(x[, 1]), as.character(x[, 2]))
        wt <- NULL
    }
    else {
        stop("\ninput data must be an edge list with 2 or 3 columns\n")
    }
    if (check.duplicates) {
        dret <- edge.duplicates(x, verbose = verbose)
        x <- dret$edges
        if (!is.null(wt)) {
            if (length(dret$inds) > 0) {
                wt <- wt[-dret$inds]
            }
        }
        rm(dret)
    }
    el <- x
    rm(x)
    len <- nrow(el)
    nnodes <- length(unique(c(as.character(el[, 1]), as.character(el[,
        2]))))
    intel <- integer.edgelist(el)
    edges <- intel$edges
    node.names <- names(intel$nodes)
    numnodes <- length(node.names)
    if (bipartite) {
        big <- graph.edgelist(as.matrix(el), directed = directed)
        bip.test <- bipartite.mapping(big)
        if (!bip.test$res) {
            stop("\nnetwork is not bi-partite; change bipartite argument to FALSE\n")
        }
        bip <- rep(1, length(bip.test$type))
        bip[which(bip.test$type == FALSE)] <- 0
        names(bip) <- V(big)$name
        bip <- bip[match(node.names, names(bip))]
        rm(big, bip.test)
    }
    else {
        bip <- 0
    }
    rm(intel)
    if (len <= edglim) {
        disk <- FALSE
        if (is.null(dist)) {
            emptyvec <- rep(1, (len * (len - 1))/2)
            if (!is.null(wt)) {
                weighted <- TRUE
            }
            else {
                wt <- 0
                weighted <- FALSE
            }
            if (!use.all.edges) {
                dissvec <- .C("getEdgeSimilarities", as.integer(edges[,
                  1]), as.integer(edges[, 2]), as.integer(len),
                  rowlen = integer(1), weights = as.double(wt),
                  as.logical(directed), as.double(dirweight),
                  as.logical(weighted), as.logical(disk), dissvec = as.double(emptyvec),
                  as.logical(bipartite), as.logical(verbose))$dissvec
            }
            else {
                dissvec <- .C("getEdgeSimilarities_all", as.integer(edges[,
                  1]), as.integer(edges[, 2]), as.integer(len),
                  as.integer(numnodes), rowlen = integer(1),
                  weights = as.double(wt), as.logical(FALSE),
                  as.double(dirweight), as.logical(weighted),
                  as.logical(disk), dissvec = as.double(emptyvec),
                  as.logical(bipartite), as.logical(verbose))$dissvec
            }
            distmatrix <- matrix(1, len, len)
            distmatrix[lower.tri(distmatrix)] <- dissvec
            colnames(distmatrix) <- 1:len
            rownames(distmatrix) <- 1:len
            distobj <- as.dist(distmatrix)
            rm(distmatrix)
        }
        else {
            if (class(dist) != "dist") {
                stop("\ndistance matrix must be of class \"dist\" (see ?as.dist)\n")
            }
            else if (attr(dist, which = "Size") != len) {
                stop("\ndistance matrix size must equal the number of edges in the input network\n")
            }
            else if (length(dist) != (len * (len - 1))/2) {
                stop("\ndistance matrix must be the lower triangular matrix of a square matrix\n")
            }
            distobj <- dist
        }
        if (verbose) {
            cat("\n   Hierarchical clustering of edges...")
        }
        hcedges <- hclust(distobj, method = hcmethod)
        hcedges$order <- rev(hcedges$order)
        rm(distobj)
        if (verbose) {
            cat("\n")
        }
    }
    else {
        disk <- TRUE
        if (!is.null(wt)) {
            weighted <- TRUE
        }
        else {
            wt <- 0
            weighted <- FALSE
        }
        if (!use.all.edges) {
            rowlen <- .C("getEdgeSimilarities", as.integer(edges[,
                1]), as.integer(edges[, 2]), as.integer(len),
                rowlen = integer(len - 1), weights = as.double(wt),
                as.logical(directed), as.double(dirweight), as.logical(weighted),
                as.logical(disk), dissvec = double(1), as.logical(bipartite),
                as.logical(verbose))$rowlen
        }
        else {
            rowlen <- .C("getEdgeSimilarities_all", as.integer(edges[,
                1]), as.integer(edges[, 2]), as.integer(len),
                as.integer(numnodes), rowlen = integer(len -
                  1), weights = as.double(wt), as.logical(FALSE),
                as.double(dirweight), as.logical(weighted), as.logical(disk),
                dissvec = double(1), as.logical(bipartite), as.logical(verbose))$rowlen
        }
        if (verbose) {
            cat("\n")
        }
        hcobj <- .C("hclustLinkComm", as.integer(len), as.integer(rowlen),
            heights = single(len - 1), hca = integer(len - 1),
            hcb = integer(len - 1), as.logical(verbose))
        if (verbose) {
            cat("\n")
        }
        hcedges <- list()
        hcedges$merge <- cbind(hcobj$hca, hcobj$hcb)
        hcedges$height <- hcobj$heights
        hcedges$order <- .C("hclustPlotOrder", as.integer(len),
            as.integer(hcobj$hca), as.integer(hcobj$hcb), order = integer(len))$order
        hcedges$order <- rev(hcedges$order)
        hcedges$method <- "single"
        class(hcedges) <- "hclust"
    }
    hh <- unique(round(hcedges$height, digits = 5))
    countClusters <- function(x, ht) {
        return(length(which(ht == x)))
    }
    clusnums <- sapply(hh, countClusters, ht = round(hcedges$height,
        digits = 5))
    numcl <- length(clusnums)
    ldlist <- .C("getLinkDensities", as.integer(hcedges$merge[,
        1]), as.integer(hcedges$merge[, 2]), as.integer(edges[,
        1]), as.integer(edges[, 2]), as.integer(len), as.integer(clusnums),
        as.integer(numcl), pdens = double(length(hh)), heights = as.double(hh),
        pdmax = double(1), csize = integer(1), as.logical(removetrivial),
        as.logical(bipartite), as.integer(bip), as.logical(verbose))
    pdens <- c(0, ldlist$pdens)
    heights <- c(0, hh)
    pdmax <- ldlist$pdmax
    csize <- ldlist$csize
    if (csize == 0) {
        stop("\nno clusters were found in this network; maybe try a larger network\n")
    }
    if (verbose) {
        cat("\n   Maximum partition density = ", max(pdens),
            "\n")
    }
    clus <- list()
    for (i in 1:csize) {
        if (verbose) {
            mes <- paste(c("   Finishing up...1/4... ", floor((i/csize) *
                100), "%"), collapse = "")
            cat(mes, "\r")
            flush.console()
        }
        clus[[i]] <- scan(file = "linkcomm_clusters.txt", nlines = 1,
            skip = i - 1, quiet = TRUE)
    }
    file.remove("linkcomm_clusters.txt")
    ecn <- data.frame()
    ee <- data.frame()
    lclus <- length(clus)
    for (i in 1:lclus) {
        if (verbose) {
            mes <- paste(c("   Finishing up...2/4... ", floor((i/lclus) *
                100), "%"), collapse = "")
            cat(mes, "\r")
            flush.console()
        }
        ee <- rbind(ee, cbind(el[clus[[i]], ], i))
        nodes <- node.names[unique(c(edges[clus[[i]], ]))]
        both <- cbind(nodes, rep(i, length(nodes)))
        ecn <- rbind(ecn, both)
    }
    colnames(ecn) <- c("node", "cluster")
    colnames(ee) <- c("node1", "node2", "cluster")
    ss <- NULL
    unn <- unique(ecn[, 2])
    lun <- length(unn)
    for (i in 1:length(unn)) {
        if (verbose) {
            mes <- paste(c("   Finishing up...3/4... ", floor((i/lun) *
                100), "%"), collapse = "")
            cat(mes, "\r")
            flush.console()
        }
        ss[i] <- length(which(ecn[, 2] == unn[i]))
    }
    names(ss) <- unn
    ss <- sort(ss, decreasing = T)
    unn <- unique(ecn[, 1])
    iecn <- as.integer(as.factor(ecn[, 1]))
    iunn <- unique(iecn)
    lunn <- length(iunn)
    nrows <- nrow(ecn)
    oo <- rep(0, lunn)
    oo <- .C("getNumClusters", as.integer(iunn), as.integer(iecn),
        counts = as.integer(oo), as.integer(lunn), as.integer(nrows),
        as.logical(verbose))$counts
    names(oo) <- unn
    if (verbose) {
        cat("\n")
    }
    pdplot <- cbind(heights, pdens)
    missnames <- setdiff(node.names, names(oo))
    m <- rep(0, length(missnames))
    names(m) <- missnames
    oo <- append(oo, m)
    all <- list()
    all$numbers <- c(len, nnodes, length(clus))
    all$hclust <- hcedges
    all$pdmax <- pdmax
    all$pdens <- pdplot
    all$nodeclusters <- ecn
    all$clusters <- clus
    all$edges <- ee
    all$numclusters <- sort(oo, decreasing = TRUE)
    all$clustsizes <- ss
    all$igraph <- graph.edgelist(el, directed = directed)
    all$edgelist <- el
    all$directed <- directed
    all$bipartite <- bipartite
    class(all) <- "linkcomm"
    if (plot) {
        if (verbose) {
            cat("   Plotting...\n")
        }
        if (len < 1500) {
            if (len < 500) {
                all <- plot(all, type = "summary", droptrivial = removetrivial,
                  verbose = verbose)
            }
            else {
                all <- plot(all, type = "summary", right = FALSE,
                  droptrivial = removetrivial, verbose = verbose)
            }
        }
        else if (len <= edglim) {
            oldpar <- par(no.readonly = TRUE)
            par(mfrow = c(1, 2), mar = c(5.1, 4.1, 4.1, 2.1))
            plot(hcedges, hang = -1, labels = FALSE)
            abline(pdmax, 0, col = "red", lty = 2)
            plot(pdens, heights, type = "n", xlab = "Partition Density",
                ylab = "Height")
            lines(pdens, heights, col = "blue", lwd = 2)
            abline(pdmax, 0, col = "red", lty = 2)
            par(oldpar)
        }
        else {
            plot(heights, pdens, type = "n", xlab = "Height",
                ylab = "Partition Density")
            lines(heights, pdens, col = "blue", lwd = 2)
            abline(v = pdmax, col = "red", lwd = 2)
        }
    }
    return(all)
}
