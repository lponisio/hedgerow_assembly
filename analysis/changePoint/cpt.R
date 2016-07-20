rm(list=ls())
setwd("~/Dropbox/hedgerow_assembly/analysis/changePoint")
load('cptPeel/baci/graphs.Rdata')
load('../../data/networks/allSpecimens.Rdata')
dats <- read.csv('cptPeel/changing_points.csv')

counts <- table(dats[,-3])

chpts.yrs <- colSums(counts)
barplot(chpts.yrs)

sites <- unique(sapply(strsplit(names(graphs), "_"), function(x)
                       x[1]))
no.chpt <- sites[!sites %in% rownames(counts)]
add.chpt <- matrix(0, nrow=length(no.chpt), ncol=ncol(counts))
rownames(add.chpt) <- no.chpt

counts <- rbind(counts, add.chpt)

chpts.sites <- rowSums(counts)

statuses <- spec$SiteStatus[match(names(chpts.sites), spec$Site)]
BACI.site <- c('Barger', 'Butler', 'Hrdy', 'MullerB', 'Sperandio')
statuses[names(chpts.sites) %in% BACI.site] <- "maturing"

layout(matrix(1:3, nrow=1))
cont <- hist(chpts.sites[statuses == "control"], prob=TRUE,
     main="Unrestored", ylim=c(0,1),  xlim=c(0,4),
     breaks=0:4)
maturing <- hist(chpts.sites[statuses == "maturing"], prob=TRUE,
     main="Maturing", ylim=c(0,1), xlim=c(0,4),
     breaks=0:4)
mature <- hist(chpts.sites[statuses == "mature"], prob=TRUE,
     main="Mature", ylim=c(0,1), xlim=c(0,4),
     breaks=0:4)

## fit poisson distributions to data
fit.cont <- fitdistr(chpts.sites[statuses == "control"],
                     densfun="Poisson")
fit.maturing <- fitdistr(chpts.sites[statuses == "maturing"],
                     densfun="Poisson")
fit.mature <- fitdistr(chpts.sites[statuses == "mature"],
                     densfun="Poisson")

## likelihood ratio text of fit of control model on maturing data,
## mature model on maturing data
lik.cont.mat <- log(prod(dpois(chpts.sites[statuses == "maturing"],
                           lambda=fit.cont$estimate)))
lik.mat.mat <- log(prod(dpois(chpts.sites[statuses == "maturing"],
                           lambda=fit.mature$estimate)))

1- pchisq(-2*(lik.cont.mat - fit.maturing$loglik), 1)

1- pchisq(-2*(lik.mat.mat - fit.maturing$loglik), 1)
