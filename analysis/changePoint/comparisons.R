rm(list=ls())
setwd("~/Dropbox/hedgerow_assembly/analysis/changePoint")
load('cptPeel/baci/graphs.Rdata')
load('../../data/networks/allSpecimens.Rdata')
library(MASS)
library(nlme)
dats <- read.csv('cptPeel/changing_points.csv')
BACI.site <- c('Barger', 'Butler', 'Hrdy', 'MullerB', 'Sperandio')
## **********************************************************
## binomial 
## **********************************************************
count.samples <- aggregate(spec$Year, list(Site=spec$Site),
                           FUN=function(x) length(unique(x)))
count.samples$x <- count.samples$x - 1

change.points.site <- tapply(dats$cp, dats$sites, length)

chpt.trial <- data.frame(chpts=change.points.site,
                         trials=count.samples$x[match(
                           rownames(change.points.site),
                           count.samples$Site)])

chpt.trial$status <- spec$SiteStatus[match(rownames(chpt.trial), spec$Site)]
chpt.trial$status[rownames(chpt.trial) %in% BACI.site] <- "maturing"

## binomial model with change points as successes
mod.chpt <- glm(cbind(chpt.trial$chpts, chpt.trial$trials - chpt.trial$chpts) ~
    chpt.trial$status, family="binomial")
summary(mod.chpt)

## maturing has more successes than mature and controls, and mature
## and controls have about the same

## **********************************************************
## poisson likelihood
## maybe not quite right given the number of trails differs
## **********************************************************
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
fit.cont
fit.maturing <- fitdistr(chpts.sites[statuses == "maturing"],
                     densfun="Poisson")
fit.maturing
fit.mature <- fitdistr(chpts.sites[statuses == "mature"],
                     densfun="Poisson")
fit.mature

## likelihood ratio text of fit of control model on maturing data,
## mature model on maturing data
lik.cont.mat <- log(prod(dpois(chpts.sites[statuses == "maturing"],
                           lambda=fit.cont$estimate)))
lik.mature.mat <- log(prod(dpois(chpts.sites[statuses == "maturing"],
                           lambda=fit.mature$estimate)))
lik.mature.cont <- log(prod(dpois(chpts.sites[statuses == "mature"],
                           lambda=fit.cont$estimate)))

## maturing is unlikely to be drawn from control
pchisq(-2*(lik.cont.mat - fit.maturing$loglik), 1, lower.tail=FALSE)

## or mature
pchisq(-2*(lik.mature.mat - fit.maturing$loglik), 1, lower.tail=FALSE)

## but mature could be draw from control
pchisq(-2*(lik.mature.cont - fit.mature$loglik), 1, lower.tail=FALSE)
