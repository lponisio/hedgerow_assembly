rm(list=ls())
library(MASS)
library(lme4)
library(nlme)

setwd("~/Dropbox/hedgerow_assembly/analysis/changePoint")

source('../../dataPrep/src/misc.R')
load('cptPeel/baci/graphs.Rdata')
load('../../data/networks/allSpecimens.Rdata')
samples <- read.csv("../../data/samples.csv")
dats <- read.csv('cptPeel/changing_points.csv')
BACI.site <- c('Barger', 'Butler', 'Hrdy', 'MullerB', 'Sperandio')

## **********************************************************
## binomial tests
## **********************************************************
chpt.trial <- aggregate(spec$Year, list(Site=spec$Site),
                           FUN=function(x) length(unique(x)))
chpt.trial$x <- chpt.trial$x - 1

chpt.trial <- chpt.trial[chpt.trial$x >= 4,]

change.points.site <- tapply(dats$cp, dats$sites, length)

chpt.trial$chpts <- change.points.site[match(chpt.trial$Site,
                                                rownames(change.points.site))]
chpt.trial$chpts[is.na(chpt.trial$chpts)] <- 0
colnames(chpt.trial) <- c("Site", "trial", "chpts")

chpt.trial$status <- spec$SiteStatus[match(chpt.trial$Site, spec$Site)]
chpt.trial$status[chpt.trial$Site %in% BACI.site] <- "maturing"

## binomial model with change points as successes
mod.chpt <- glm(cbind(chpt.trial$chpts, chpt.trial$trial - chpt.trial$chpts) ~
    chpt.trial$status, family="binomial")
summary(mod.chpt)

exp(cbind(coef(mod.chpt), confint(mod.chpt)))

# controls with change points
nrow(chpt.trial[chpt.trial$status == "control" &
                chpt.trial$chpts != 0,])/
  nrow(chpt.trial[chpt.trial$status == "control",])

nrow(chpt.trial[chpt.trial$status == "mature" &
                chpt.trial$chpts != 0,])/
  nrow(chpt.trial[chpt.trial$status == "mature",])

nrow(chpt.trial[chpt.trial$chpts != 0,])/
  nrow(chpt.trial)


## maturing has more successes than mature and controls, and mature
## and controls have about the same


## **********************************************************
## some methodological tests: are years with larger differences in the
## number of samples more likely to have a change point?
## **********************************************************

diff.samp <- matrix(NA, nrow=nrow(samples),
                    ncol=ncol(samples) -2)
colnames.prep <- character(ncol(diff.samp))

for(i in 3:ncol(samples)){
    diff.samp[,i -2] <- abs(samples[, i] -  samples[, i -1])
    colnames.prep[i-2] <- paste(colnames(samples)[i -1], "-",
                                colnames(samples)[i], sep="")
}

rownames(diff.samp) <- samples$X
colnames.prep <- gsub("X", "", colnames.prep)

colnames(diff.samp) <- colnames.prep

diff.dats <- comm.mat2sample(diff.samp)

diff.dats$cp <- match(paste(diff.dats$Site,
                                        diff.dats$Date),
                      paste(dats$sites, dats$cp))

diff.dats.cp <- diff.dats[diff.dats$Site %in% chpt.trial$Site,]
diff.dats.cp <- diff.dats.cp[diff.dats.cp$Samp !=0,]

diff.dats.cp$cp[is.na(diff.dats.cp$cp)] <- 0
diff.dats.cp$cp[diff.dats.cp$cp > 1] <- 1

mod.samp <- glmer(cp ~ Samp + (1|Site), data=diff.dats.cp, family="binomial")
summary(mod.samp)

## probability of change point does not increase with the difference
## in the number of samples between years

## **********************************************************
## poisson likelihood
## maybe not quite right given the number of trials differs, not
## included in ms but had similar results
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
