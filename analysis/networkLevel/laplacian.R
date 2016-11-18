rm(list=ls())
library(lme4)
library(lmerTest)
setwd('~/Dropbox/hedgerow_assembly/analysis/networkLevel')

load('../../data/networks/allSpecimens.Rdata')
f.path <- "../changePoint/cptPeel/baci"
load(file=file.path(f.path, "graphs.Rdata"))

source('src/laplacian_functions.R')

fig.path <- '../speciesLevel/Figures'
sites <- sapply(strsplit(names(nets), "[.]"), function(x) x[1])
years <- sapply(strsplit(names(nets), "[.]"), function(x) x[2])

status.table <- table(spec$Site, spec$SiteStatus)
status.table <- as.data.frame(cbind(rownames(status.table),
                                    colnames(status.table)[apply(status.table,
                                                            1, which.max)]))
colnames(status.table)<- c("Site", "SiteStatus")

all.alg.Con <- t(do.call(cbind.data.frame, lapply(nets, algCone)))
all.alg.Con <- as.data.frame(cbind(sites,years, all.alg.Con))
colnames(all.alg.Con) <- c("Site", "Year","Ncomp", "AlgCon", "EigenRatio")

## add status
all.alg.Con.status <- merge(x=all.alg.Con,
                            y= status.table,
                            by.x="Site",
                            by.y="Site")

## add ypr
all.alg.Con.status$ypr <- spec$ypr[match(paste(all.alg.Con.status$Site,
                                               all.alg.Con.status$Year),
                                         paste(spec$Site, spec$Year))]

baci.sites <- c("MullerB", "Sperandio", "Barger", "Butler", "Hrdy")

all.alg.Con.status$ypr[!all.alg.Con.status$Site %in% baci.sites] <- NA

## change AlgCon to numeric (it is a factor for some reason..)
all.alg.Con.status$AlgCon <- as.numeric(as.character(
  all.alg.Con.status$AlgCon))
all.alg.Con.status$EigenRatio <- as.numeric(as.character(
  all.alg.Con.status$EigenRatio))
all.alg.Con.status$Ncomp <- as.numeric(as.character(
  all.alg.Con.status$Ncomp))

alg.con.mod <- lmer(AlgCon ~ SiteStatus +
                    (1|Site) + (1|Year),    
                    data=all.alg.Con.status)

summary(alg.con.mod)

# Algebrain connectivity x Ypr
alg.con.mod.ypr <- lmer(AlgCon ~ ypr + 
                        (1|Site) + (1|Year),    
                        data=all.alg.Con.status)
summary(alg.con.mod.ypr)

print(summary(alg.con.mod.ypr))

# Eigenvalue ratio  x Ypr
EigenRatio.mod.ypr <- lmer(EigenRatio ~ ypr +  Ncomp +
                          (1|Site) + (1|Year),    
                        data=all.alg.Con.status)
summary(EigenRatio.mod.ypr)

save(alg.con.mod.ypr, all.alg.Con.status, EigenRatio.mod.ypr,
     file="saved/mods/AlgCon.Rdata")

## look at correlation between AlgCon and species richness

net.mets <- sapply(nets, function(x) sum(dim(empty(x))))

all.alg.Con.status$tot.rich <- net.mets[match(names(net.mets),
                                              paste(all.alg.Con.status$Site,
                                                    all.alg.Con.status$Year,
                                                    sep="."))]

cor.test(all.alg.Con.status$tot.rich, all.alg.Con.status$AlgCon)

cor.test(all.alg.Con.status$tot.rich[!is.na(all.alg.Con.status$ypr)],
         all.alg.Con.status$AlgCon[!is.na(all.alg.Con.status$ypr)])
