rm(list=ls())
setwd('~/Dropbox/hedgerow_assembly/analysis/networkLevel')
source('src/initialize.R')
N <- 999

## ************************************************************
## create site by year neworks and calculate statistics
## (in terminal)
## ************************************************************
## create pp matrix for each site, year
nets <- break.net(spec)
save(nets, file='../../data/networks/networks_years.Rdata')

## ************************************************************
## individuals nulls
## ************************************************************
nulls <- lapply(nets, vaznull, N=N)
save(nulls, file='saved/nulls/all.Rdata')

## ************************************************************
## calculate metrics and zscores
## ************************************************************
load(file='saved/nulls/all.Rdata')

mets <- lapply(nets, calc.metric)
null.mets <- rapply(nulls, calc.metric, how="replace")
null.mets <- lapply(null.mets, function(x) do.call(rbind, x))
load(file='saved/nullMets.Rdata')

cor.mets <- mapply(function(a, b)
                   cor.metrics(true.stat= a,
                               null.stat= b,
                               N=N),
                   a=mets,
                   b=null.mets,
                   SIMPLIFY=FALSE)

cor.dats <- prep.dat(cor.mets,  spec)

## ************************************************************
## niche overlap
## ************************************************************
no <- t(sapply(nets, networklevel, index="niche overlap"))

cor.dats$niche.overlap.pol <- no[, "niche.overlap.HL"][match(rownames(no),
                                     paste(cor.dats$Site,
                                     cor.dats$Year, sep="."))]
cor.dats$niche.overlap.plants <- no[, "niche.overlap.LL"][match(rownames(no),
                                     paste(cor.dats$Site,
                                     cor.dats$Year, sep="."))]

save(cor.dats, file='saved/corMets.Rdata')

## distribution is niche overlap

dis.mats <- lapply(lapply(nets, t), vegdist, method="chao")

layout(matrix(1:6, nrow=2))
cols <- rainbow(length(unique(cor.dats$Year)))
lapply(unique(cor.dats$Site), function(x){
  this.mats <- dis.mats[cor.dats$Site == x]
  plot(NA, ylim=c(0,10), xlim=c(0,1.5),
       ylab="Frequency",
       xlab="Niche Overlap",
       main= x)
  for(i in 1:length(this.mats)){
    points(density(this.mats[[i]]), col=cols[i], type="l", lwd=2)
  }
})

plot(NA, ylim=c(0,1), xlim=c(0,1), xaxt="n", yaxt="n", ylab="", xlab="")
legend("center", col=cols, lwd="2",
       legend=sort(unique(cor.dats$Year)),
       bty="n")


## ************************************************************
## effect of years post restoration
## ************************************************************
load(file='saved/corMets.Rdata')

## nestedness
baci.nodf.mod <- lmer(zNODF ~ scale(ypr) +
                 (1|Site) + (1|Year),
                 data=cor.dats[!is.na(cor.dats$ypr),])
summary(baci.nodf.mod)
save(baci.nodf.mod, file='saved/mods/baci_nodf.Rdata')

## modularity
baci.mod.mod <- lmer(zmodD ~ scale(ypr) +
                 (1|Site) + (1|Year),
                 data=cor.dats[!is.na(cor.dats$ypr),])
summary(baci.mod.mod)
save(baci.mod.mod, file='saved/mods/baci_mod.Rdata')

## h2
baci.h2.mod <- lmer(zH2 ~ scale(ypr) +
                 (1|Site) + (1|Year),
                 data=cor.dats[!is.na(cor.dats$ypr),])
summary(baci.h2.mod)
save(baci.h2.mod, file='saved/mods/baci_h2.Rdata')

## niche overlap pollinators
baci.no.pol.mod <- lmer(niche.overlap.pol ~ scale(ypr) +
                 (1|Site) + (1|Year),
                 data=cor.dats[!is.na(cor.dats$ypr),])
summary(baci.no.pol.mod)
save(baci.no.pol.mod, file='saved/mods/baci_no_pol.Rdata')

## niche overlap plants
baci.no.plant.mod <- lmer(niche.overlap.plants ~ scale(ypr) +
                 (1|Site) + (1|Year),
                 data=cor.dats[!is.na(cor.dats$ypr),])
summary(baci.no.plant.mod)
save(baci.no.plant.mod, file='saved/mods/baci_no_plant.Rdata')


source('plotting/baci.R')
