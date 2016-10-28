rm(list=ls())
setwd('~/Dropbox/hedgerow_assembly/analysis/networkLevel')
source('src/initialize.R')
load('../../data/networks/baci_networks_years.Rdata')
N <- 99

## ## ************************************************************
## ## calculate metrics and zscores
## ## ************************************************************
## mets <- lapply(nets, network.metrics, N)

## cor.dats <- prep.dat(mets,  spec)

## cor.dats$tot.rich <- cor.dats$number.of.species.LL +
##   cor.dats$number.of.species.HL

## save(cor.dats, file='saved/corMets.Rdata')

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
baci.mod.mod <- lmer(zmod.met.D ~ scale(ypr) +
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


## species richness pol
baci.rich.hl.mod <- glmer(number.of.species.HL ~ scale(ypr) +
                 (1|Site) + (1|Year), family="poisson",
                 data=cor.dats[!is.na(cor.dats$ypr),])
summary(baci.rich.hl.mod)
save(baci.rich.hl.mod, file='saved/mods/baci_rich_hl.Rdata')

## species richness plants
baci.rich.ll.mod <- glmer(number.of.species.LL ~ scale(ypr) +
                 (1|Site) + (1|Year), family="poisson",
                 data=cor.dats[!is.na(cor.dats$ypr),])
summary(baci.rich.ll.mod)
save(baci.rich.ll.mod, file='saved/mods/baci_rich_ll.Rdata')

## total species richness
baci.rich.tot.mod <- glmer(tot.rich ~ scale(ypr) +
                 (1|Site) + (1|Year), family="poisson",
                 data=cor.dats[!is.na(cor.dats$ypr),])
summary(baci.rich.tot.mod)
save(baci.rich.tot.mod, file='saved/mods/baci_rich_tot.Rdata')

## connectance
baci.conn.mod <- lmer(connectance ~ scale(ypr) +
                 (1|Site) + (1|Year),
                 data=cor.dats[!is.na(cor.dats$ypr),])
summary(baci.conn.mod)
save(baci.conn.mod, file='saved/mods/baci_conn.Rdata')

## weighted connectance
baci.wconn.mod <- lmer(weighted.connectance ~ scale(ypr) +
                 (1|Site) + (1|Year),
                 data=cor.dats[!is.na(cor.dats$ypr),])
summary(baci.wconn.mod)
save(baci.wconn.mod, file='saved/mods/baci_wconn.Rdata')

## niche overlap pollinators
baci.no.pol.mod <- lmer(niche.overlap.HL ~ scale(ypr) +
                 (1|Site) + (1|Year),
                 data=cor.dats[!is.na(cor.dats$ypr),])
summary(baci.no.pol.mod)
save(baci.no.pol.mod, file='saved/mods/baci_no_pol.Rdata')

## niche overlap plants
baci.no.plant.mod <- lmer(niche.overlap.LL ~ scale(ypr) +
                 (1|Site) + (1|Year),
                 data=cor.dats[!is.na(cor.dats$ypr),])
summary(baci.no.plant.mod)
save(baci.no.plant.mod, file='saved/mods/baci_no_plant.Rdata')


source('plotting/baci.R')


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
