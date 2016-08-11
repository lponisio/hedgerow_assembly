source('plotting/src/predictIntervals.R')
source('plotting/src/CIplotting.R')
source('plotting/src/plotPanels.R')
load(file='saved/corMets.Rdata')

dd <- expand.grid(ypr=seq(from= min(cor.dats$ypr, na.rm=TRUE),
                          to= max(cor.dats$ypr, na.rm=TRUE),
                          length=10))

## ************************************************************
## nodf
## ************************************************************
load(file='saved/mods/baci_nodf.Rdata')

dd.nodf <- cbind(dd, zNODF=0)

nodf.pi <- predict.int(mod= baci.nodf.mod,
                        dd=dd.nodf,
                        y="zNODF",
                        family="gaussian")

## plot.predict.ypr(new.dd=nodf.pi,
##                  ylabel="Nestedness",
##                  dats=cor.dats,
##                  y1="zNODF")

## ************************************************************
## modularity
## ************************************************************
load(file='saved/mods/baci_mod.Rdata')

dd.mod <- cbind(dd, zmod.met.D=0)

mod.pi <- predict.int(mod= baci.mod.mod,
                        dd=dd.mod,
                        y="zmod.met.D",
                        family="gaussian")

## plot.predict.ypr(new.dd=mod.pi,
##                  ylabel="Modularity",
##                  dats=cor.dats,
##                  y1="zmodD")

## ************************************************************
## specialization
## ************************************************************
load(file='saved/mods/baci_h2.Rdata')

dd.h2 <- cbind(dd, H2=0)

h2.pi <- predict.int(mod= baci.h2.mod,
                        dd=dd.h2,
                        y="H2",
                        family="gaussian")

## plot.predict.ypr(new.dd=h2.pi,
##                  ylabel="Specialization",
##                  dats=cor.dats,
##                  y1="zH2")


## ************************************************************
## connectance
## ************************************************************
load(file='saved/mods/baci_conn.Rdata')

dd.conn <- cbind(dd, connectance=0)

conn.pi <- predict.int(mod= baci.conn.mod,
                        dd=dd.conn,
                        y="connectance",
                        family="gaussian")

## ************************************************************
## plotting
## ************************************************************
plot.panels()

