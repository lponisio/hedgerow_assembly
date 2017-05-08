
## setwd('~/Dropbox/hedgerow_assembly/analysis/networkLevel')
setwd('networkLevel')
source('src/initialize.R')
load('../../data/networks/all_networks_years.Rdata')

## either "abund" or "degree"
 extinction.method <- "degree"

## **********************************************************
## robustness
## **********************************************************
## simulate plant extinction

res <- simExtinction(nets, extinction.method, spec)

save(res, file=file.path(save.path,
            sprintf('resilience_%s.Rdata',
                    extinction.method)))

## no change in robustness by site status
## not included in ms jsut interesting
mod.status <- lmer(Robustness ~ SiteStatus
             + (1|Site) + (1|Year),
             data=res)
summary(mod.status)
save(mod.status, file=file.path(save.path,
            sprintf('mods/resilience_status_%s.Rdata',
                    extinction.method)))


## no effect of ypr on robustness
mod.ypr <- lmer(Robustness ~ ypr
             + (1|Site) + (1|Year),
             data=res[!is.na(res$ypr),])
summary(mod.ypr)
save(mod.ypr, file=file.path(save.path,
            sprintf('mods/resilience_ypr_%s.Rdata',
                    extinction.method)))

print(summary(mod.ypr))
