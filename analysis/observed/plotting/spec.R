library(RColorBrewer)
source('plotting/src/predictIntervals.R')
source('plotting/src/CIplotting.R')
source('plotting/src/plotPanels.R')
load(file="saved/spec.Rdata")

lapply(metrics, plotSpecs,
       diff.spec = diff.spec[diff.spec$speciesType ==
                                      "pollinator",],
       type="pol",
       path.f="figures/specialization")



lapply(metrics, plotSpecs,
       diff.spec = diff.spec[diff.spec$speciesType ==
                                      "plant",],
       type="plant",
       path.f="figures/specialization")
