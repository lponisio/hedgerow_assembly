rm(list=ls())
## prep data
## 1) drop non bee/syrphid data
## 2) create networks
## 3) calculate network metrics for whole landscape
## 4) add trait data to spec sheet
source("~/Dropbox/hedgerow_assembly/dataPrep/dataPrep.R")

##***************************************************************
## change points
##***************************************************************
## create graphs for change point analysis
source('~/Dropbox/hedgerow_assembly/analysis/changePoint/dataPrep.R')

## run python code for calculating network change points
## plot networks and network communities
source('~/Dropbox/hedgerow_assembly/analysis/changePoint/plotting/networks.R')

## binomial analysis in comparisons.R
## exploratory analysis in  core.R

##***************************************************************
## species level metrics
##***************************************************************
## calcualte species level metrics, run linear models
source('~/Dropbox/hedgerow_assembly/analysis/speciesLevel/specMets.R')
## make species level metric plots
source('~/Dropbox/hedgerow_assembly/analysis/speciesLevel/plotting/specMets.R')

##***************************************************************
## network level metrics
##***************************************************************
## network metrics through time
source('~/Dropbox/hedgerow_assembly/analysis/networkLevel/baci.R')

## robustness
## either "abund" or "degree"
extinction.method <- "degree"
source('~/Dropbox/hedgerow_assembly/analysis/networkLevel/resilience.R')
extinction.method <- "abund"
source('~/Dropbox/hedgerow_assembly/analysis/networkLevel/resilience.R')
## sensitivity to perturbation
source('~/Dropbox/hedgerow_assembly/analysis/networkLevel/laplacian.R')

## plotting
source('~/Dropbox/hedgerow_assembly/analysis/networkLevel/plotting/baci.R')
source('~/Dropbox/hedgerow_assembly/analysis/networkLevel/plotting/resilience.R')

##***************************************************************
## turnover metrics
##***************************************************************
## coefficient of variation
source('~/Dropbox/hedgerow_assembly/analysis/variability/cv.R')
source('~/Dropbox/hedgerow_assembly/analysis/variability/plotting/cv.R')

## nulls for computing stats
## pols, int or plants
type <- "plants"
source('~/Dropbox/hedgerow_assembly/analysis/variability/nulls.R')
## species turnover through years
source('~/Dropbox/hedgerow_assembly/analysis/variability/beta-div.R')

type <- "ints"
source('~/Dropbox/hedgerow_assembly/analysis/variability/nulls.R')
## species turnover through years
source('~/Dropbox/hedgerow_assembly/analysis/variability/beta-div.R')

type <- "pols"
source('~/Dropbox/hedgerow_assembly/analysis/variability/nulls.R')
## species turnover through years
source('~/Dropbox/hedgerow_assembly/analysis/variability/beta-div.R')
plotDistPanels()

## "link" turnover
source('~/Dropbox/hedgerow_assembly/analysis/variability/beta-link.R')

## species turnover through years
source('~/Dropbox/hedgerow_assembly/analysis/variability/beta-div.R')

## interaction partner variation
## source('~/Dropbox/hedgerow_assembly/analysis/variability/beta-disp.R')

## plotting
source('~/Dropbox/hedgerow_assembly/analysis/variability/plotting/beta-int.R')
