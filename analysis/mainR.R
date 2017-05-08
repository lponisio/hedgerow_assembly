rm(list=ls())
analysis.dir <- getwd()

##***************************************************************
## change points
##***************************************************************

## data is preped and change point analysis is run through the bash
## script mainChangePoint.sh

## binomial analysis in comparisons.R, no resulting plots
source(file.path(analysis.dir,
                 'changePoint/comparisons.R'))
print("change point comparisons completed")
## exploratory analysis in  src/core.R

##***************************************************************
## closeness CV
##***************************************************************

## coefficient of variation
source(file.path(analysis.dir,
                 'variability/cv.R'))
print("variability analysis completed")
source(file.path(analysis.dir,
                 'variability/plotting/cv.R'))

##***************************************************************
## species level metrics
##***************************************************************

## calcualte species level metrics, run linear models
source(file.path(analysis.dir,
                 'speciesLevel/specMets.R'))
## make species level metric plots
source(file.path(analysis.dir,
                 'speciesLevel/plotting/specMets.R'))


##***************************************************************
## turnover metrics
##***************************************************************

## nulls for computing stats
## pols, int or plants
type <- "plants"
source(file.path(analysis.dir,
                 'variability/nulls.R'))
## species turnover through years
source(file.path(analysis.dir,
                 'variability/beta-div.R'))

type <- "ints"
source(file.path(analysis.dir,
                 'variability/nulls.R'))
## species turnover through years
source(file.path(analysis.dir,
                 'variability/beta-div.R'))

type <- "pols"
source(file.path(analysis.dir,
                 'variability/nulls.R'))
## species turnover through years
source(file.path(analysis.dir,
                 'variability/beta-div.R'))
plotDistPanels())

## "link" turnover
source(file.path(analysis.dir,
                 'variability/beta-link.R'))

## species turnover through years
source(file.path(analysis.dir,
                 'variability/beta-div.R'))

## plotting
source(file.path(analysis.dir,
                 'variability/plotting/beta-int.R'))


##***************************************************************
## network level metrics
##***************************************************************
## network metrics through time
source(file.path(analysis.dir,
                 'networkLevel/baci.R'))

## robustness
## either "abund" or "degree"
extinction.method <- "degree"
source(file.path(analysis.dir,
                 'networkLevel/resilience.R'))
extinction.method <- "abund"
source(file.path(analysis.dir,
                 'networkLevel/resilience.R'))

## sensitivity to perturbation
source(file.path(analysis.dir,
                 'networkLevel/laplacian.R'))

## plotting
source(file.path(analysis.dir,
                 'networkLevel/plotting/baci.R'))
source(file.path(analysis.dir,
                 'networkLevel/plotting/resilience.R'))
