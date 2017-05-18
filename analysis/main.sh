#!/usr/bin/env bash

## All paths are relative to the analysis folder within the github
## repo

## prep data
## only needed for original analysis, the saved .Rdata files should
## all in in github

## Rscript ../dataPrep/dataPrep.R

##***************************************************************
## changepoint point analysis
##***************************************************************
bash changePoint/mainChangePoint.sh
Rscript changePoint/plotting/networks.R
Rscript changePoint/comparisons.R

##***************************************************************
## coefficient of variation of closeness
##***************************************************************
Rscript variability/cv.R
Rscript variability/plotting/cv.R

##***************************************************************
## species level metrics
##***************************************************************
Rscript speciesLevel/specMets.R
Rscript speciesLevel/plotting/specMets.R

##***************************************************************
## turnover
##***************************************************************
## pollinator turnover
Rscript variability/nulls.R "pols"
Rscript variability/beta-div.R "pols"

## interation turnover
Rscript variability/nulls.R "ints"
Rscript variability/beta-div.R "ints"

## plant turnover
Rscript variability/nulls.R "plants"
Rscript variability/beta-div.R "plants"

## weighted link turnover
Rscript variability/beta-link.R
Rscript variability/plotting/beta-int.R

##***************************************************************
## network level analysis
##***************************************************************
## network metrics through assembly
Rscript networkLevel/baci.R
Rscript networkLevel/plotting/baci.R

## resilence through assembly
Rscript networkLevel/resilience.R "degree"
Rscript networkLevel/resilience.R "abund"

## algebraic connectance through time
Rscript networkLevel/laplacian.R

Rscript networkLevel/plotting/resilience.R

