#!/usr/bin/env bash

## @knitr external_chunk_1
## prep data
## only needed for original analysis, the saved .Rdata files should
## all in in github

## RScript ../dataPrep/dataPrep.R

##***************************************************************
## changepoint point analysis
##***************************************************************
bash changePoint/mainChangePoint.sh
python changePoint/convertfiles.py
RScript changePoint/plotting/networks.R

RScript changePoint/plotting/networks.R
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
## network level metrics
##***************************************************************
## pollinator turnover
Rscript variability/nulls.R "pols"
Rscript variability/beta-div.R "pols"

## interation turnover
Rscript variability/nulls.R "int"
Rscript variability/beta-div.R "int"

## plant turnover
Rscript variability/nulls.R "plants"
Rscript variability/beta-div.R "plants"
