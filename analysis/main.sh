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
RScript changePoint/plotting/networks.R
Rscript changePoint/comparisons.R

##***************************************************************
## coefficient of variation of closeness
##***************************************************************
Rscript variability/cv.R
Rscript variability/plotting/cv.R

##***************************************************************
## coefficient of variation of closeness
##***************************************************************
Rscript speciesLevel/specMets.R



Rscript speciesLevel/plotting/specMets.R

python changePoint/convertfiles.py
RScript changePoint/plotting/networks.R

## all of the other analyses, can be run one at a time in the mainR.R
## file

RScript mainR.R
