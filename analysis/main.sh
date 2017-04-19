#!/usr/bin/env bash

## prep data
## only needed for original analysis, the saved .Rdata files should
## all in in github
RScript ../dataPrep/dataPrep.R

## changepoint analysis
RScript changePoint/dataPrep.R
python changePoint/hedgerows.py
RScript changePoint/prepChangePointOutput.R
python changePoint/postChangePoint.py saved/consensus.txt
python changePoint/convertfiles.py
RScript changePoint/plotting/networks.R

## all of the other analyses, can be run one at a time in the mainR.R
## file

RScript mainR.R
