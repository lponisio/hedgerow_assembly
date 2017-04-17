#!/usr/bin/env bash
RScript dataPrep.R
python hedgerows.py
RScript prepChangePointOutput.R
python postChangePoint.py saved/consensus.txt
