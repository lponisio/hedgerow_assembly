#!/usr/bin/env bash

## compile the python functions for the change point analysis
cd changePoint/cptPeel/fitHRG_GPL_Bayes/
make cleanall; make
cd ../consensusHRG_GPL_Bayes/
make cleanall; make
cd ../../


##
RScript dataPrep.R

for i in `seq 100 120`; do
    ## runs in parallel on 2 cores
    python hedgerows.py

    ## convert data to a helpful form
    filename=run_${i}.csv
    Rscript prepChangePointOutput.R $filename --save

    mv cptPeel/LogLs_4.txt  cptPeel/LogLs_4.{$i}.txt
    mv cptPeel/results_4.txt  cptPeel/results_4.{$i}.txt
    rm cptPeel/LogLs_4.txt cptPeel/results_4.txt
done


## create consensus trees
Rscript consensusChangePoints.R
python postChangePoint.py saved/consensus.txt
# python convertConsensusTrees.py saved/lastyr_consensus.txt

## *************************************************************
## the above works for running on a single computer, for the cluster
## we had trouble with the file processing so is in done in a batch
## here

# Rscript prepChangePointOutput_cluster.R

## then run from create consensus trees
