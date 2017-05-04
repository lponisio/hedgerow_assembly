#!/usr/bin/env bash

for i in `seq 1 10`; do
    rm -rf cptPeel/baci
    mkdir -p cptPeel/baci

    RScript dataPrep.R
    ## remove old resuts because the change point analysis appends them

    ## otherwise
    rm cptPeel/LogLs_4.txt cptPeel/results_4.txt
    ## runs in parallel on 2 cores
    python hedgerows.py

    ## convert data to a helpful form
    fileshort=test_${i}.csv
    RScript prepChangePointOutput.R $fileshort --save
done





## binomial comparisons
# RScript comparisons.R

# ## create consensus trees
# python postChangePoint.py saved/consensus.txt
# python convertConsensusTrees.py saved/lastyr_consensus.txt
# RScript plotting/networks.R
