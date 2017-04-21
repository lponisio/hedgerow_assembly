#!/usr/bin/env bash
RScript dataPrep.R
## remove old resuts because the change point analysis appends them
## otherwise
rm cptPeel/LogLs_4.txt cptPeel/results_4.txt
## runs in parallel on 2 cores
python hedgerows.py
RScript prepChangePointOutput.R
python postChangePoint.py saved/consensus.txt
python convertConsensus.py saved/lastyr_consensus.txt
# python convertfiles.py
RScript plotting/networks.R
