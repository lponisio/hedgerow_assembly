#!/usr/bin/env python

import os
import sys
import subprocess as sp
import csv

# paradigms_file = sys.argv[1]
paradigms_file = "saved/lastyr_consensus.txt"
fh = open(paradigms_file, 'r')


fh =[i[0].strip() for i in csv.reader(fh, delimiter=' ') if i != ['', '']]

paradigms = [l.strip() for l in fh]

os.chdir('cptPeel')
import convertConsensus as cc

for paradigm in paradigms:
    if os.path.isfile('%s' %paradigm):
        cc.writeToFile('%s'  %paradigm, 'gml')
