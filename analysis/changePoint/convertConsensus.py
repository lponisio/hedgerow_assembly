#!/usr/bin/env python
import os
import sys
import subprocess as sp
import convertConsensus as cc

paradigms_file = sys.argv[1]
fh = open(paradigms_file, 'r')

paradigms = [l.strip() for l in fh]

os.chdir('cptPeel')

for paradigm in paradigms:
    command = "python cc.writeToFile {} {} 'gml'"
    command = command.format(len(paradigm.split()), paradigm)
    job = sp.Popen(command, shell=True)
    job.wait()

