#!/usr/bin/env python

import os
import sys
import subprocess as sp
import csv


paradigms_file = sys.argv[1]
fh = open(paradigms_file, 'r')


fh =[i[0].strip() for i in csv.reader(fh, delimiter=' ') if i != ['', '']]

paradigms = [l.strip() for l in fh]


# for paradigm in paradigms:

#     print paradigm
#     # command = "python cc.writeToFile {} 'gml'"
#     command = "python cc.writeToFile %s 'gml'" %paradigm
#     # command = command.format(len(paradigm.split()), paradigm)
#     job = sp.Popen(command, shell=True)
#     job.wait()

os.chdir('cptPeel')
import convertConsensus as cc

for paradigm in paradigms:
    if os.path.isfile('%s' %paradigm):
        cc.writeToFile('%s'  %paradigm, 'gml')
