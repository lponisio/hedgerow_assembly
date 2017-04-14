#!/usr/bin/env python
import sys
import subprocess as sp

paradigms_file = sys.argv[1]
fh = open(paradigms_file, 'r')

paradigms = [l.strip() for l in fh]

for paradigm in paradigms:
    command = "python runConsensus.py graph-names.lut {} '{}' -p 'baci/'"
    command = comamnd.format(len(paradigm.split()), paradigm)
    job = sp.Popen(command, shell=True)
    job.wait()
