
## if this were in R I would list all the files then divide by using
## grep on the names.... then loop over the name groups 

import glob
from subprocess import call

## find .pair files
l = glob.glob("baci/*_*.pairs")

prefixes = []
for file in l:
    s = file.split("_")
    fname = s[0]
    fname = fname[5:]
    prefixes.append(fname)

prefixes = set(prefixes)

## run change point analysis for each site 
for prefix in prefixes:
    query = "baci/" + prefix + "_*.pairs"
    pairlist0 = glob.glob(query)
    pairlist = []
    for pair in pairlist0:
        pairlist = pair[5:]
    line = "runNetworkChangePoint.py graph-names.lut 4 "
    for f in pairlist:
        line += f
        line += " "
    line = line + "-p 'baci/'"
    print(line)
    call(line)

