
## if this were in R I would list all the files then divide by using
## grep on the names.... then loop over the name groups

import os
import glob
import subprocess as sp
#from subprocess import call
from multiprocessing import Pool

def runNetworkChangePoint(queries):
    line = "python runNetworkChangePoint.py graph-names.lut 4 '{}' -p 'baci/'".format(' '.join(queries))
    job = sp.Popen(line, shell=True)
    job.wait()

def mapPool(threads, args):
    # args = [(function, (arg1, arg2), {keyword1:arg1, keyword2:arg2}), nextjob, ...]
    p = Pool(threads)
    output = p.map(star_func, args)
    p.close()
    p.join()
    return output

## cpaste
def star_func(arglist):
    func = arglist[0]
    args = arglist[1]
    try:
        kwargs = arglist[2]
        try:
            return func(*args, **kwargs)
        except:
            raise Exception("".join(traceback.format_exception(*sys.exc_info())))
    except IndexError:
        try:
            return func(*args)
        except:
            raise Exception("".join(traceback.format_exception(*sys.exc_info())))

## --

os.chdir(os.path.expanduser('~/Dropbox/hedgerow_assembly/analysis/changePoint/cptPeel'))

ncores = 2

## find .pair files
l = glob.glob("baci/*_*.pairs")
l = [os.path.basename(f) for f in l]

prefixes = [f.split('_')[0] for f in l]
prefixes = [p for p in prefixes if prefixes.count(p) >=5]
prefixes = set(prefixes)

## run change point analysis for each site
jobs = [(runNetworkChangePoint, [f for f in l if prefix in f]) for prefix in prefixes]
mapPool(ncores, jobs)


