from subprocess import call
import os

def runConsensus(networkFileSequence,namesFile,window,path):
    namesFile=os.path.join(path,namesFile)
    if window>1:
        paramStr="-n %s -m %i" % (namesFile,len(networkFileSequence))
    else: 
        paramStr="-n %s -f" % namesFile

    for networkFile in networkFileSequence:
        paramStr+=" %s" % (os.path.join(path,networkFile))

    #~ paramStr+=" > /dev/null"
    outputStr=" > /dev/null"
    #~ outputStr=" > log.txt"
    #~ paramStr+=" > log.txt"

    if not os.path.isfile(os.path.join(path,networkFileSequence[-1].strip(".pairs")+"_joint_best-dendro.hrg")):
        print "Run: ./fitHRG_GPL_Bayes/fitHRG -t joint %s%s" % (paramStr,outputStr)
        call("./fitHRG_GPL_Bayes/fitHRG -t joint %s%s" % (paramStr,outputStr), shell=True)

    #~ if not os.path.isfile(os.path.join(path,networkFileSequence[-1].strip(".pairs")+"_joint_s1-consensus.tree")):
    if window==1:
        paramStr=paramStr.split('-f')[0]
    print "Run: ./consensusHRG_GPL_Bayes/consensusHRG %s -f %s_joint_best-dendro.hrg" % (paramStr,os.path.join(path,networkFileSequence[-1].strip(".pairs")))
    call("./consensusHRG_GPL_Bayes/consensusHRG %s -f %s_joint_best-dendro.hrg %s" % (paramStr,os.path.join(path,networkFileSequence[-1].strip(".pairs")),outputStr), shell=True)







if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Infer GHRG model for a sequence of networks")
    
    parser.add_argument("nodenamesfile", help="Input node names file e.g. names.lut")
    parser.add_argument("windowsize", help="Length of sliding window")
    parser.add_argument("networkfilesequence", help='Input sequence of network files: e.g. "network1.pairs network2.pairs network3.pairs"')
    parser.add_argument("-p","--path", help="Path to files (if not in current directory)")
    args = parser.parse_args()




    import changepointDetection

    if args.path:
        path=args.path
    else:
        path="."

    networkFileSequence=args.networkfilesequence.split()
    
    #~ cpDetector = changepointDetection.anomalyDetection()
    window = int(args.windowsize)
    runConsensus(networkFileSequence,args.nodenamesfile,window,path)
