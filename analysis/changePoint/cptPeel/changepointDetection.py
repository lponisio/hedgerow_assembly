"""change-pointDetection.py - detects change-points within sequences of graphs
    Copyright (C) 2014 Leto Peel

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA"""

from sys import stdout
import dendropy
from dendropy.calculate import treemeasure
import numpy as np
import os
from scipy.special import gammaln
from multiprocessing import Pool
from subprocess import call
from scipy.stats import norm
from collections import Counter
import time


#~ def unwrap_self_sampleNullModel(arg, **kwarg):
    #~ return sampleNullModel(*arg, **kwarg)


class Filenames(object):
    def __init__(self,consensusFilename,path,currentNormalNetworks,namesFile):
        self.consensusFilename = consensusFilename
        self.path = path
        self.currentNormalNetworks = currentNormalNetworks
        self.namesFile = namesFile

class CountModel(object):
    def __init__(self):
        self.alpha = {}
        self.beta = {}


class anomalyDetection(object):

    def __init__(self):
        self.taxon_namespace = None
        self.namesLUT = None
        self.nsamples = 1000
        #~ self.nsamples = 1
        self.window=None


    def detectAnomaliesInSequence(self,networkFileSequence,namesFile,window,path="."):
        self.namesFile=os.path.join(path,namesFile)
        self.window=window
        self.importNames(self.namesFile)
        self.alpha=0.95
        t=0
        normal=True
        print "start"
        while t<len(networkFileSequence)-window:
            #~ print "test", test,"window at t=", t
            startindx=t
            endindx=t+window

            paramStr="-n %s -m %i" % (self.namesFile,window)

            for gi in xrange(startindx,endindx):
                paramStr+=" %s" % (os.path.join(path,networkFileSequence[gi]))

            #~ paramStr+=" > /dev/null"
            outputStr=" > /dev/null"
            #~ paramStr+=" > log.txt"

            if not os.path.isfile(os.path.join(path,networkFileSequence[gi].strip(".pairs")+"_joint_best-dendro.hrg")):
                print "Run: ./fitHRG_GPL_Bayes/fitHRG -t joint %s%s" % (paramStr,outputStr)
                call("./fitHRG_GPL_Bayes/fitHRG -t joint %s%s" % (paramStr,outputStr), shell=True)

            if not os.path.isfile(os.path.join(path,networkFileSequence[gi].strip(".pairs")+"_joint_s1-consensus.tree")):
                print "Run: ./consensusHRG_GPL_Bayes/consensusHRG %s -f %s_joint_best-dendro.hrg" % (paramStr,os.path.join(path,networkFileSequence[gi].strip(".pairs")))
                call("./consensusHRG_GPL_Bayes/consensusHRG %s -f %s_joint_best-dendro.hrg" % (paramStr,os.path.join(path,networkFileSequence[gi].strip(".pairs"))), shell=True)

            #then run null model from ad.
            self.consensusFilename="%s_joint_s1-consensus.tree" % os.path.join(path,networkFileSequence[gi].strip(".pairs"))
            self.currentNormalNetworks=[networkFileSequence[gi] for gi in xrange(startindx,endindx)]
            self.createNullModel(path)
            self.testNetwork(path,self.currentNormalNetworks)

            changeDetected = self.calculateP_values()
            t+=1
            #~ normal =  (t < experimentLength) and not changeDetected
        #~ return tau, r, time.time()-stime




    def calculateP_values(self):
        stdout.write("\rCalculating p-value..." )
        stdout.flush()

        index = np.argmax(self.test_model)

        loglikelihoods = np.append(np.sort(self.null_model),np.inf)

        test_order = np.argsort(self.test_model)
        test_loglikeL = np.sort(self.test_model)


        m = 0
        curr_LL = -np.inf
        pvals = np.empty(len(test_loglikeL))

        for t in xrange(len(test_loglikeL)):
            stdout.write("\rCalculating p-value... \t%i%% " % (t*100/len(test_loglikeL)))
            stdout.flush()
            if not curr_LL>test_loglikeL[t]:
                m,curr_LL = ((c,ll) for c,ll in enumerate(loglikelihoods[m:],start=m) if ll>test_loglikeL[t]).next()
            pvals[test_order[t]] = m

        pvals/=self.nsamples

        self.pvals=pvals
        if self.window is not None:
            with open("results_%i.txt" % self.window, "a") as f:
                f.write("%s %s" % (self.currentNormalNetworks[0].strip(".pairs"),self.currentNormalNetworks[-1].strip(".pairs")))
                for p in pvals:
                    f.write(" %f" % p)
                f.write("\n")
            with open("LogLs_%i.txt" % self.window, "a") as f:
                f.write("%s %s" % (self.currentNormalNetworks[0].strip(".pairs"),self.currentNormalNetworks[-1].strip(".pairs")))
                for p in self.test_model:
                    f.write(" %f" % p)
                f.write("\n")

        stdout.write("\rCalculating p-value... \tDONE. \n" )
        stdout.flush()
        print pvals
        return pvals[index]>self.alpha



    def testNetwork(self,path,networkFileSequence):
        stdout.write("\rTesting networks... " )
        stdout.flush()
        test_loglikeL=0
        G=self.edge_count

        ng=len(networkFileSequence)
        loglikelihoodR = np.zeros(ng-1)
        logL_null = sum(self.calculateLikelihood(self.normal_model,G,g) for g in xrange(ng))

        sample_model_0 = CountModel()
        sample_model_1 = CountModel()

        for c in xrange(1,ng):

            self.calculatePriors(sample_model_0,G,ng=c,gis=range(c))
            self.calculatePriors(sample_model_1,G,ng=ng-c,gis=range(c,ng))

            S_c  =sum(self.calculateLikelihood(sample_model_0,G,g) for g in xrange(c))
            S_c+=sum(self.calculateLikelihood(sample_model_1,G,g) for g in xrange(c,ng)) - logL_null
            loglikelihoodR[c-1]=S_c

        stdout.write("\tDONE.                     \n" )
        stdout.flush()
        self.test_model=loglikelihoodR


    def fitNullModel(self,path):

        stdout.write("Fitting null model %s -- %s (%i networks)... " % (self.currentNormalNetworks[0].strip(".pairs"),self.currentNormalNetworks[-1].strip(".pairs"),len(self.currentNormalNetworks)))
        stdout.flush()

        self.fitJointNetworks(path)
        self.consensusJointNetworks(path)
        self.consensusFilename="%s_joint_s1-consensus.tree" % os.path.join(path,self.currentNormalNetworks[-1].strip(".pairs"))

        stdout.write("\tDONE.\n")
        stdout.flush()


    def createNullModel(self,path):
        self.normal_tree = self.createDendroTreeFromConsensus(self.consensusFilename)
        edgelist=[]
        G=[]
        for networkFile in self.currentNormalNetworks:
            with open(os.path.join(path,networkFile)) as f:
                edges=[edge.split() for edge in f.readlines()]
            edgelist.extend(edges)
            G.append(edges)

        self.mapDendro(self.normal_tree,G)

        self.normal_model = CountModel()

        self.calculatePriors(self.normal_model,self.edge_count,ng=len(G))

        filenames = Filenames(self.consensusFilename,path,self.currentNormalNetworks,self.namesFile)

        stdout.write("\rCreating null distribution... \t0% ")
        stdout.flush()

        nsamples = self.nsamples
        self.null_model=np.empty(nsamples)

        for ni in xrange(nsamples):
            stdout.write("\rCreating null distribution... \t%.1f%% " % (ni*100./self.nsamples))
            stdout.flush()
            loglikelihood = self.sampleNullModel(filenames,ni)
            self.null_model[ni]=loglikelihood
        stdout.write("\rCreating null distribution... \tDONE. \n")
        stdout.flush()




    def sampleNullModel(self,filenames,label):
        ng=len(filenames.currentNormalNetworks)

        #generate sample graphs
        G=self.generateGraph(self.normal_model,ng)

        #infer normal model
        normal_model = CountModel()
        self.calculatePriors(normal_model,G,ng=ng)

        loglikelihoodR = -np.inf
        logL_null = sum(self.calculateLikelihood(normal_model,G,g) for g in xrange(ng))


        max_c=-1

        sample_model_0 = CountModel()
        sample_model_1 = CountModel()

        for c in xrange(1,ng):

            self.calculatePriors(sample_model_0,G,ng=c,gis=range(c))
            self.calculatePriors(sample_model_1,G,ng=ng-c,gis=range(c,ng))

            S_c  =sum(self.calculateLikelihood(sample_model_0,G,g) for g in xrange(c))
            S_c+=sum(self.calculateLikelihood(sample_model_1,G,g) for g in xrange(c,ng)) - logL_null
            if S_c > loglikelihoodR:
                max_c=c
                loglikelihoodR = S_c

        return loglikelihoodR


    def createDendroTreeFromConsensus(self,file):
        ## Creates a dendropy Tree from the consensus file
        try:
            with open(file) as f:
                rows=f.readlines()
        except IOError:
            return None

        tree = dendropy.Tree()
        tree.is_rooted = True

        if len(rows) > 0:
            parent=self._addnodes(tree,tree.seed_node,rows,-1)

        notin=[]

        for taxon in self.taxon_namespace:
            if not tree.taxon_namespace.has_taxon_label(label=str(taxon)):
                notin.append(taxon)
                node = tree.seed_node.new_child(taxon=tree.taxon_namespace.require_taxon(label=str(taxon)),edge_length=1)

        if len(notin)>0:
            print notin

        return tree


    def mapDendro(self,tree,G):
        print "MAP"
        pdm = treemeasure.PatristicDistanceMatrix(tree)

        tns = tree.taxon_namespace

        #~ tree.edge_mrca = {}
        node_map = {}
        self.edge_count ={}
        self.mrca_count=Counter([])

        for i in xrange(self.N):
            v=self.namesLUTr[str(i)]
            nodes=pdm._mrca[tns.get_taxon(v)]
            self.mrca_count += Counter([ancestor for taxon,ancestor in nodes.iteritems()])

        n_internals=len(self.mrca_count.keys())
        node_map= dict((node,id) for id,node in enumerate(self.mrca_count.iterkeys()))

        for node in self.mrca_count.keys():
            self.mrca_count[node_map[node]] = self.mrca_count.pop(node)
            self.edge_count[node_map[node]]=np.zeros(len(G))


        for gi,edgeList in enumerate(G):
            for edge in edgeList:
                emrca = pdm.mrca(tns.get_taxon(edge[0]),tns.get_taxon(edge[1]))
                self.edge_count[node_map[emrca]][gi]+=1.
        print "DENDRO"


    def _addnodes(self,tree,parent,consensusTree,index,hrg=False):
        ##Recursively add nodes to tree. Function used by createDendroTreeFromConsensus
        row=consensusTree[index]
        if hrg:
            leafGenerator= (row.split()[i] for i in [4,5,7,8])
        else:
            leafGenerator= (val for val in row.split("\t")[3].split()[2::])
        for leaf in leafGenerator:
            nodeType = leafGenerator.next().strip("()")
            if nodeType == "D":
                node = parent.new_child(edge_length=1)
                self._addnodes(tree,node,consensusTree,int(leaf),hrg)
            else:
                node = parent.new_child(taxon=tree.taxon_namespace.require_taxon(label=leaf),edge_length=1)
        return parent


    def calculatePriors(self,tree,edge_count,ng=1,gis=None):
        #~ print "calculatePriors"
        tree.alpha = {}
        tree.beta = {}

        if gis is None:
            for node in self.mrca_count.iterkeys():
                tree.alpha[node] = 1 + edge_count[node].sum()
                tree.beta[node] = 1 + self.mrca_count[node]*ng - edge_count[node].sum()
        else:
            for node in self.mrca_count.iterkeys():
                tree.alpha[node] = 1 + edge_count[node][gis].sum()
                tree.beta[node] = 1 + self.mrca_count[node]*ng - edge_count[node][gis].sum()

                #~ print asda


    def calculatePvals(self,tree,edgeList=None,ng=1,setPriors=False):
        st=time.time()
        #~ print "Pvals"

        pdm = treemeasure.PatristicDistanceMatrix(tree)

        tns = tree.taxon_namespace

        #~ print "Init pvals",len(tree.nodes()),time.time()-st
        for node in tree.nodes():
            if node.is_internal():
                node.ei=0.0

        tree.edge_mrca={}

        if edgeList is not None:
            #~ print "Calc ei",len(edgeList),time.time()-st
            for edge in edgeList:
                emrca = pdm.mrca(tns.get_taxon(edge[0]),tns.get_taxon(edge[1]))
                emrca.ei += 1
                tree.edge_mrca[str(edge)] = emrca
                #~ self.mrca(pdm,tns.get_taxon(edge[0]),tns.get_taxon(edge[1])).ei += 1

        #~ print "Calc p",len(tree.nodes()),time.time()-st
        for node in tree.nodes():
            if node.is_internal():
                split_sizes = [len(b.leaf_nodes()) for b in node.child_nodes()]
                node.ni = np.sum(split_sizes*(np.sum(split_sizes)-split_sizes))        # number of possible links (split_sizes vector times reverse cumsum split_sizes)
                #~ node.p = (node.ei+1) / (node.ni+2) #Expected value and not actually used
                node.alpha = 1  #set priors - NB: these are overwritten for the null model
                node.beta = 1
                if setPriors:
                    node.alpha+=node.ei
                    node.beta+=(node.ni*ng-node.ei)
        #~ print "Fin",time.time()-st

    def importNames(self,namesFile):
        self.namesFile=namesFile
        with open(namesFile) as f:
            rows = f.readlines()[1:]

        self.namesLUT = dict([tuple(row.strip().split()[::-1]) for row in rows])
        self.namesLUTr = dict([tuple(row.strip().split()) for row in rows])

        self.taxon_namespace = self.namesLUT.keys()
        self.N = len(self.taxon_namespace)


    def calculateLikelihood(self,tree,edge_count,g=0,ng=1,disp=False):
        #~ print "calculating likelihood..."
        L=0.0

        for node in self.mrca_count.iterkeys():
            dL = gammaln(edge_count[node][g]+tree.alpha[node]) + gammaln(self.mrca_count[node]*ng - edge_count[node][g] + tree.beta[node])
            dL += - gammaln(self.mrca_count[node]*ng +tree.alpha[node] + tree.beta[node])
            dL += gammaln(tree.alpha[node] + tree.beta[node]) - gammaln(tree.alpha[node]) - gammaln(tree.beta[node])
            L+=dL

        #~ print "logL:",L
        return L



    def generateGraph(self,tree,ng):
        #~ print "generateGraph"
        """
        Generates edge counts - no actual graph needed.
        """
        np.random.seed()

        G=[]
        edge_count = {}
        for node in self.mrca_count.iterkeys():
            try:
                edge_count[node]=np.random.binomial(1,np.random.beta(tree.alpha[node],tree.beta[node],size=(self.mrca_count[node],ng))).sum(0)
            # deals with single possible edge case
            except AttributeError:
                edge_count[node]=np.random.binomial(1,np.random.beta(tree.alpha[node],tree.beta[node]))


        return edge_count
