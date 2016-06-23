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
import numpy as np
import os
from scipy.special import gammaln
from multiprocessing import Pool
from subprocess import call
from scipy.stats import norm


def unwrap_self_sampleNullModel(arg, **kwarg):
    return sampleNullModel(*arg, **kwarg)


class Filenames(object):
    def __init__(self,consensusFilename,path,currentNormalNetworks,namesFile):
        self.consensusFilename = consensusFilename
        self.path = path
        self.currentNormalNetworks = currentNormalNetworks
        self.namesFile = namesFile
    


def sampleNullModel(filenames,label):
    ad = anomalyDetection()
    ad.importNames(filenames.namesFile)
    ng=len(filenames.currentNormalNetworks)
    normal_model = ad.createDendroTreeFromConsensus(filenames.consensusFilename)####
    edgelist=[]
    for networkFile in filenames.currentNormalNetworks:###
        with open(os.path.join(filenames.path,networkFile)) as f:
            edgelist.extend([edge.split() for edge in f.readlines()])
    ad.calculatePvals(normal_model,edgelist,ng,setPriors=True)
    
    #generate sample graphs
    G=[]
    for gi in xrange(ng):
        G.append(ad.generateGraph(normal_model))
    
    ad.calculatePvals(normal_model,[edge for g in G for edge in g],ng,setPriors=True)
    
    loglikelihoodR = -np.inf
    logL_null = sum(ad.calculateLikelihood(normal_model,g) for g in G)
    #~ logL_null = ad.calculateLikelihood(normal_model,[edge for g in G for edge in g])
    max_c=-1
    for c in xrange(1,ng):
        sample_model_0 = ad.createDendroTreeFromConsensus(filenames.consensusFilename)
        sample_model_1 = ad.createDendroTreeFromConsensus(filenames.consensusFilename)
        G_0=[edge for g in G[:c] for edge in g]
        G_1=[edge for g in G[c:] for edge in g]
        ad.calculatePvals(sample_model_0,G_0,c,setPriors=True)
        ad.calculatePvals(sample_model_1,G_1,ng-c,setPriors=True)
        S_c  =sum(ad.calculateLikelihood(sample_model_0,g) for g in G[:c])
        #~ S_c  =ad.calculateLikelihood(sample_model_0,G_0)
        S_c+=sum(ad.calculateLikelihood(sample_model_1,g) for g in G[c:]) - logL_null
        #~ S_c+=ad.calculateLikelihood(sample_model_1,G_1) - logL_null
        if S_c > loglikelihoodR:
            max_c=c
            loglikelihoodR = S_c
    
    return loglikelihoodR
    

class anomalyDetection(object):
    
    def __init__(self):
        self.taxon_namespace = None
        self.namesLUT = None
        self.nsamples = 1000
        self.window=None
    
    
    def detectAnomaliesInSequence(self,networkFileSequence,namesFile,window,path="."):
        self.namesFile=os.path.join(path,namesFile)
        self.window=window
        self.importNames(self.namesFile)
        self.alpha=0.95
        t=0
        normal=True
        while t<len(networkFileSequence)-window:
            #~ print "test", test,"window at t=", t
            startindx=t
            endindx=t+window
            
            paramStr="-n %s -m %i" % (self.namesFile,window)

            for gi in xrange(startindx,endindx):
                paramStr+=" %s" % (os.path.join(path,networkFileSequence[gi]))
            
            paramStr+=" > /dev/null"
            
            call("./fitHRG_GPL_Bayes/fitHRG -t joint %s" % paramStr, shell=True)
            call("./consensusHRG_GPL_Bayes/consensusHRG %s -f %s_joint_best-dendro.hrg" % (paramStr,os.path.join(path,networkFileSequence[gi].strip(".pairs"))), shell=True)
            
            #then run null model from ad.
            self.consensusFilename="%s_joint_s1-consensus.tree" % os.path.join(path,networkFileSequence[gi].strip(".pairs"))
            self.currentNormalNetworks=[networkFileSequence[gi] for gi in xrange(startindx,endindx)]
            self.createNullModel(path)
            self.testNetwork(path,self.currentNormalNetworks)
            
            changeDetected = self.calculateP_values()
            if changeDetected:
                tau=t+window
                t+=np.argmax(self.pvals)
            t+=1
            #~ normal =  (t < experimentLength) and not changeDetected
        #~ return tau, r, time.time()-stime
    
    
    
    
    def calculateP_values(self):
        stdout.write("\rCalculating p-value..." )
        stdout.flush()
        #~ with open("approx/null_%s_%i.model" % (self.currentNormalNetworks[-1].strip(".pairs"),method)) as f:
            #~ null_model = [row.strip() for row in f.readlines()]
        
        #~ with open("approx/test_%s_%i.model" % (self.currentNormalNetworks[-1].strip(".pairs"),method)) as f:
            #~ test_model= [row.strip() for row in f.readlines()]
        
        #~ null_model = np.float64(null_model)
        #~ test_model = np.float64(test_model)
        
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
            with open("results%i.txt" % self.window, "a") as f:
                f.write("%s %s" % (self.currentNormalNetworks[0].strip(".pairs"),self.currentNormalNetworks[-1].strip(".pairs")))
                for p in pvals:
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
        G=[]
        for networkFile in networkFileSequence:
            with open(os.path.join(path,networkFile)) as f:
                G.append([edge.split() for edge in f.readlines()])
        
        ng=len(networkFileSequence)
        loglikelihoodR = np.empty(ng-1)
        logL_null = sum(self.calculateLikelihood(self.normal_model,g) for g in G)
        #~ logL_null = self.calculateLikelihood(self.normal_model,[edge for g in G for edge in g])
        #~ max_c=-1
        for c in xrange(1,ng):
            sample_model_0 = self.createDendroTreeFromConsensus(self.consensusFilename)
            sample_model_1 = self.createDendroTreeFromConsensus(self.consensusFilename)
            G_0=[edge for g in G[:c] for edge in g]
            G_1=[edge for g in G[c:] for edge in g]
            self.calculatePvals(sample_model_0,G_0,c,setPriors=True)
            self.calculatePvals(sample_model_1,G_1,ng-c,setPriors=True)
            S_c  =sum(self.calculateLikelihood(sample_model_0,g) for g in G[:c])
            #~ S_c  =self.calculateLikelihood(sample_model_0,G_0)
            S_c+=sum(self.calculateLikelihood(sample_model_1,g) for g in G[c:]) - logL_null
            #~ S_c+=self.calculateLikelihood(sample_model_1,G_1) - logL_null
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
    
    def createNullBaselineModel(self,sequence):
        ng=len(sequence)
        
        null_mean = np.mean(sequence)
        null_std = np.std(sequence)
        #~ print "null", null_std, null_mean
        
        loglikelihoods = np.zeros(self.nsamples)-np.inf
        for i in xrange(self.nsamples):
            sample = norm.rvs(null_mean,null_std,size=ng)
            sample_std = sample.std()
            logL_null = np.sum(norm.logpdf(sample,sample.mean(),sample_std))
            
            max_c=-1
            for c in xrange(1,ng):
                S_c = np.sum(norm.logpdf(sample[:c],sample[:c].mean(),sample_std))
                S_c += np.sum(norm.logpdf(sample[c:],sample[c:].mean(),sample_std)) - logL_null
                if S_c > loglikelihoods[i]:
                    max_c=c
                    loglikelihoods[i] = S_c
        
        self.null_model = loglikelihoods
        #~ loglikelihoods = np.append(np.sort(loglikelihoods),np.inf)
        
        logL_null = np.sum(norm.logpdf(sequence,null_mean,null_std))
        
        test_loglikeL = np.zeros(ng-1)
        max_c=-1
        for c in xrange(1,ng):
            test_loglikeL[c-1] = np.sum(norm.logpdf(sequence[:c],sequence[:c].mean(),null_std))
            test_loglikeL[c-1] += np.sum(norm.logpdf(sequence[c:],sequence[c:].mean(),null_std)) - logL_null
        
        self.test_model = test_loglikeL
        
        
    
    def createNullModel(self,path):
        self.normal_model = self.createDendroTreeFromConsensus(self.consensusFilename)
        edgelist=[]
        for networkFile in self.currentNormalNetworks:
            with open(os.path.join(path,networkFile)) as f:
                edgelist.extend([edge.split() for edge in f.readlines()])
        self.calculatePvals(self.normal_model,edgelist,len(self.currentNormalNetworks),setPriors=True)
        
        
        filenames = Filenames(self.consensusFilename,path,self.currentNormalNetworks,self.namesFile)
        
        with open("nProcessors.txt") as f:
            nProcessors=int(f.readlines()[0])
        
        stdout.write("\rCreating null distribution... \t0% ")
        stdout.flush()
        
        nsamples = self.nsamples
        self.null_model=np.empty(nsamples)
        
        if nProcessors>1:
        
            sampler = Pool(processes=nProcessors)
            
            already_done = 0
            
            nsamples-=already_done
            if nsamples>0:
                loglikelihoods = sampler.imap(unwrap_self_sampleNullModel, zip([filenames]*nsamples, range(nsamples)),chunksize=(nsamples/(nProcessors*10)))
                #~ loglikelihoods = sampler.imap(unwrap_self_sampleNullModel, [filenames]*nsamples,chunksize=(nsamples/(nProcessors)))
                
                for i,loglikelihood in enumerate(loglikelihoods):
                    stdout.write("\rCreating null distribution... \t%i%% " % ((i+1+already_done)*100/self.nsamples))
                    stdout.flush()
                    if loglikelihood is not None:
                        self.null_model[i]=loglikelihood
                sampler.close()
                sampler.join()
                
        else:
            for ni in xrange(nsamples):
                loglikelihood = unwrap_self_sampleNullModel((filenames,ni))
                self.null_model[ni]=loglikelihood
        stdout.write("\rCreating null distribution... \tDONE. \n")
        stdout.flush()
    
    
    
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
                node = tree.seed_node.new_child(taxon=tree.taxon_namespace.require_taxon(label=str(taxon)))
        
        if len(notin)>0:
            print notin
        
        #~ tree.reroot_at_node(tree.seed_node)
        
        return tree
        


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
                node = parent.new_child()
                self._addnodes(tree,node,consensusTree,int(leaf),hrg)
            else:
                node = parent.new_child(taxon=tree.taxon_namespace.require_taxon(label=leaf))
        return parent
    
    
    def calculatePvals(self,tree,edgeList=None,ng=1,setPriors=False):
        for node in tree.nodes():
            if node.is_internal():
                node.ei=0.0
        
        if edgeList is not None:
            for edge in edgeList:
                tree.mrca(taxon_labels=[edge[0],edge[1]]).ei += 1
        
        for node in tree.nodes():
            if node.is_internal():
                split_sizes = [len(b.leaf_nodes()) for b in node.child_nodes()]
                node.ni = np.sum(split_sizes*(np.sum(split_sizes)-split_sizes))        # number of possible links (split_sizes vector times reverse cumsum split_sizes)
                node.p = (node.ei+1) / (node.ni+2) #Expected value and not actually used
                node.alpha = 1  #set priors - NB: these are overwritten for the null model
                node.beta = 1
                if setPriors:
                    node.alpha+=node.ei
                    node.beta+=(node.ni*ng-node.ei)
                    
    
    def importNames(self,namesFile):
        self.namesFile=namesFile
        with open(namesFile) as f:
            rows = f.readlines()[1:]
        
        self.namesLUT = dict([tuple(row.strip().split()[::-1]) for row in rows])
        self.namesLUTr = dict([tuple(row.strip().split()) for row in rows])
        
        self.taxon_namespace = self.namesLUT.keys()
        self.N = len(self.taxon_namespace)
    
    
    def calculateLikelihood(self,tree,edgeList,ng=1,disp=False):
        #~ print "calculating likelihood..."
        for node in tree.nodes():
            if node.is_internal():
                node.ei=0.0
        
        
        for edge in edgeList:
            tree.mrca(taxon_labels=[edge[0],edge[1]]).ei += 1
        
        L = 0.0
        for node in tree.nodes():
            if node.is_internal():
                dL = gammaln(node.ei+node.alpha) + gammaln(node.ni*ng - node.ei + node.beta) - gammaln(node.ni*ng + node.alpha + node.beta)
                dL += gammaln(node.alpha+node.beta) - gammaln(node.alpha) - gammaln(node.beta)
                L+=dL
                if disp:
                    print "a=",node.alpha, "b=",node.beta,"ei=",node.ei,"ni=",node.ni,"LL=",dL
        if disp:
            print L
        return L
    
    
    
    def generateGraph(self,tree):
        np.random.seed()
        edgeList = []
        for i in range(self.N):
            for j in range((i+1),self.N):
                v=self.namesLUTr[str(i)]
                u=self.namesLUTr[str(j)]
                node=tree.mrca(taxon_labels=[v,u])
                p = np.random.beta(node.alpha,node.beta)
                if np.random.rand() < p:
                    edgeList.append([v,u])
                    edgeList.append([u,v])
        return edgeList
