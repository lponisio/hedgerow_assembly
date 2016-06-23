"""activelearn.py - module to run active discovery of network roles
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

import supervisedBM 
from random import sample
import numpy as np
import os


def expInstance(args):
    np.random.seed()
    networkfile,classlist,iters,no = args
    exp=os.path.split(networkfile)[1].split(".")[0]
    print exp
    #~ print asdsad
    #~ networkfile="../data/%s_edges.txt" % exp
    #~ classlist ="../data/%s_labels.txt" % exp
    inG = supervisedBM.nx.read_adjlist(networkfile, create_using=supervisedBM.nx.DiGraph())
    
    with open(classlist) as f:
        labels = f.readlines()
    class_lut=dict((label.split()[0],label.split()[1].strip()) for label in labels)
    
    #select one of each class
    classes=np.unique(class_lut.values())
    print classes
    for c in classes:
        nodes_of_class_c = [k for k,v in class_lut.iteritems() if v==c]
        v=nodes_of_class_c[int(np.floor(np.random.rand()*len(nodes_of_class_c)))]
        inG.node[v]["label"]=class_lut[v]
    
    
    #~ sb = supervisedBM.MixedMembership(inG)
    sb = supervisedBM.MaxBM(inG,K=(2*len(classes)))
    sb.id=no
    
    for i in range(iters):
        sb.infer()
        
        #X-validation step
        maxacc=0
        maxr=64
        print "*************"
        training_nodes = dict((node,nodeData["cindx"]) for node,nodeData in sb.G.nodes_iter(data=True) if nodeData["cindx"] is not None)
        if len(training_nodes)>5:
            for r in [64,49,36,25,16,9,4,1,81,100]:
                xv_count=0
                sb.reg=r
                accr=0
                tn=np.floor(len(training_nodes)/5.)
                tn+=np.random.rand()<(len(training_nodes)/5. - tn)
                print tn
                for xv in range(5):
                    tnodes=sample(training_nodes.keys(),int(tn))
                    test_labels=[]
                    for node in training_nodes.keys():
                        if node in tnodes:
                            sb.G.node[node]["cindx"] = training_nodes[node]
                        else:
                            sb.G.node[node]["cindx"] = None
                            test_labels.append(node+"\t"+class_lut[node])
                    try:
                        sb.update_eta()
                        sb.predict()
                        #~ print "labels",len(test_labels)
                        accr+=validate(sb.G,test_labels)
                        xv_count+=1
                        #~ print "ac",accr,xv_count
                    except:
                        pass
                    #~ print accr
                if xv_count>0:
                    #~ print accr, xv_count
                    accr/=xv_count
                    print "reg",r,accr
                    if accr>maxacc:
                        maxacc=accr
                        maxr=r
        #return training set
        for node in training_nodes.keys():
            sb.G.node[node]["cindx"] = training_nodes[node]
        
        print "REG=",maxr
        sb.reg=maxr
        sb.update_eta()
        sb.predict()
        acc=validate(sb.G,labels)
        v=H(sb.G)
        sb.G.node[v]["label"]=class_lut[v]
        sb.setupClasses()
        acc=validate(sb.G,labels)
        with open("result_%s_%i.txt" % (exp,no),"a") as f:
            f.write("%i:%s %s %f %f; \n" % (i,v,class_lut[v],sb.reg,acc))
    with open("result_%s_%i.txt" % (exp,no),"a") as f:
        f.write("\n")
    return sb


def H(inG):
    Hmax=0
    Mmin=np.inf
    for node,v in inG.nodes_iter(data=True):
        if v["cindx"] is None:
            try:
                p=v["prob"]
                Hv=-np.sum(supervisedBM.xlogy(p,p))
                #~ print Hv,p
                if Hv > Hmax:
                    Hmax=Hv
                    vmax=node
            except:
                if v["margin"] < Mmin:
                    Mmin= v["margin"] 
                    vmax=node
                #~ print "new max",node,Hmax
    return vmax


def validate(G,labels):
    tt=0
    yt=0
    for label in labels:
        v,l = label.split()
        #~ print l.strip(),G.node[v]["prediction"],G.node[v]["cindx"] is None
        if G.node[v]["cindx"] is None:
            tt+=1.
            yt+=l.strip()==G.node[v]["prediction"]
    print tt,yt,yt/tt
    return yt/tt


if __name__=="__main__":
    print "supervisedBM v0.6"
    import argparse
    parser = argparse.ArgumentParser()
    networkfiletype = parser.add_mutually_exclusive_group()
    networkfiletype.add_argument("-a", "--adjlist", action="store_true", help="adjacency list network format (default)")
    #~ networkfiletype.add_argument("-q", "--quiet", action="store_true") # we can add functionality to use other networkx input methods here.
    parser.add_argument("networkfile", help="Input network file")
    parser.add_argument("classlist", help="File of class labels list")
    parser.add_argument("iterations", help="Number of active learning iterations")
    #~ parser.add_argument("-o","--outputFile", help="Output File")
    #~ parser.add_argument("-p","--cen", help="Censor prob")
    args = parser.parse_args()
    
    expInstance((args.networkfile,args.classlist,int(args.iterations),0))