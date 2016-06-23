"""supervisedBM.py - supervised blockmodel module for classification of network nodes
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


import networkx as nx
#~ from prettyplotlib import plt
#~ import prettyplotlib as ppl
#~ from sklearn import manifold
import numpy as np
from scipy.special import gammaln
from scipy.cluster.vq import kmeans,whiten,vq
from operator import mul
from scipy.optimize import fmin_cg
from sys import stdout
#~ import _pcolormesh
import json
from networkx.readwrite import json_graph

def num(s):
    try:
        return int(s)
    except ValueError:
        try:
            return float(s)
        except ValueError:
            return s

try:
    import numexpr as ne
    def xlogy(x, y):
        return ne.evaluate("log(y+1e-200)*x")
except:
    #safe log function
    def xlogy(x, y):
        return x*np.log(y+1e-200)


class Supervised_Blockmodel(object):
    def __init__(self, G=None,K=None):
        self.id=0
        if G is not None:
            self.G=G
        params={}
        with open("params.txt") as file:
            for line in file.read().splitlines():
                items = line.split(":")
                key, values = items[0], items[1].strip()
                params[key] = num(values)
        self.maxiterations=params["maxiterations"]
        self.miniterations=2
        self.convtol = 1e-6
        self.fp_tol = 1e-5
        if K is None:
            self.K = params["K"]
        else:
            self.K = K
        self.I = len(self.G.edges())
        self.alpha = np.ones((self.K,self.K))*params["alpha"]#+(np.eye(self.K)+np.eye(self.K,self.K,1))*self.I/(self.K*2.-1)
        self.beta = params["beta"]
        self.disp = bool(params["disp"])
        print self.disp
        self.reg = params["reg"]
        self.N = len(self.G)
        self.initialiseZ()
        self.setupClasses()
        
    def setupClasses(self):
        print "****************************************************************************"
        #set lookup for classes
        self.cindx_lut={}
        self.cindx_lutr={}
        cindx=0
        for name,node in self.G.node.iteritems():
            #~ print node
            if not node.has_key("label"):
                print name,node
                node["label"]=""
            if node["label"]=="":
                node["cindx"] = None
            else:
                if not self.cindx_lut.has_key(node["label"]):
                    self.cindx_lut[node["label"]]=cindx
                    self.cindx_lutr[cindx]=node["label"]
                    cindx+=1
                node["cindx"] = self.cindx_lut[node["label"]]
        self.C = cindx
    
    def infer(self):
        t=0
        self.ELBO = 1
        self.ELBOdiff = 1
        lastt = self.maxiterations
        self.fopt=0
        self.eta_updates=self.N/float(self.I)
        self.eta_counts=self.I
        edgelist = self.G.edges()
        order = np.arange(self.I)
        #~ print self.G.edges(),self.I
        self.update_eta()
        while t < self.maxiterations:
            np.random.shuffle(order)
            self.predict()
            for i,eidx in enumerate(order) :
                if self.disp:
                    stdout.write("********Iteration %i : %i%% \r" %(t,i*100/self.I))
                    stdout.flush()
                e=edgelist[eidx]
                self.updateZ(e)
            self.update_eta()
            self.oldELBO = self.ELBO
            self.calculateELBO(t)
            t += 1
            if ((self.ELBOdiff < self.convtol))*(t>self.miniterations):
                lastt = t
                t = self.maxiterations
        self.predict()


class MixedMembership(Supervised_Blockmodel):
    def __init__(self, G=None,K=None):
        print len(G)
        super(MixedMembership, self).__init__(G,K)
        self.etamaxiter = 1000
        self.G.graph["sum_nv"]=np.zeros(self.K)
        for node,nodeData in self.G.nodes_iter(data=True):
            fv= np.array(sum(v["Zs"] for v in self.G.succ[node].itervalues()))
            gv= np.array(sum(v["Zr"] for v in self.G.pred[node].itervalues()))
            nodeData["nv"]=fv+gv
            self.G.graph["sum_nv"]+=nodeData["nv"]
        self.G.graph["sum_nv"]+=self.N*self.beta
        self.eta = np.random.randn(self.K,self.C)*self.reg + np.eye(self.K,self.C)
        self.getnormh()
    
    
    def initialiseZ_kmeans(self):
        nodes=self.G.nodes()
        node_lut=dict(((n,v) for v,n in enumerate(nodes)))
        X1=np.eye(self.N)
        X2=np.eye(self.N)
        for v,node in enumerate(nodes):
            X1[v,[node_lut[sv] for sv in self.G.successors(node)]]=1
            X2[v,[node_lut[pv] for pv in self.G.predecessors(node)]]=1
        X=np.append(X1,X2,1)+(np.random.rand(self.N,self.N*2)*1e-100)
        W=whiten(X)
        cm,zv=kmeans(W,self.K,iter=100)
        zv,d_=vq(W,cm)
        zinit=dict(zip(nodes,zv))
        self.G.graph["c"]=np.zeros(np.shape(self.alpha))
        for link in self.G.edges_iter(data=True):
            edge=link[2]
            Z=np.random.mtrand.dirichlet(self.alpha.flatten()).reshape(self.K,self.K)
            Z[zinit[link[0]],:]+=10
            Z[:,zinit[link[1]]]+=10
            Z/=np.sum(Z)
            edge["Z"]=Z
            edge["Zs"]=np.sum(Z,1)
            edge["Zr"]=np.sum(Z,0)
            self.G.graph["c"]+=Z
        
        
    
    def initialiseZ(self):
        K2 = self.K*self.K
        self.G.graph["c"]=np.zeros(np.shape(self.alpha))
        for link in self.G.edges_iter(data=True):
            edge=link[2]
            edge["Z"]=np.random.mtrand.dirichlet(self.alpha.flatten()).reshape(self.K,self.K)
            Z=edge["Z"]
            edge["Zs"]=np.sum(Z,1)
            edge["Zr"]=np.sum(Z,0)
            self.G.graph["c"]+=Z
    
    
    def calculateELBO(self,t):
        self.ELBO = -self.fopt + np.sum(self.eta*self.eta*self.reg/2.)
        self.ELBO += gammaln(np.sum(self.alpha)) - np.sum(gammaln(self.alpha))
        self.ELBO += np.sum(gammaln(self.G.graph["c"]+self.alpha)) - gammaln(np.sum(self.G.graph["c"]+self.alpha))
        
        self.ELBO += self.K*(gammaln(self.beta*self.N) - gammaln(self.beta)*self.K*self.N)
        self.ELBO += sum(sum(gammaln(v["nv"]+self.beta) for node,v in self.G.nodes_iter(data=True)))
        self.ELBO -= np.sum(gammaln(self.G.graph["sum_nv"]+self.N*self.beta))
        
        
        self.ELBO -= np.sum(sum(xlogy(edge[2]["Z"],edge[2]["Z"]) for edge in self.G.edges_iter(data=True)))
        self.ELBOdiff = np.abs((self.ELBO-self.oldELBO)/self.ELBO)
        #~ self.ELBOdiff = (self.ELBO-self.oldELBO)
        update="***"
        if self.disp:
            #~ print " f1=%.2f(%.2f) acc=%.2f(%.2f)  t,self.f1_tr,self.f1_te,self.p_train, self.p_test, )
            print "%i:ELBO = %f  ELBO diff = %f" % (t,self.ELBO, self.ELBOdiff,)
    
    
    def updateZ(self,link):
        
        #remove counts including current node v
        G=self.G
        s=G.node[link[0]]
        r=G.node[link[1]]
        edge=G.edge[link[0]][link[1]]
        #~ print edge["Zr"]
        zs=edge["Zs"]
        zr=edge["Zr"]
        z=edge["Z"]
        
        G.graph["c"] -= z
        
        s["nv"] -= zs
        r["nv"] -= zr
        G.graph["sum_nv"]-=zs
        G.graph["sum_nv"]-=zr
        
        z_v = (G.graph["c"]+self.alpha)*np.dot((s["nv"]+self.beta)[:,np.newaxis],(r["nv"] +self.beta)[np.newaxis,:]) 
        z_v /= (np.dot(G.graph["sum_nv"][:,np.newaxis],G.graph["sum_nv"][np.newaxis,:]) + np.eye(self.K))
        
        
        ys = s["cindx"]
        yr = r["cindx"]
        
        
        if ys is not None:
            expetaNs = self.exp_etaNv[self.G.degree()[link[0]]]
            expetaNsy = expetaNs[:,ys]
            z_v*=expetaNsy[:,np.newaxis]
        if yr is not None:
            expetaNr = self.exp_etaNv[self.G.degree()[link[1]]]
            expetaNry = expetaNr[:,yr]
            z_v*=expetaNry[np.newaxis,:]
        
        
        
        fp = 1
        newZi = np.copy(z)
        hnorm_si = np.copy(s["hnorm"])
        hnorm_ri = np.copy(r["hnorm"])
        while fp > self.fp_tol: 
            oldZi = np.copy(newZi)
            newZi = np.copy(z_v)
            if ys is not None:
                hsdotz = np.sum(hnorm_si)
                hnorm_si /= np.dot(zs,expetaNs)
                hs = np.sum(expetaNs * hnorm_si, 1)
                newZi *= np.exp(-(hs/hsdotz)[:,np.newaxis])
            if yr is not None:
                hrdotz = np.sum(hnorm_ri)
                hnorm_ri /= np.dot(zr,expetaNr)
                hr = np.sum(expetaNr * hnorm_ri, 1)
                newZi *= np.exp(-(hr/hrdotz)[np.newaxis,:])
            
            newZi /= np.sum(newZi)
            
            zs = np.sum(newZi,1)
            zr = np.sum(newZi,0)
            
            if ys is not None:
                hnorm_si *= np.dot(zs,expetaNs)
            if yr is not None:
                hnorm_ri *= np.dot(zr,expetaNr)
            
            fp = np.sum(np.abs(oldZi - newZi))/(self.K*self.K)
        
        s["hnorm"] = hnorm_si
        r["hnorm"] = hnorm_ri
        
        edge["Zs"]=zs
        edge["Zr"]=zr
        edge["Z"]=newZi
        self.G.graph["sum_nv"]+=zr
        self.G.graph["sum_nv"]+=zs
        s["nv"] += zs
        r["nv"] += zr
        G.graph["c"] += newZi
    
    def getnormh(self):
        self.hnorm = {}
        degrees=self.G.degree()
        uniqueDegrees=np.unique(degrees.values())
        #cache values for e^(\eta/d)
        exp_etaNv = dict((d,np.exp(self.eta/d)) for d in uniqueDegrees)
        C=self.C
        Cones=np.ones(C)
        for node,nodeData in self.G.nodes_iter(data=True):
            expetaN = exp_etaNv[degrees[node]]
            nodeData["hnorm"] = Cones
            nodeData["hnorm"] = reduce(mul,(np.dot(v["Zs"],expetaN) for v in self.G.succ[node].itervalues()),nodeData["hnorm"]) 
            nodeData["hnorm"] = reduce(mul,(np.dot(v["Zr"],expetaN) for v in self.G.pred[node].itervalues()),nodeData["hnorm"])
        self.exp_etaNv = exp_etaNv
        
    def update_eta(self):
        
        eta_prime = np.zeros((self.K,self.C))
        for node,nodeData in self.G.nodes_iter(data=True):
            d = self.G.degree()[node]
            #~ nodeData["Z_hat"] = sum(v["Zs"] for v in self.G.succ[node].itervalues()) / d    ##change:update z_hat with update z?
            #~ nodeData["Z_hat"] += sum(v["Zr"] for v in self.G.pred[node].itervalues()) / d
            nodeData["Z_hat"] = nodeData["nv"]/d
            if nodeData["cindx"] is not None:
                eta_prime[:,nodeData["cindx"]] += nodeData["Z_hat"]
        
        self.eta_prime=eta_prime
        
        eta = self.eta.reshape(self.eta.size)
        fopt_b = self.L_eta(eta)
        eta, fopt, func_calls, grad_calls, warnflag = fmin_cg(self.L_eta, eta, self.L_eta_prime, maxiter=self.etamaxiter, disp=0,full_output=1)
        
        #only accept changes that improve ELBO
        if (fopt_b-fopt)>self.fp_tol or self.fopt<=0:
            if self.disp:
                print "ETA ACCEPT",func_calls, grad_calls, warnflag, fopt_b, fopt, fopt_b-fopt
            self.eta = eta.reshape(self.K, self.C)
            self.fopt = fopt
            self.getnormh()
        else:
            if self.disp:
                print "eta reject",func_calls, grad_calls, warnflag, fopt_b, fopt, fopt_b-fopt
            self.fopt=fopt_b
        #~ self.etamaxiter = 1
        #~ self.etamaxiter = np.min([30,self.etamaxiter+1])

    
    def L_eta(self, eta) :
        
        L_eta = 0
        eta = eta.reshape(self.K, self.C)
        
        degrees=self.G.degree()
        uniqueDegrees=np.unique(degrees.values())
        #cache values for e^(\eta/d)
        exp_etaNv = dict((d,np.exp(eta/d)) for d in uniqueDegrees) 
        Cones=np.ones(self.C)
        for node,nodeData in self.G.nodes_iter(data=True):
            y = nodeData["cindx"]
            if y is not None:
                expetaN = exp_etaNv[degrees[node]]
                L_eta += np.dot(nodeData["Z_hat"],eta[:,y])
                sCones=Cones
                rCones=Cones
                L_eta += -xlogy(1.0,np.sum(reduce(mul,(np.dot(v["Zs"],expetaN) for v in self.G.succ[node].itervalues()),sCones) * reduce(mul,(np.dot(v["Zr"],expetaN) for v in self.G.pred[node].itervalues()),rCones) ))
        #~ print "L_eta", L_eta
        return -(L_eta - np.sum(eta*eta*self.reg/2.))
        

    def L_eta_prime(self, eta) :
        
        K = self.K
        C = self.C
        
        eta = eta.reshape(K, C)
        L_eta_prime = np.zeros((K,C))
        
        degrees=self.G.degree()
        uniqueDegrees=np.unique(degrees.values())
        exp_etaNv = dict((d,np.exp(eta/d)) for d in uniqueDegrees)
        
        
        for node,nodeData in self.G.nodes_iter(data=True):
            y = nodeData["cindx"]
            if y is not None:
                d = degrees[node]
                expetaN = exp_etaNv[d]
                
            
                z1s = np.array([v["Zs"] for v in self.G.succ[node].itervalues()])
                z2r = np.array([v["Zr"] for v in self.G.pred[node].itervalues()])
                
                #c_level dims: 1 x C
                #catch instances where there is no in- or out- degree
                try:
                    c_level = np.prod(np.dot(z1s,expetaN), 0) #u
                    #this_k_levels dims: ns x k2 x C
                    this_k_levels = z1s[:,:,np.newaxis] * expetaN#w
                    this_k_levels /= (d * np.sum(this_k_levels,1)[:,np.newaxis,:])#x
                except ValueError:
                    this_k_levels=0
                    try:
                        c_level = np.prod(np.dot(z2r,expetaN), 0)
                        this_k_levelr = z2r[:,:,np.newaxis] * expetaN#w
                        this_k_levelr /= (d * np.sum(this_k_levelr,1)[:,np.newaxis,:])#x
                    except ValueError:
                        c_level = np.zeros(C)
                        this_k_levelr=0
                else:
                    try:
                        c_level*= np.prod(np.dot(z2r,expetaN), 0)
                        this_k_levelr = z2r[:,:,np.newaxis] * expetaN#w
                        this_k_levelr /= (d * np.sum(this_k_levelr,1)[:,np.newaxis,:])#x
                    except ValueError:
                        this_k_levelr=0
                c_level /= np.sum(c_level)#v
            
                
                
                L_eta_prime += -c_level.reshape(1,C) * (np.sum(this_k_levels,0) + np.sum(this_k_levelr,0))
            
        L_eta_prime -= self.reg*eta
        return -(self.eta_prime+L_eta_prime).reshape(L_eta_prime.size)
        
    def predict(self):
        for node,v in self.G.nodes_iter(data=True):
            v["prob"]=np.exp(np.dot(v["Z_hat"], self.eta))
            v["prob"]/=np.sum(v["prob"])
            v["prediction"]=self.cindx_lutr[np.argmax(np.exp(np.dot(v["Z_hat"], self.eta)))]
    
    
    
class MaxBM(MixedMembership):
    
    def getnormh(self):
        pass
    
    def __init__(self, G,K=None):
        super(MaxBM, self).__init__(G,K)
        from subprocess import call
        self.call = call
        self.svm = 0
        
    
    def updateZ(self,link):
        
        #remove counts including current node v
        G=self.G
        s=G.node[link[0]]
        r=G.node[link[1]]
        edge=G.edge[link[0]][link[1]]
        #~ print edge["Zr"]
        zs=edge["Zs"]
        zr=edge["Zr"]
        z=edge["Z"]
        
        G.graph["c"] -= z
        
        s["nv"] -= zs
        r["nv"] -= zr
        G.graph["sum_nv"]-=zs
        G.graph["sum_nv"]-=zr
        eta=self.eta
        
        z_v = (G.graph["c"]+self.alpha)*np.dot((s["nv"]+self.beta)[:,np.newaxis],(r["nv"] +self.beta)[np.newaxis,:]) 
        z_v /= (np.dot(G.graph["sum_nv"][:,np.newaxis],G.graph["sum_nv"][np.newaxis,:]) + np.eye(self.K))
        
        
        ys = s["cindx"]
        yr = r["cindx"]
        
        #~ print "nv",s["nv"],r["nv"]
        
        if ys is not None:
            mu_idx = self.mu[self.tnode_dict[link[0]],:]
            zt = (mu_idx[np.newaxis,:,np.newaxis]*((eta.T)[ys,:] - (eta.T))).sum(1)/self.G.degree()[link[0]]
            z_v *= np.exp(zt-zt.max()).T
        if yr is not None:
            mu_idx = self.mu[self.tnode_dict[link[1]],:]
            zt = (mu_idx[np.newaxis,:,np.newaxis]*((eta.T)[yr,:] - (eta.T))).sum(1)/self.G.degree()[link[0]]
            z_v *= np.exp(zt-zt.max())
        
        
            
        newZi = z_v/np.sum(z_v)
        
        zs = np.sum(newZi,1)
        zr = np.sum(newZi,0)
        
        edge["Zs"]=zs
        edge["Zr"]=zr
        edge["Z"]=newZi
        self.G.graph["sum_nv"]+=zr
        self.G.graph["sum_nv"]+=zs
        s["nv"] += zs
        r["nv"] += zr
        G.graph["c"] += newZi
    
    
    
    
    def update_eta(self):
        
        v_train=0
        self.tnode_dict={}
        with open("z%i.dat" % self.id,"w") as f:
            for node,nodeData in self.G.nodes_iter(data=True):
                d = self.G.degree()[node]
                nodeData["Z_hat"] = nodeData["nv"]/d
                if nodeData["cindx"] is not None:
                    self.tnode_dict[node]=v_train
                    v_train+=1
                
                    f.write("%i" % (nodeData["cindx"]+1))
                    for k in xrange(self.K):
                        f.write(" %i:%f" % ((k+1),nodeData["Z_hat"][k]))
                    f.write("\n")
        
        self.call(["./svm_multiclass_learn", "-c", "%f" % self.reg, "-w", "0","-p","2", "-o", "1", "-v", "-1", "z%i.dat" % self.id, "model%i" % self.id])
        
        
        with open("model%i_w.dat" % self.id) as f:
            w = f.readlines()[0].strip(" \n")
            w = np.float64([wi for wi in w.split()])
        
        with open("model%i_mu.dat" % self.id) as f:
            alphas = [[ai for ai in a.strip("\n").split()] for a in f.readlines()]
        
        mu = np.zeros(self.C*len(self.tnode_dict.keys()))
        for a in alphas:
            #~ mu[int(a[1])] = float(a[2])
            mu[int(a[1])] = float(a[2])
            
        #~ print "w",w,"mu",mu,alphas
        
        w = w[1:]
        
        self.eta = w.reshape(self.C, self.K).T
        self.mu = mu.reshape(len(self.tnode_dict.keys()),self.C)
    

    def predict(self) :
        nodeOrder=[]
        with open("z%i.dat" % self.id,"w") as f:
            for node,nodeData in self.G.nodes_iter(data=True):
                nodeOrder.append(node)
                d = self.G.degree()[node]
                #~ nodeData["Z_hat"] = sum(v["Zs"] for v in self.G.succ[node].itervalues()) / d    ##change:update z_hat with update z?
                #~ nodeData["Z_hat"] += sum(v["Zr"] for v in self.G.pred[node].itervalues()) / d
                nodeData["Z_hat"] = nodeData["nv"]/d
                
                f.write("1")# % (nodeData["cindx"]+1))
                for k in xrange(self.K):
                    f.write(" %i:%f" % ((k+1),nodeData["Z_hat"][k]))
                f.write("\n")
        self.call(["./svm_multiclass_classify", "-v", "-1", "z%i.dat" % self.id, "model%i" % self.id, "predictions%i.dat" % self.id])
        
        with open("predictions%i.dat" % self.id) as f:
            #~ pred = np.int32([p.split()[0] for p in f.readlines()])
            pred = np.float64([p.split() for p in f.readlines()])
        
        diff = (pred[:,1:].max(1))[:,np.newaxis]-pred[:,1:]
        mindiff = [d[d!=0].min() for d in diff]
        
        train_correct=0
        for i,v in enumerate(nodeOrder):
            self.G.node[v]["prediction"]=self.cindx_lutr[int(pred[i,0]-1)]
            self.G.node[v]["margin"]=mindiff[i]
            if self.G.node[v].has_key("label"):
                train_correct+=self.G.node[v]["prediction"]==self.G.node[v]["label"]
        #~ print train_correct
    
    
    
    ###################################################
    
        
        



class Display(object):
    def __init__(self,blockmodel,exp,Korder=None,nn=5):
        self.blockmodel=blockmodel
        self.plotRoles(nn,exp=exp,Korder=Korder)
    
    def saveModel(self,filename):
        
        for link in self.blockmodel.G.edges_iter(data=True):
            link[2]["group"]=np.argmax(link[2]["Zs"])
        
        data = json_graph.node_link_data(self.blockmodel.G)
        jsonFile=filename+".json"
        print jsonFile
        with open(jsonFile,"w") as f:
            f.write(json.dumps(data,default=str))
    
    def plotRoles(self,n_neighbors=5,n_components=1,exp="seafoodweb-feeding",Korder=None,nodelabels=True):
        
        yticklabels = {"seafoodweb-feeding" :["Primary","Carniv","Carn/Nec","Detriv","Herb/Det","Omniv"],
                                "karate" : ["President","Instructor"],
                                "word" : ["Adjective","Noun"]}[exp]
        blockmodel = self.blockmodel
        
        if Korder is None:
            Korder=range(blockmodel.K)
        K=len(Korder)
        
        Z_hat = np.empty((blockmodel.N,K))
        
        classlist="../data/%s_labels.txt" % exp
        with open(classlist) as f:
            labels = f.readlines()
        class_lut=dict((label.split()[0],label.split()[1].strip()) for label in labels)
        
        classes=np.unique(class_lut.values())
        
        Ytrue=[]
        nodes=[]
        v=0
        #~ words=[]
        for node,nodeData in blockmodel.G.nodes_iter(data=True):
            if not class_lut[node]=="0":
                Z_hat[v,:] = nodeData["nv"][Korder]/np.sum(nodeData["nv"])
                Ytrue.append(class_lut[node])
                nodes.append(node)
                if np.max(Z_hat[v,:])>0.5:
                    if (np.argmax(Z_hat[v,:])+1==2) or (np.argmax(Z_hat[v,:])+1==3):
                        print int(node)-1,np.argmax(Z_hat[v,:])+1, class_lut[node],np.max(Z_hat[v,:])
                v+=1
        nodes=np.array(nodes)
        #~ print classes,Ytrue
        #~ for c in classes:
            #~ nodes_of_class_c = [k for k,v in class_lut.iteritems() if v==c]
            #~ v=nodes_of_class_c[int(np.floor(np.random.rand()*len(nodes_of_class_c)))]
        
        cmap=plt.cm.bone_r
        #~ Ytrue = blockmodel.Ytrue
        #~ Z_hat = blockmodel.Z_hat
        plt.ion()        
        #~ fig=plt.figure()
        #~ plt.figure(figsize=(8,6))
        fig, ax = plt.subplots(1,figsize=(6,7))
        interactions=blockmodel.G.graph["c"][np.ix_(Korder,Korder)]/blockmodel.I
        #~ interactions/=(interactions+.01)
        interactions=np.log(interactions*blockmodel.I)
        interactions[interactions<0]=0
        #~ plt.pcolor(interactions,cmap=cmap,edgecolors="k")
        _pcolormesh.pcolormesh(fig,plt.gca(),interactions,edgecolors=(0.5,0,0))
        if exp=="karate":
            #~ plt.text(K+0.7,K/2,"log (number of links)", va="center", rotation=90, size=24)
            plt.text(K/2,-2,"log (number of links)", ha="center", size=24)
            plt.xticks(np.arange(K)+.5,np.arange(K)+1,fontsize=24)
            plt.yticks(np.arange(K)+.5,np.arange(K)+1,fontsize=24)
        elif exp=="word":
            plt.text(K/2,-2,"log (number of links)", ha="center", size=24)
            plt.xticks(np.arange(K)+.5,np.arange(K)+1,fontsize=24)
            plt.yticks(np.arange(K)+.5,np.arange(K)+1,fontsize=24)
        else:
            plt.text(K+2,K/2,"log (number of links)", va="center", rotation=90, size=24)
            plt.xticks(np.arange(K)+.5,np.arange(K)+1,fontsize=20)
            plt.yticks(np.arange(K)+.5,np.arange(K)+1,fontsize=20)
        
        
        plt.tight_layout(pad=0.1)
        plt.savefig("%sblocks.pdf" % exp,bbox_inches='tight')
        Ytrue=np.array(Ytrue)
        
        order = np.argsort(Ytrue)
        Z_hat = Z_hat[order]
        nodes = nodes[order]
        
        startclass=0
        plt.figure(figsize=(8,8))
        
        yticklocs = []
        plt.plot([0,K],[startclass,startclass],"b--",lw=6,alpha=.6)
        #~ for c in classes[1:]:
        for c in classes:
            print c
            endclass = startclass+np.sum(Ytrue==c)
            print startclass,endclass
            
            Y = manifold.Isomap(n_neighbors, n_components,eigen_solver='dense').fit_transform(Z_hat[startclass:endclass,:])
            
            order[startclass:endclass] = startclass + np.argsort(Y[:,0])
            plt.plot([0,K],[endclass,endclass],"b--",lw=6,alpha=.6)
            yticklocs.append(startclass+(endclass-startclass)/2.)
            startclass=endclass
        Z_hat = Z_hat[order]
        nodes= nodes[order]
        plt.pcolor(Z_hat,cmap=cmap)# determines empty positions
        if exp=="seafoodweb-feeding":
            plt.fill_betweenx([0,blockmodel.N],[1,1], [5,5], facecolor=(0,.5,0),alpha=0.3)
            plt.fill_betweenx([0,blockmodel.N],[5,5], [K,K], facecolor=(0,0,.5),alpha=0.3)
            plt.text(3,blockmodel.N+5,"eats-plants",ha="center",color=(0,.5,0),size=22)
            plt.text(8.5,blockmodel.N+5,"eats-animals",ha="center",color=(0,0,.5),size=22)
            plt.yticks(yticklocs,yticklabels,fontsize=22)#,rotation=45)
        elif exp=="karate":
            plt.yticks(np.arange(blockmodel.N)+.5,nodes,fontsize=18)
            plt.fill_betweenx([0,blockmodel.N],[1,1], [2,2], facecolor=(44/255.,160/255.,44/255.),alpha=0.5)
            plt.fill_betweenx([0,blockmodel.N],[0,0], [1,1], facecolor=(1,127/255.,14/255.),alpha=0.5)
            plt.fill_betweenx([0,blockmodel.N],[2,2], [3,3], facecolor=(31/255.,119/255.,180/255.),alpha=0.5)
            plt.text(1.5,blockmodel.N+.2,"Inter-faction",ha="center",color=(44/255.,160/255.,44/255.),size=20)
            plt.text(0.5,blockmodel.N+.2,"Instructor",ha="center",color=(1,127/255.,14/255.),size=22)
            plt.text(2.5,blockmodel.N+.2,"President",ha="center",color=(31/255.,119/255.,180/255.),size=22)
        elif exp=="word":
            plt.yticks(yticklocs,yticklabels,fontsize=22,rotation=90)
        plt.xticks(np.arange(K)+.5,np.arange(K)+1,fontsize=24)
        plt.ylim(0,blockmodel.N)
        plt.colorbar()
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(plt.gca())
        cax = divider.append_axes("top", "4%", pad="0%")
        #~ print dir(cax)
        plt.tight_layout(pad=0.1)
        cax.set_visible(False)
        plt.savefig("%smems.pdf" % exp,bbox_inches='tight')


if __name__=="__main__":
    print "supervisedBM v0.6"
    import argparse
    parser = argparse.ArgumentParser()
    networkfiletype = parser.add_mutually_exclusive_group()
    networkfiletype.add_argument("-a", "--adjlist", action="store_true", help="adjacency list network format (default)")
    #~ networkfiletype.add_argument("-q", "--quiet", action="store_true") # we can add functionality to use other networkx input methods here.
    parser.add_argument("networkfile", help="Input network file")
    parser.add_argument("-c","--classlist", help="File of class labels list")
    parser.add_argument("-o","--outputFile", help="Output File")
    parser.add_argument("-p","--cen", help="Censor prob")
    args = parser.parse_args()
    
    
    if args.adjlist:
        inG = nx.read_adjlist(args.networkfile, create_using=nx.DiGraph())
    else:
        pass
    if args.classlist:
        with open(args.classlist) as f:
            labels = f.readlines()
        for label in labels:
            v,l = label.split()
            if args.cen:
                if np.random.rand()<float(args.cen):
                    inG.node[v]["label"]=l.strip()
            else:
                inG.node[v]["label"]=l.strip()
    

    sb = MixedMembership(inG)
    sb.infer()
    tt=0
    yt=0
    for label in labels:
        v,l = label.split()
        print l.strip(),inG.node[v]["prediction"],inG.node[v]["cindx"] is None
        if inG.node[v]["cindx"] is None:
            tt+=1.
            yt+=l.strip()==inG.node[v]["prediction"]
    print tt,yt,yt/tt
    if args.outputFile:
        with open(args.outputFile,"w") as f:
            for node,v in sb.G.nodes_iter(data=True):
                f.write("%s\t%s\n" % (str(node),str(v["prediction"])))