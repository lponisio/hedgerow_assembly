import numpy as np
import networkx as nx


def writeToFile(file,format):
    tree=createNxGraphFromConsensus(file)
    if format=='edgelist':
        nx.write_edgelist(tree, file.split('.')[0]+'.edgelist')
    elif format=='gml':
        nx.write_gml(tree, file.split('.')[0]+'.gml')
    elif format=='graphml':
        nx.write_graphml(tree, file.split('.')[0]+'.graphml')
    elif format=='pajek':
        nx.write_pajek(tree, file.split('.')[0]+'.net')
    
    
def createNxGraphFromConsensus(file):
        ## Creates a networkx graph from the consensus file
        with open(file) as f:
            rows=f.readlines()
        #~ rows.reverse()
        
        firstrow=-1
        hrg=False
        print file
        if file.split(".")[1]=="hrg":
            hrg=True
            firstrow=0
            rows=rows[1:]
        
        tree = nx.DiGraph()
        tree.add_node("root")
        
        if len(rows) > 0:
            mg=_addnodes(tree,"root",rows,firstrow,hrg,0)
        
        maxgroups=len(rows)
        minNodes=np.zeros(maxgroups)+np.inf
        
        for node,nodeData in tree.nodes_iter(data=True):
            if not (node.startswith("D") or node.startswith("r")):
                #~ print node, nodeData
                if minNodes[nodeData["group"]]>int(node):
                    minNodes[nodeData["group"]]=int(node)
        
        newGroups = dict((n,i) for i,n in enumerate(np.argsort(minNodes)))
        #~ print minNodes,newGroups
        for node,nodeData in tree.nodes_iter(data=True):
            if not (node.startswith("D") or node.startswith("r")):
                #~ print node, nodeData["group"],newGroups[nodeData["group"]]
                nodeData["group"]=newGroups[nodeData["group"]]
        
        return tree


def _addnodes(tree,parent,consensusTree,index,hrg=False,group=0):
    ##Recursively add nodes to tree. Function used by createDendroTreeFromConsensus
    nextgroup=group+1
    row=consensusTree[index]
    if hrg:
        leafGenerator= (row.split()[i] for i in [4,5,7,8])
    else:
        leafGenerator= (val for val in row.split("\t")[3].split()[2::])
    for leaf in leafGenerator:
        nodeType = leafGenerator.next().strip("()") 
        if nodeType == "D":
            #~ node = parent.new_child(taxon=tree.taxon_set.require_taxon(label=leaf+nodeType))
            node = nodeType+leaf
            tree.add_node(node)
            tree.add_edge(parent,node)
            _addnodes(tree,node,consensusTree,int(leaf),hrg,nextgroup)
            nextgroup+=1
        else:
            #~ print leaf
            tree.add_node(leaf)
            tree.add_edge(parent,leaf)
            if index==-1:
                tree.node[leaf]["group"]=len(consensusTree)-1
            else:
                tree.node[leaf]["group"]=index
            #~ node = parent.new_child(taxon=tree.taxon_set.require_taxon(label=leaf))
    return nextgroup