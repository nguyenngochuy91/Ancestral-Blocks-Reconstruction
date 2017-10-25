#!/usr/bin/env python
''' Author  : Huy Nguyen
    Program : Boostrapping
    Start   : 09/10/2017
    End     : 09/18/2016
'''
from itertools import chain, combinations
from ete3 import Tree,TextFace,TreeStyle
import argparse
from findParent_local import setOfBlocks,setOfGene
from findParent_global import set_inner_genes,minimize_del,initialize_block_number,minimize_split,find_dup,minimize_dup
import os

# traverse and get the file
def traverseAll(path):
    res=[]
    for root,dirs,files in os.walk(path):
        for f in files:
            res.append(root+f)
    return res
    
def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputDir","-i", help="input directory name(reconstruction_global)")
    parser.add_argument("--outDir","-o", help="Output Directory for the bootstrap")
#    parser.add_argument("--ref","-r", help="reference genome (ncbi accession number)")
    parser.add_argument("--group","-g", help="Group by text file(result/group.txt)")
    args = parser.parse_args()
    return args


'''@function: parsing the genes name into a set
   @input   : textfile
   @output  : set
'''
def parsingMap(infile):
    infile = open(infile,'r')
    line = infile.readlines()
    line = line[0].split()
    genes  = set()
    for info in line:
        genes.add(info.split(',')[1])
    return genes

### generate the powerset of a given set
def powerset(iterable):
    """
    powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)
    """
    xs = list(iterable)
    # note we return an iterator rather than a list
    return chain.from_iterable(combinations(xs,n) for n in range(1,len(xs)))
    
    ### given the reference gene block, the current initial gene block of the inner node, generate the sample to test
    

### find the relevant gene block
'''@function: given initial block and its 2 children block,calculate the relevant gene block
   @input   : block, intersection
   @output  : string
'''     
def relevant(block,intersection):
    string = ""
    for letter in block:
        if letter == "|":
            string+=letter
        else:
            if letter in intersection:
                string+=letter
    return string

### calculate split
def getSplit(string):
    count = 0
    block = setOfBlocks(string)
    for item in block:
        if len(block)!= 0:
            count+=1
    return count
### check if the new gene block is valid,
def isValid(block):
    blocks = setOfBlocks(block)
    check  = False
    for item in blocks:
        if len(item)>=2:
            return True
    return check
### if the block is valid, reformat the block
def reformat(block):
    blocks = setOfBlocks(block)
    res = []

    for item in blocks:
        if len(item)>0:
            res.append(''.join(sorted(item)))

    return '|'.join(res)
### get the duplicated gene set of a given block
def getDuplication(string):
    dup = set()
    block = setOfBlocks(string)
    for item in block:
        dic = {}
        for letter in item:
            if letter not in dic:
                dic[letter]=1
            else:
                dup.add(letter)
    return dup
### class block to store geneblock info
class Block(object):
    def __init__(self,geneBlock,deletion,duplication ,split):
        self.geneBlock = geneBlock
        self.deletion = deletion
        self.duplication = duplication
        self.split  = split
    
    def calculateDistance(self,geneBlock):
        genes1 = setOfGene(self.geneBlock)
        genes2 = setOfGene(geneBlock)
        del_distance = len(genes1.symmetric_difference(genes2))
        intersection = genes1.intersection(genes2)
        deletion = [del_distance,self.getDeletion()]
            
        # remove gene that is not in itnersection
        string1 = reformat(relevant(geneBlock,intersection))
        string2 = reformat(relevant(self.geneBlock,intersection))
        
        # duplciation
        dup1 = getDuplication(string1)
        dup2 = getDuplication(string2)
        duplication = [abs(len(dup1)-len(dup2)),self.getDupication()]
        
        # split
        count1 = getSplit(string1)
        count2 = getSplit(string2)
        split = [abs(count1-count2),self.getSplit()]     
        
        distance = [deletion,duplication,split]      
        return distance
    
    def getDeletion(self):
        return int(self.deletion.split('|')[1])

    def getDupication(self):
        return int(self.duplication.split('|')[1])

    def getSplit(self):
        return int(self.split.split('|')[1])        
'''@function: given ref block, and current innitial, generate possible suboptimal sets
   @input   : string1, string2
   @output  : set of strings
'''   
def generateSample(node):
    res = {}
    geneInit = setOfGene(node.initial)
    children = [child for child in node.get_children()]
    childrenBlock  = []
    for child in children:
        if child.is_leaf():
            newBlock = Block(child.gene_block,child.deletion,child.duplication,child.split)
            childrenBlock.append(newBlock)
        else:
            newBlock = Block(child.initial,child.deletion,child.duplication,child.split)
            childrenBlock.append(newBlock)
    geneChild = [setOfGene(child.geneBlock) for child in childrenBlock]

    unionGenes = geneChild[0].union(geneChild[1])
    intersectionGenes = geneChild[0].intersection(geneChild[1])
    if geneInit == unionGenes:
        
        for gene in geneInit:
        
            temp = (node.initial).replace(gene,"")
            if not isValid(temp):
                continue
            if temp in res or temp == node.initial:
                continue
            temp = reformat(temp)
            distance1 = childrenBlock[0].calculateDistance(temp)
            distance2 = childrenBlock[1].calculateDistance(temp)

            distance = []
            for i in range(3):
                dist = []
                dist.append(distance1[i][0]+distance2[i][0])
                dist.append(dist[0]+distance1[i][1]+distance2[i][1])
                distance.append(dist)
            res[temp] = distance
    else:
        onlyOne = unionGenes-intersectionGenes
        powerSet = powerset(onlyOne)
        for subset in powerSet:
            # add to our initial
            toAdd = []
            toRemove = []
            for gene in subset:
                if gene not in node.initial:
                    toAdd.append(gene)
                else:
                    toRemove.append(gene)
            # add new gene to initial, this needs more logic
            if len(toAdd) >0:
                temp = node.initial
                for gene in toAdd:
                    temp +=gene
                if temp != node.initial:
                    
                    distance1 = childrenBlock[0].calculateDistance(temp)
                    distance2 = childrenBlock[1].calculateDistance(temp)
                    distance = []
                    for i in range(3):
                        dist = []
                        dist.append(distance1[i][0]+distance2[i][0])
                        dist.append(dist[0]+distance1[i][1]+distance2[i][1])
                        distance.append(dist)
                    res[temp] = distance
            # remove the gene from initial
            if len(toRemove) >0:
                temp = node.initial
                for gene in toRemove:
                    temp.replace(gene,"")
                if not isValid(temp) or temp == node.initial:
                    continue
                distance1 = childrenBlock[0].calculateDistance(temp)
                distance2 = childrenBlock[1].calculateDistance(temp)
                distance = []
                for i in range(3):
                    dist = []
                    dist.append(distance1[i][0]+distance2[i][0])
                    dist.append(dist[0]+distance1[i][1]+distance2[i][1])
                    distance.append(dist)
                res[temp] = distance
            
    return res

### newick file, get the reference gene blockres.add((temp,0))
'''@function: given a tree file, find the reference initial gene block 
   @input   : tree
   @output  : string (ref gene block)
'''   
def getRef(tree,ref):
    for node in tree.get_leaves():
        if ref in node.name:
            return node.gene_block
### newick file, storing info for each inner node about its cost so far, as well as computed a sample of sub optimal initial
'''@function: given a tree file, for each inner node, get the set of sub optimal value 
   @input   : tree
   @output  : tree
'''
def parseTree(tree):
    for node in tree.iter_descendants("postorder"):
        if not node.is_leaf():
            # create face contain initial set info
            node.add_features(sample= generateSample(node))
#            print (node.sample)
#    print (tree)
    return tree

'''@function: Reconstruct the newick tree file with gene block info for inner 
              node using local GLOBAL scheme
   @input   : tree in nwk format,and a dictionary between specie name and gene block for leaf, and set of genes
   @output  : tree in nwk format,gene g and a string of the info
'''
def reconstruct_global(tree,genes):
    leaves           =  tree.get_leaves() # get leave data so dont have to keep on calling 
    tree             =  minimize_del(tree,genes) # globally minimize deletion events, provide gene set for each inner node
    tree             =  initialize_block_number(tree,leaves) # using the gene set to get relevant gene block for each leaf
    tree             =  minimize_split(tree)
    check,tree,genes = find_dup(tree,leaves) 
    if check:
        tree  = minimize_dup(tree,genes)
    return tree
    
'''@function: set the genes that will appear in each node as dictionary, and
              the deletion, split, duplication events, and the genes set for each 
              inner node
   @input   : tree in nwk format
   @output  : tree in nwk format
'''

def set_inner_genes_special(rooted_tree,genes,name,distance):
    for node in rooted_tree.traverse("levelorder"):
        node.add_features(data={})
        node.add_features(genes=set())
        if node.name == name:
#            print (node.name)
            node.deletion = distance[0]
            node.duplication = distance[1]
            node.split = distance[2]
#            print (node.deletion,node.duplication,node.split)
        else:
            node.deletion = [0,0]
            node.duplication = [0,0]
            node.split = [0,0]
        if node.is_leaf():    
            for gene in genes:
                if gene in node.gene_block:
                    node.data[gene] = {1}
                else:
                    node.data[gene] = {0}
        else:
            for gene in genes:
                node.data[gene] = {0}
    return rooted_tree
    
### quick function to sum all the cost
def getTotalDistanceString(tree): 
    children= []
    for child in tree.get_children():
        children.append(child)
    deletion_total = 0
    duplication_total = 0
    split_total = 0
    
    for child in children:
        deletion_total+= int(child.deletion.split('|')[1])
        duplication_total+= int(child.duplication.split('|')[1])
        split_total+= int(child.split.split('|')[1])
    return (deletion_total,duplication_total,split_total)
    
### quick function to sum all the cost
def getTotalDistanceList(tree): 
    children= []
    for child in tree.get_children():
        children.append(child)
    deletion_total = 0
    duplication_total = 0
    split_total = 0
    
    for child in children:
        deletion_total+= (child.deletion[1])
        duplication_total+= (child.duplication[1])
        split_total+= (child.split[1])
    return (deletion_total,duplication_total,split_total)
if __name__ == "__main__":

    args     = get_arguments()
    inputDir = args.inputDir
    res      = traverseAll(inputDir)
    outputsession= args.outDir
    group    = args.group
    # try to create the boostrap directory
    try:
        os.mkdir(outputsession)
    except:
        print ("outdir already created")
    for operon in sorted(res):
        if "mapping" in operon:
            continue
        else:
            operonName = operon.split('/')[-1]
            print ("Boostraping operon:",operonName)
            mapping    = operon+"_mapping"
            # try to create a dir
            operonDir = outputsession+"/"+operonName
            try:
                os.mkdir(operonDir)
            except:
                print ("This operon directory is already created")
            data = operonDir +"/data"
            try:
                os.mkdir(data)
            except:
                print ("This data operon directory is already created")
                
            visualization = operonDir + "/visualization"
            try:
                os.mkdir(visualization)
            except:
                print ("This visualization operon directory is already created")
            tree     = Tree(operon)
            total1   = getTotalDistanceString(tree)
            genes    = parsingMap(operon+"_mapping")
            
            # create an info file that generate the information from all the reconstruction
            
            # get the gene block in reference genomes, generate the sameple
            tree = parseTree(tree) 
            lower = 1
            count =1
            outfile = open(data+"/analysis",'w')
            outfile.write("Our reconstruction cost:"+str(total1)+"\n")
            better =[] 
            # from the sample for each inner node, prune the tree and run the reconstruction, then generate the normal tree.
            for node in tree.iter_descendants("postorder"):
                if not node.is_leaf():
                    # get the name
                    name = node.name
                    sample = node.sample
                    if len(sample) == 0:
                        continue
        #            print ("Node to sample:",node.name)
        #            print ("sample set:",sample)
                    # pick out a sample from the sample list
#                    print (name,sample)
                    for candidate in sample:
                        
                        count+=1
        #                print ("String to check:",candidate)
                        sampleTree = Tree(operon)

#                        sampleTree.show()                
                        nodeInSample = sampleTree&name
                        nodeInSample.add_features(gene_block= None)
                        nodeInSample.gene_block = candidate
                        # detach the children of this nodeInSample
                        children = nodeInSample.get_children()
                        child1   = children[0].detach()
                        child2   = children[1].detach()

#                        sampleTree.show()
                        # get the distance of this in the sampleTree
                        distance = sample[candidate]
#                        print ("new distance:",distance)
                        sampleTree = set_inner_genes_special(sampleTree,genes,name,distance)
#                        print (nodeInSample.name,nodeInSample.deletion,nodeInSample.duplication,nodeInSample.split)
                        sampleTree = reconstruct_global(sampleTree,genes)

                        nodeInSample.add_child(child1)
                        nodeInSample.add_child(child2)
                        for node in sampleTree.iter_descendants("postorder"): 
                            if node.name == name:
                                node.add_features(modified= 1)
                            else:
                                node.add_features(modified= 0)
        #                sampleTree.show()outfile.write("Sample recosntruction:"+str(total2))
                        dataOutfile  = data+"/"+operonName+"_"+str(count)
                        total2 = getTotalDistanceList(sampleTree)
                        if sum(total1) <=sum(total2):
                            lower+=1
                        else:
                            print (dataOutfile)
                            better.append(operonName+"_"+str(count))
#                        print (nodeInSample.name,nodeInSample.deletion,nodeInSample.duplication,nodeInSample.split)
                        sampleTree.write(format=2, outfile=dataOutfile,features=['name',
                'initial','gene_block','deletion','duplication','split','modified'])
                        visualOutfile = visualization+ "/"+operonName+"_"+str(count)
                        cmd11 = './show_boostrap.py -i {} -g {} -o {} -m {}'.format(dataOutfile,group,visualOutfile,mapping)
                        os.system(cmd11)
            outfile.write("% that our reconstruction is better:"+str(lower/float(count)*100)+"\n")
            outfile.write("Reconstruction files with candidate in sample that have lower cost: \n")
            for item in better:
                outfile.write(item +"\n")
            outfile.close()
            

            
        
        
        
        
