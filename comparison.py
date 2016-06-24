#!/usr/bin/env python
from __future__ import division
from ete3 import *
import argparse
import os
from findParent import *
import random


'''
The algorithm idea is to do a global scheme as follow
1. Travel bottom up and do the 
'''

def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--Reconstruction","-i", help="Reconstruction directory")
    args = parser.parse_args()
    return args

def traverseAll(path):
    res=[]
    for root,dirs,files in os.walk(path):
        for f in files:
            res.append(root+'/'+f)
    return res
    
    
if __name__ == "__main__":
    args = get_arguments()
    res = traverseAll(args.Reconstruction)
    outfile=open('New_Comparison','w')
    for r in res:
        # check for DS.Store
        root,f = os.path.split(r)
        if "DS_Store" in f or "mapping" in f:
            continue
        else:
            rooted_tree= Tree(r)
    
            '''comparing part'''
            children= []
            for child in rooted_tree.get_children():
                children.append(child)
            deletion_cost1 = (children[0].deletion).split('|')[1]
            deletion_cost2 = (children[1].deletion).split('|')[1]
            deletion_total = int(deletion_cost1)+int(deletion_cost2)
            for node in rooted_tree.traverse("levelorder"):
              if not node.is_leaf():
                  node.add_features(genes=set())
            total_count = 0 
            reference = rooted_tree.search_nodes(name='Escherichia_coli_NC_000913')
            leaves= rooted_tree.get_leaves()
            reference_block = reference[0].gene_block
            genes = setOfGene(reference_block)
            for gene in genes: # iterate through all genes of reference
                for node in rooted_tree.traverse("levelorder"):
                    # Set data as None on node
                    node.add_features(data=set())
                for leaf in leaves:
                    if gene in leaf.gene_block: # if gene appear in the genome gene block
                        (leaf.data).add(1)
                    else:
                        (leaf.data).add(0)
                                                        
                # using of Fitch algorithm 
                # traverse Tree in post-order
                for node in rooted_tree.traverse('postorder'):
                    if not node.is_leaf():
                        children = node.get_children()
                        intersect = (children[0].data).intersection(children[1].data)
                        if len(intersect) == 0:
                            node.data = (children[0].data).union(children[1].data)
                        else:
                            node.data = intersect
                # traverse top-down 
                
                for node in rooted_tree.traverse('levelorder'):
                    if node.is_root(): # for the root 
                        # if the root has 2 candidate, randomly choose 1, and get the numeric value
                        node.data = (random.sample(node.data,1))[0] 
                    else:
                        # for children node, first check the data from the ancestor
                        ancestors = node.get_ancestors() # get the list of ancestor
                        data = ancestors[0].data # get the data from the parent
                        if data in node.data:# check if the node.data has value equal to its parent data
                            node.data =data
                        else:
                            node.data = (random.sample(node.data,1))[0]
                    if node.data ==1 and not node.is_leaf(): # include the gene into set of genes if data =1
                        (node.genes).add(gene)
                        
                count =0
                for node in rooted_tree.traverse('postorder'):
                    if not node.is_leaf() and not node.is_root():
                        for child in node.get_children():
                            count += abs(node.data-child.data)
                    #    print 'inner: ', node.name, node.data, node.total_sum, node.total_leaf
                    #else:
                    #    print 'leaf: \t', node.name, node.data
                total_count += count
            '''
            if 'paa' in f:
                print total_count
                tree_style = TreeStyle()
                tree_style.show_leaf_name = False
                tree_style.title.add_face(TextFace('total_count: '+str(total_count)),column=0)
                for node in rooted_tree.traverse('postorder'):
                    if node.is_leaf(): # display gene block if leaf node
                        info = TextFace(node.gene_block)
                        node.add_face(info,column=0,position = 'branch-right')
                    else: # display gene set if inner node
                        info = TextFace(node.genes)
                        node.add_face(info,column=0,position = 'branch-top')
                rooted_tree.show(tree_style=tree_style) '''
            mystring = f+':'+str(deletion_total)+' '+str(total_count)+'\n'
            outfile.write(mystring)
    outfile.close()
    
