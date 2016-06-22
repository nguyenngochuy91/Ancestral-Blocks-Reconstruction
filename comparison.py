#!/usr/bin/env python
from __future__ import division
from ete3 import *
import argparse
import os
from findParent import *



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
              # Set data as None on node
              node.add_features(data=None)
            total_count = 0 
            reference = rooted_tree.search_nodes(name='Escherichia_coli_NC_000913')
            leaves= rooted_tree.get_leaves()
            reference_block = reference[0].gene_block
            genes = setOfGene(reference_block)
            for gene in genes: # iterate through all genes of reference
                for leaf in leaves:
                    if gene in leaf.gene_block: # if gene appear in the genome gene block
                        leaf.data = 1
                    else:
                        leaf.data = 0
                                                        
                for node in rooted_tree.traverse('postorder'):
                    if not node.is_leaf():                
                        children = node.get_children()
                        small_sum =0
                        number_leaf =0 
                        for child in children:
                            if child.is_leaf():
                                number_leaf +=1
                                small_sum += child.data
                            else:
                                number_leaf += child.total_leaf
                                small_sum += child.total_sum
                        node.add_features(total_sum= small_sum, total_leaf=number_leaf)
                
                if (rooted_tree.total_sum/rooted_tree.total_leaf) >=.5:
                    rooted_tree.data = 1      
                for node in rooted_tree.traverse('levelorder'):
                    if not node.is_leaf():
                        if node.total_sum/node.total_leaf < .5:
                           node.data = 0  
                        elif node.total_sum/node.total_leaf >.5:
                            node.data = 1
                        elif node.total_sum/node.total_leaf ==.5:
                            ancestors = node.get_ancestors()
                            if len(ancestors) !=0:
                                node.data = ancestors[0].data
                    node.add_face(TextFace(node.data), column =0, position ="branch-top")        
                    
                count =0
                for node in rooted_tree.traverse('postorder'):
                    if not node.is_leaf() and not node.is_root():
                        for child in node.get_children():
                            count += abs(node.data-child.data)
                    #    print 'inner: ', node.name, node.data, node.total_sum, node.total_leaf
                    #else:
                    #    print 'leaf: \t', node.name, node.data
                total_count += count
                tree_style = TreeStyle()
                tree_style.show_leaf_name = False
                tree_style.title.add_face(TextFace('count: '+str(count)),column=0)
                # rooted_tree.show(tree_style=tree_style)
            mystring = f+':'+str(deletion_total)+' '+str(total_count)+' '+str(percent)+'\n'
            outfile.write(mystring)
    outfile.close()
    