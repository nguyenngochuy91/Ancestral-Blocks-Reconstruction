#!/usr/bin/env python
''' Author : Huy Nguyen
    Project: provide global approach to minimize deletion event

    Start  : 19/06/2016
    End    : 15/07/2016
'''
from __future__ import division
from ete3 import *
import argparse
import os
from findParent import *
'''@function: typical function to run by command
   @input   : 
   @output  : arguments
'''
def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--Operon","-i", help="Operon file name")
    args = parser.parse_args()
    return args
    
###############################################################################
## Helper function
###############################################################################
'''@function: set initial value data for each node as None:
   @input   : tree in nwk format
   @output  : tree in nwk format
'''
def set_initial_value(rooted_tree):
    for node in rooted_tree.traverse("levelorder"):
        node.add_features(data=None)
    return rooted_tree
    
'''@function: set initial value data for leaf based on whether it has the gene g
   @input   : tree in nwk format,and gene g
   @output  : tree in nwk format
'''

def set_leaf_data(rooted_tree,gene):
    leaves= rooted_tree.get_leaves()
    for leaf in leaves:
        if gene in leaf.gene_block: # if gene appear in the genome gene block
            leaf.data = 1
        else:
            leaf.data = 0
    return rooted_tree
    
'''@function: Transform some of the inner node to be counted as leaf, and remove
              some leaf node out of calculation
   @input   : tree in nwk format,and gene g
   @output  : tree in nwk format
'''

def set_new_leaf(rooted_tree):
    # serve to reduce amount of time to check, as well as provide more accurate
    # global optimization
    for node in rooted_tree.traverse('postorder'): # set the parent data to 1 if both child data is 1, and 
                                           # 0 if both child data is 0
        if not node.is_leaf():                
            children = node.get_children()
            if children[0].data == children[1].data and children[0].data!=None: # check children data
                node.data = children[0].data
    return rooted_tree
###############################################################################
## Main function
###############################################################################
if __name__ == "__main__":
    args = get_arguments()
    rooted_tree= Tree(args.Operon)
    rooted_tree = set_initial_value(rooted_tree) #set initial data value as None
    total_count = 0 
    reference = rooted_tree.search_nodes(name='Escherichia_coli_NC_000913')
    reference_block = reference[0].gene_block
    genes = setOfGene(reference_block) 

    for gene in genes: # iterate through all genes of reference
    
        # set data for the leaves
        rooted_tree = set_leaf_data(rooted_tree,gene)
        # set new leafs
        rooted_tree = set_new_leaf(rooted_tree)
        # have to check whether the leaf to account for (As well as inner node)
        for node in rooted_tree.traverse('postorder'):
            if not node.is_leaf() and node.data == None:
                children = node.get_children()
                small_sum =0
                number_leaf =0 
                for child in children:
                    if child.data != None: #this means that either this child is an intial leaf, or becomes the leaf after 
                                           # previous transformation
                        number_leaf +=1
                        small_sum += child.data
                    else:
                        number_leaf += child.total_leaf
                        small_sum += child.total_sum
                node.add_features(total_sum= small_sum, total_leaf=number_leaf)
        if (rooted_tree.total_sum/rooted_tree.total_leaf) >=.5:
            rooted_tree.data = 1      
        for node in rooted_tree.traverse('levelorder'):
            if not node.is_leaf() and node.data == None:
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
        print 'gene :',gene 
        print 'count :',count
    # rooted_tree.show(tree_style=tree_style)
    print(total_count)