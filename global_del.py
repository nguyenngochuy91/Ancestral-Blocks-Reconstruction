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
import random
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
'''@function: set the genes that will appear in each inner node
   @input   : tree in nwk format
   @output  : tree in nwk format
'''
def set_inner_genes(rooted_tree):
    for node in rooted_tree.traverse("levelorder"):
        if not node.is_leaf():
            node.add_features(genes=set())
    return rooted_tree
    
    
'''@function: set initial value data for each node as empty Set
   @input   : tree in nwk format
   @output  : tree in nwk format
'''
def set_initial_value(rooted_tree):
    for node in rooted_tree.traverse("levelorder"):
        node.add_features(data=set())
    return rooted_tree
    
'''@function: set initial value data for leaf based on whether it has the gene g
   @input   : tree in nwk format,and gene g
   @output  : tree in nwk format
'''

def set_leaf_data(rooted_tree,gene):
    leaves= rooted_tree.get_leaves()
    for leaf in leaves:
        if gene in leaf.gene_block: # if gene appear in the genome gene block
            (leaf.data).add(1)
        else:
            (leaf.data).add(0)
    return rooted_tree
    
'''@function: Transform some of the inner node to be counted as leaf, and remove
              some leaf node out of calculation
   @input   : tree in nwk format,and gene g
   @output  : tree in nwk format
'''

    
'''@function: Display the tree, each inner node has info about genes to include.
              Leaf node has gene block
   @input   : tree in nwk format,and gene g
   @output  : tree in nwk format
'''
def display(rooted_tree):
    for node in rooted_tree.traverse('postorder'):
        if node.is_leaf(): # display gene block if leaf node
            info = TextFace(node.gene_block)
            node.add_face(info,column=0,position = 'branch-right')
        else: # display gene set if inner node
            info = TextFace(node.genes)
            node.add_face(info,column=0,position = 'branch-top')
    return rooted_tree
###############################################################################
## Main function
###############################################################################
if __name__ == "__main__":
    args = get_arguments()
    rooted_tree= Tree(args.Operon)
    rooted_tree = set_inner_genes(rooted_tree) # set inner node's gene set
    total_count = 0 
    reference = rooted_tree.search_nodes(name='Escherichia_coli_NC_000913')
    reference_block = reference[0].gene_block
    genes = setOfGene(reference_block) 

    for gene in genes: # iterate through all genes of reference
    
        rooted_tree = set_initial_value(rooted_tree) #set initial data value as empty Set
        # set data for the leavesrooted_tree.data = 1
        rooted_tree = set_leaf_data(rooted_tree,gene)
        
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
        print 'gene :',gene 
        print 'count :',count
    print(total_count)
    rooted_tree = display(rooted_tree)
    tree_style = TreeStyle()
    tree_style.show_leaf_name = False
    tree_style.title.add_face(TextFace('total_count: '+str(total_count)),column=0)
    rooted_tree.show(tree_style=tree_style)
