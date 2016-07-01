#!/usr/bin/env python
''' Author : Huy Nguyen
    Project: provide global approach to minimize deletion event
             and finding local split event based prioritizing deletion event.
    Start  : 19/06/2016
    End    : 15/07/2016
'''
from __future__ import division
from ete3 import *
import argparse
import os
import random
from file_handle import *

###############################################################################
# Helper function to set data structure
###############################################################################

'''@function: set the genes that will appear in each node as dictionary, and
              the deletion, split, duplication events.
   @input   : tree in nwk format
   @output  : tree in nwk format
'''
def set_inner_genes(rooted_tree,genes):
    for node in rooted_tree.traverse("levelorder"):
        node.add_features(deletion=[0,0],duplication=[0,0],split=[0,0])
        node.add_features(data={})
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
    
    
   
'''@function: set gene block data for each leaf node, and name for inner node
   @input   : tree in nwk format,genomes dictionaru
   @output  : tree in nwk format
'''   
def set_leaf_gene_block(rooted_tree,genomes):
    count = 0
    for node in rooted_tree.traverse("postorder"):
        if node.name == '':
            count +=1
            node.name = 'Node'+ ' ' + str(count)
        else:
            mylist= (node.name).split('_')
            name = mylist[2]+'_'+mylist[3]
            node.add_features(gene_block=genomes[name])
    return rooted_tree
    
    
'''@function: Given gene block, remove the gene that not in the set of genes 
              given from minimize deletion method
   @input   : gene blocks
   @output  : list of gene blocks
'''
def reduce_gene(gene_block,genes):
    result=[]
    for block in gene_block: # iterate trhough each gene block in the genome
        new_block =''
        for gene in block: # iterate through each gene in the gene block
            if gene in genes:
                new_block +=gene # add gene to the string
        if len(new_block)>0:
            result.append(''.join(sorted(new_block))) # add the sorted to the result
    return result

### display the info
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
            info = TextFace(node.data)
            node.add_face(info,column=0,position = 'branch-right')
        else: # display gene set if inner node
            info = TextFace(node.split)
            node.add_face(info,column=0,position = 'branch-top')
            info = TextFace(node.initial)
            node.add_face(info,column=0,position = 'branch-top')
    return rooted_tree
    
#######################################################################################
# Helper functions to calculate the edit distance
#######################################################################################
### Deletion event helper
'''@function: Calculate deletion distance for all gene at each inner node
   @input   : tree in nwk formatgene_set
   @output  : tree in nwk format
'''
def del_distance(rooted_tree,genes):
    for node in rooted_tree.traverse('postorder'):
        if not node.is_leaf():# and not node.is_root():
            for child in node.get_children():
                value = 0
                for gene in genes:
                    value = abs(node.data-child.data)
                    node.deletion[0]+=value
    return rooted_tree
    
'''@function: Calculate accumulation deletion distance for each gene at each inner node
   @input   : tree in nwk format
   @output  : tree in nwk format
'''    
def accumulate_del(rooted_tree):
    for node in rooted_tree.traverse('postorder'):
        if not node.is_leaf():# and not node.is_root():
            node.deletion[1]+=node.deletion[0]
            for child in node.get_children():                
                if not child.is_leaf():
                    node.deletion[1]+=child.deletion[1]
    return rooted_tree
    
### Duplication event helper
    
'''@function: Calculate accumulation deletion distance for each gene at each inner node
   @input   : tree in nwk format
   @output  : tree in nwk format
'''  
def find_Dup():
    return 
    
### Split event helper
'''@function: From the set of genes of inner node, get the related gene block 
              in leaf nodes. Assign this the processed gene_block to attribute 
              proccessed block of the leaf node
   @input   : tree in nwk format, leave nodes
   @output  : tree in nwk format, matrix score of number of block
'''   
def initialize_block_number(rooted_tree,leaves):
    for leaf in leaves:
        ancestors = leaf.get_ancestors()
        gene_set  = ancestors[0].genes # get the gene set of the ancestor
        processed = leaf.gene_block.split('|') 
        processed = reduce_gene(processed,gene_set)
        leaf.add_features(processed_block=processed) # assign the processed block
        # print leaf.name, leaf.gene_block, gene_set
        leaf.data.add(len(processed)) # assign the number of block of the processed
    return rooted_tree
    
'''@function: Fitch algorithm
   @input   : tree in nwk format, set of genes
   @output  : tree in nwk format
'''       
def Fitch(rooted_tree,genes):
    # using of Fitch algorithm 
    
    # traverse Tree in post-order
    for node in rooted_tree.traverse('postorder'):
        if not node.is_leaf():
            children = node.get_children()
            for gene in genes:
                intersect = (children[0].data[gene]).intersection(children[1].data[gene])
                if len(intersect) == 0:
                    node.data[gene] = (children[0].data[gene]).union(children[1].data[gene])
                else:
                    node.data[gene] = intersect
                    
    # traverse top-down     
    for node in rooted_tree.traverse('levelorder'):
        for gene in genes:
            if node.is_root(): # for the root 
                # if the root has 2 candidatnode.add_features()e, randomly choose 1, and get the numeric value
                node.data[gene] = (random.sample(node.data[gene],1))[0] 
            else:
                # for children node, first check the data from the ancestor
                ancestors = node.get_ancestors() # get the list of ancestor
                data = ancestors[0].data[gene] # get the data from the parent
                if data in node.data[gene]:# check if the node.data has value equal to its parent data
                    node.data[gene] = data
                else:
                    node.data[gene] = (random.sample(node.data[gene],1))[0]
    return rooted_tree
'''@function: Calculate accumulation split distance for each gene at each inner node
   @input   : tree in nwk format
   @output  : tree in nwk format
'''    
def accumulate_split(rooted_tree):
    for node in rooted_tree.traverse('postorder'):
        if not node.is_leaf():# and not node.is_root():
            node.split[1]+=node.split[0]
            for child in node.get_children():                
                if not child.is_leaf():
                    node.split[1]+=child.split[1]
    return rooted_tree
                    
###############################################################################
# Functions to minimze deletion, duplication and  split cost.
###############################################################################
'''@function: Globablly minimize deletion cost by provide a set of genes
              to be included in inner node., and output the info of deletion
              distances at each inner node
   @input   : tree in nwk format,and set of genes and leaves node
   @output  : tree in nwk format, and total_count of deletion event
'''
def minimize_del(rooted_tree,genes):
    rooted_tree = Fitch(rooted_tree,genes)
    # calculate accumulation deletion distances
    rooted_tree = accumulate_del(rooted_tree)
    return rooted_tree
    
    
'''@function: pseudo globally minimize split costusing the set of genes
              to be included in inner node., and output the total split events
   @input   : tree in nwk format,and set of genes
   @output  : tree in nwk format, and total_count of split event
'''
def minimize_split(rooted_tree):
    rooted_tree = Fitch(rooted_tree)
    for node in rooted_tree.traverse('postorder'):
        if not node.is_leaf():
            node.add_features(initial=[])       
            children = node.get_children()
            children_blocks=[]
            for child in children:
                gene_set = node.genes
                if child.is_leaf(): # check if the child is a leaf
                    children_blocks.append(child.processed_block)
                else:
                    children_blocks.append(reduce_gene(child.initial,gene_set))
            # check if the 2 related set of gene blocks has same number of block:
            # also provide a count for split cost

            # set the split distances, check with the parent data
            number_of_block_parent = node.data
            number_of_block1 = len(children_blocks[0])
            number_of_block2 = len(children_blocks[1])
            if number_of_block1 == number_of_block2:
                node.initial = random.sample(children_blocks,1)[0]

            else:
                if number_of_block1 !=0 and number_of_block2 !=0:
                    node.split[0] = abs(number_of_block1-number_of_block2)
                abs_1 = abs(number_of_block_parent-number_of_block1)
                abs_2 = abs(number_of_block_parent-number_of_block2)
                dic={number_of_block1:children_blocks[0],
                     number_of_block2:children_blocks[1]}  
                if abs_2 > abs_1:
                    node.initial = dic[number_of_block1]    
                else:
                    node.initial = dic[number_of_block2]
    rooted_tree = accumulate_split(rooted_tree)    
    return rooted_tree
###############################################################################
# Main function to reconstruct
###############################################################################
if __name__ == "__main__":
    args = get_arguments()
    rooted_tree= Tree(args.TreeFile)
    mapping,genomes = parsing(args.InputDataDirectory)
    genes = set()
        # get the genes of the operon
    for key in mapping:
        genes.add(key)
    rooted_tree = set_leaf_gene_block(rooted_tree,genomes) # assign the gene block infofor the leaf from the genomes dic
    rooted_tree = set_inner_genes(rooted_tree,genes) # set inner node's gene set dictionary, distances, and block number
    leaves      = rooted_tree.get_leaves()
    rooted_tree =  minimize_del(rooted_tree,genes)
    rooted_tree = initialize_block_number(rooted_tree,leaves)
    rooted_tree =  minimize_split(rooted_tree)
    rooted_tree = display(rooted_tree)
    tree_style = TreeStyle()
    tree_style.show_leaf_name = False
    rooted_tree.show(tree_style=tree_style)
