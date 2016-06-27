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
from findParent import *
import random

'''@function: typical function to run by command
   @input   : 
   @output  : arguments
'''
def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--Operon","-i", help="Operon file name, which contains gene block for each genome")
    parser.add_argument("--TreeFile","-t", help="Tree file name")
    args = parser.parse_args()
    return args
    
'''@function: parse the file that has infor of gene blocks for each leaf
   @input   : file
   @output  : dictionary, key is genome name, value is the gene block
'''
def parsing(file):
    genomes={}
    myfile = open(file,'r')
    for line in myfile.readlines():
        if line[0]=='N':
            item = line.split(':')
            name = item[0]
            gene_blocks= item[1].split('\n')[0]
            genomes[name]=gene_blocks
    return genomes
###############################################################################
# Helper function
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
    
    
'''@function: set initial value data for each node as empty Set, and set name for
              initial node
   @input   : tree in nwk format
   @output  : tree in nwk format
'''
def set_initial_value(rooted_tree):
    count = 0
    for node in rooted_tree.traverse("postorder"):
        node.add_features(data=set())
        if node.name == '':
            count +=1
            node.name = 'Node'+ ' ' + str(count)
    return rooted_tree
   
'''@function: set gene block data for each leaf node
   @input   : tree in nwk format,genomes dictionaru
   @output  : tree in nwk format
'''   
def set_leaf_gene_block(rooted_tree,genomes):
    leaves= rooted_tree.get_leaves()
    for leaf in leaves:
        mylist= (leaf.name).split('_')
        name = mylist[2]+'_'+mylist[3]
        leaf.add_features(gene_block=genomes[name])
    return rooted_tree
    
'''@function: set initial value data for leaf based on whether it has the gene 
              in the gene block
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
            info = TextFace(node.initial)
            node.add_face(info,column=0,position = 'branch-top')
    return rooted_tree
###############################################################################
# Functions to minimze deletion, duplication and  split cost.
###############################################################################
'''@function: Globablly minimize deletion cost by provide a set of genes
              to be included in inner node., and output the total deletion events.
   @input   : tree in nwk format,and set of genes
   @output  : tree in nwk format, and total_count of deletion event
'''
def minimize_del(rooted_tree,genes):
    total_count = 0     
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
            if not node.is_leaf():# and not node.is_root():
                for child in node.get_children():
                    count += abs(node.data-child.data)
        print 'gene :',gene 
        print 'count :',count
            #    print 'inner: ', node.name, node.data, node.total_sum, node.total_leaf
            #else:
            #    print 'leaf: \t', node.name, node.data
        total_count += count
    return (rooted_tree,total_count)
'''@function: Locally minimize split costusing the set of genes
              to be included in inner node., and output the total split events
   @input   : tree in nwk format,and set of genes
   @output  : tree in nwk format, and total_count of split event
'''
def minimize_split(rooted_tree):
    total_count = 0
    for node in rooted_tree.traverse('postorder'):
        if not node.is_leaf():
            node.add_features(initial=[])
            children = node.get_children()
            children_blocks=[]
            for child in children:
                gene_set = node.genes
                if child.is_leaf(): # check if the child is a leaf
                    processed = (child.gene_block).split('|')
                    children_blocks.append(reduce_gene(processed,gene_set))
                else:
                    children_blocks.append(reduce_gene(child.initial,gene_set))
            # check if the 2 related set of gene blocks has same number of block:
            # also provide a count for split cost
            count = 0
            number_of_block1 = len(children_blocks[0])
            number_of_block2 = len(children_blocks[1])
            if number_of_block1 == number_of_block2:
                node.initial = random.sample(children_blocks,1)[0]

            else:
                if number_of_block1 !=0 and number_of_block2 !=0:
                    count = abs(number_of_block1-number_of_block2)
                dic={number_of_block1:children_blocks[0],
                     number_of_block2:children_blocks[1]}  
                if number_of_block1 > number_of_block2:
                    node.initial = dic[number_of_block1]
                else:
                    node.initial = dic[number_of_block2]
            print node.name, node.initial, children_blocks[0], children_blocks[1]
            total_count += count # the minimum count will always be the differences
            
            '''
            ancestor = node.get_ancestors()
            if len(ancestor) == 0 : # this mean it is the root
                node.initial = random.sample(node.children_blocks,1)
            else:
                extra_children_info = ancestor[0].get_children()  '''
    return (rooted_tree,total_count)
###############################################################################
# Main function to reconstruct
###############################################################################
if __name__ == "__main__":
    args = get_arguments()
    rooted_tree= Tree(args.TreeFile)
    genomes = parsing(args.Operon)
    rooted_tree = set_leaf_gene_block(rooted_tree,genomes) # assign the gene block infofor the leaf from the genomes dic
    rooted_tree = set_inner_genes(rooted_tree) # set inner node's gene set
    total_count = 0 
    reference = rooted_tree.search_nodes(name='Escherichia_coli_NC_000913')
    reference_block = reference[0].gene_block
    genes = setOfGene(reference_block) 
    rooted_tree,total_count_del   =  minimize_del(rooted_tree,genes)
    rooted_tree,total_count_split =  minimize_split(rooted_tree)
    print 'total_count_del: ',total_count_del
    print 'minimize_split: ',total_count_split
    rooted_tree = display(rooted_tree)
    tree_style = TreeStyle()
    tree_style.show_leaf_name = False
    tree_style.title.add_face(TextFace('total_count_del: '+str(total_count_del)+ '  '
                                       'total_count_split:'+str(total_count_split)),column=0)
    rooted_tree.show(tree_style=tree_style)
