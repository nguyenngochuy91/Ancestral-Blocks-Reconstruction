#!/usr/bin/env python
"""
Created on Wed Jul  6 16:14:23 2016

@author: huyn
@purpose: testing professor John output
"""
from ete3 import *
#######################################################################################
# Helper functions to calculate the edit distance
#######################################################################################
''' @function   : return set of gene in a Genome
    @input      : list (this list comprises of genes, and '|' as split event)
    @output     : set of gene
'''
def setOfGene(Genome):
    Genome_gene= set()
    for gene in Genome:
        if gene != '|':
            Genome_gene.add(gene)
    return Genome_gene
    
def del_distance(initial1,initial2,initial):
    count =0
    
    initial_gene = set()
    for block in initial:
        gene_set= setOfGene(block)
        initial_gene= initial_gene.union(gene_set)
        
    initial1_gene = set()
    for block in initial1:
        gene_set= setOfGene(block)
        initial1_gene= initial1_gene.union(gene_set)
        
    initial2_gene = set() 
    for block in initial2:
        gene_set= setOfGene(block)
        initial2_gene= initial2_gene.union(gene_set)
        
    count += (len(initial_gene-initial1_gene)+len(initial1_gene-initial_gene)+
    len(initial_gene-initial2_gene)+len(initial2_gene-initial_gene))
    return count,initial_gene

def reduction_gene(initial, gene_set):
    result= set()
    for element in initial:
        # create a copy for each element starting as empty
        copy =''
        for gene in element:
            # adding the gene in the element if the count is 2
            if gene in gene_set:
                copy += gene
        if len(copy) >0:
            result.add(copy)
    return result
    
def split_distance(initial1,initial2,initial): 
    if initial ==['']:
        return 0
    if len(initial1) == 0:
        count= abs(len(initial)-len(initial2))
    else:
        if len(initial2) == 0:
            count= abs(len(initial)-len(initial1))
        else:
            count= abs(len(initial)-len(initial1))+abs(len(initial)-len(initial2))
            
    return count
### testing
def parse(myfile):
    data = open(myfile,'r')
    data = data.readlines()
    tree = data[0]
    tree = Tree(tree,format = 1)
    dic  = {}
    for index in range(1,len(data)):
        array = data[index].split(':')
        name = array[0]
        info = array[1].split('(')
        if len(info) >1:
            gene_block = info[0].split(' ')[1]
            dic[name]=gene_block
        else:
            dic[name]=''

    for node in tree.traverse('postorder'):
        name = node.name
        node.add_features(gene_block = dic[name],distances=[0,0,0])
    return dic,tree

def calculate_distance(tree):
    for node in tree.traverse('postorder'):
        if not node.is_leaf():
            gene_block = TextFace(node.gene_block)
            children=[]
            for child in node.get_children():
                children.append(child)
            accumulate_distance0 = children[0].distances
            accumulate_distance1 = children[1].distances
    
            # calculate  deletion distance
            deletion,gene_set   = del_distance(children[0].gene_block,children[1].gene_block,node.gene_block)
            new_deletion = deletion+ accumulate_distance0[0]+accumulate_distance1[0]
            # calcualte  split distance
            # reduce by the gene_set
            initial1  = (children[0].gene_block).split('|')
            initial1  = reduction_gene(initial1,gene_set)
            initial2  = (children[1].gene_block).split('|')
            initial2  = reduction_gene(initial2,gene_set)
            initial   = (node.gene_block).split('|')
            split     = split_distance(initial1,initial2,initial)
            new_split = split + accumulate_distance0[2]+accumulate_distance1[2]
            # new_distance
            node.distances = [new_deletion,0,new_split]
    return tree
    
def display(tree):
    for node in tree.traverse('postorder'):
        string = node.name +': '+node.gene_block
        gene_block = TextFace(string)
        distances = TextFace(node.distances)
        if node.is_leaf():            
            node.add_face(gene_block, column=0, position = "branch-right")
        else:
            node.add_face(gene_block, column=0, position = "branch-top")
            node.add_face(distances, column=0, position = "branch-top")
    return tree
    
dic,tree = parse('sample_out.txt')
tree = calculate_distance(tree)
tree = display(tree)
tree_style = TreeStyle()
tree_style.show_leaf_name = False
tree_style.min_leaf_separation = 5
tree_style.extra_branch_line_type = 0
tree_style.draw_guiding_lines=True
tree_style.guiding_lines_type = 1
tree.render('output_check.pdf',dpi=1000,tree_style=tree_style)
tree.show(tree_style=tree_style)
