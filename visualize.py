#!/usr/bin/env python
"""
Created on Wed Jul  6 16:14:23 2016

@author: huyn
@purpose: visualize professor John output
"""
from ete3 import *
from findParent_local import *
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
            distances  = info[1].split(')')[0]
            distances  = distances.split(',')
            deletion = int(distances[0])
            duplication = int(distances[1])
            split = int(distances[2])
            distances = [deletion,duplication,split]
            dic[name] ={}
            dic[name]['gene_block']=gene_block
            dic[name]['distances'] = distances
        else:
            dic[name]= {}
            dic[name]['gene_block']=''
            dic[name]['distances']=[]
    
    for node in tree.traverse('postorder'):
        name = node.name
        node.add_features(gene_block = dic[name]['gene_block'],distances=dic[name]['distances'])
    return dic,tree

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
# print dic
tree = display(tree)
tree_style = TreeStyle()
tree_style.show_leaf_name = False
tree_style.min_leaf_separation = 5
tree_style.extra_branch_line_type = 0
tree_style.draw_guiding_lines=True
tree_style.guiding_lines_type = 1
tree.render('output_prof.pdf',dpi=1000,tree_style=tree_style)
tree.show(tree_style=tree_style)
