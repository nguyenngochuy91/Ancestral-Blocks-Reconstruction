#!/usr/bin/env python
"""
Created on Wed Jul  6 16:14:23 2016

@author: huyn
@purpose: testing professor John output
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
            gene_block = info[0].split(' ')[0]
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
    # print dic
parse('sample_out.txt')