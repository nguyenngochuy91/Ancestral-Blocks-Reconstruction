#!/usr/bin/env python
''' Author  : Huy Nguyen
    Program : Reconstruct the ancestral gene blocks given the newick tree and the leaf gene blocks
    Start   : 05/04/2016
    End     : 05/05/2016
'''

from findParent import *
from ete3 import *

# function to parse each file in the directory and try to assign the genes block for each genome
# this return a tuple, in this tuple contains mapping gene to alphabet, and the genomes and its gene blocks
def parsing(file):
        mapping ={}
        genomes={}
        myfile = open(file,'r')
        for line in myfile.readlines():
                if line[0]!='N':
                        mylist= line.split('\t')[:-1]
                        for item in mylist:
                                tupple= item.split(',')
                                mapping[tupple[1]]=tupple[0]
                        # print (mapping)
                else:
                        item = line.split(':')
                        name = item[0]
                        gene_blocks= item[1].split('\n')[0]
                        genomes[name]=gene_blocks
        return (mapping,genomes)
result=parsing('./new_result/astCADBE')
genomes=result[1]
mapping=result[0]
tree= Tree('muscle.ph')
# print (tree)
count =0
for node in tree.iter_descendants("postorder"):
	if node.name == '':
		count +=1
		node.name = 'Node'+ str(count)
		# include new feature to create gene block
		node.add_feature('gene_block','')
		leaf = 0
		for children in node.get_children():
                        if children.is_leaf():
                                leaf+=1
                # create feature to know whether to apply which function
                if leaf == 1:                        
                        node.add_feature('node_type', 'SG')
                if leaf == 0:                        
                        node.add_feature('node_type' , 'GG')
                if leaf == 2:
                        node.add_feature('node_type' , 'SS')
        else:
                node.add_feature('gene_block',)   

# this serve to find which node has which children
for node in tree.iter_descendants("postorder"):
	if not node.is_leaf():
		string = str(node.name)+':'
		for children in node.get_children():
			string += children.name+'\t'
		print (string)
print (tree.write(features=['gene_block','node_type']))
                            

