#!/usr/bin/env python
''' Author  : Huy Nguyen
    Program : Reconstruct the ancestral gene blocks given the newick tree and the leaf gene blocks
    Start   : 05/04/2016
    End     : 05/05/2016
'''
import os
import argparse
import time
import uuid
from findParent_local import *
from ete3 import *
from file_handle import *

###############################################################################
## Reconstruct method
###############################################################################
'''@function: Reconstruct the newick tree file with gene block info for inner 
              node.
   @input   : tree in nwk format,and a mapping between alphabet and gene name
   @output  : tree in nwk format,gene g and a string of the info
'''
def reconstruct(genomes,tree):
    # traverse the first time to assign name, and the node_type
    count =0
    for node in tree.traverse('postorder'):
        if node.name == '':
            count +=1
            node.name = 'Node'+ ' ' + str(count)
            leaf = 0
            for children in node.get_children():
                if children.is_leaf():
                    leaf+=1
                # create feature to apply which function of 3 functions
                if leaf == 1:                        
                    node.add_features(node_type= 'SG')
                if leaf == 0:                        
                    node.add_features(node_type= 'SS')
                if leaf == 2:
                    node.add_features(node_type= 'GG')
        else:
            mylist= node.name.split('_')
            name = mylist[2]+'_'+mylist[3]
            # print(genomes[name])
            node.add_features(gene_block=genomes[name])   
     
    # use findParent function to find geneblock
    for node in tree.traverse('postorder'):
        if not node.is_leaf():
            if node.node_type == 'GG':
                # get the info from the 2 genomes
                genome_list=[]
                for children in node.get_children():
                    genome_list.append(children.gene_block)
                split1= countSplit(genome_list[0])
                split2=countSplit(genome_list[1])
                # run the GG function
                if split1 == 0:
                    mytuple= findSetInitial_GG(genome_list[0],genome_list[1],split1,
                                               split2)
                else:
                    mytuple= findSetInitial_GG(genome_list[1],genome_list[0],split2,
                                               split1)
                               
            if node.node_type == 'SG':            
                for children in node.get_children():
                    if children.is_leaf():
                        genome=children.gene_block # get the genome info
                    else: # it is a set, extract the info into a tuple
                        initial_tuple = (children.initial,children.elementCount,
                                   children.count,children.dup,children.deletion,
                                   children.duplication,children.split)
                split1= countSplit(genome)
                # run the SG function
                mytuple=findSetInitial_SG(initial_tuple,genome,split1)
                
            if node.node_type == 'SS': 
                genome_list=[]
                for children in node.get_children():
                    # append the whole node first to the genome_list
                    genome_list.append(children)
                tuple1=(genome_list[0].initial,genome_list[0].elementCount,
                                   genome_list[0].count,genome_list[0].dup,
                                   genome_list[0].deletion,genome_list[0].duplication,
                                   genome_list[0].split)
                tuple2=(genome_list[1].initial,genome_list[1].elementCount,
                                   genome_list[1].count,genome_list[1].dup,
                                   genome_list[1].deletion,genome_list[1].duplication,
                                   genome_list[1].split)
                # run the SS function
                mytuple=findSetInitial_SS(tuple1,tuple2)
            # at the end, from mytuple, insert the info into the internal node
            node.add_features(initial=mytuple[0],elementCount=mytuple[1],
                                  count = mytuple[2],dup=  mytuple[3],
                                  deletion = mytuple[4], duplication=mytuple[5],
                                  split = mytuple[6])
                                       
    return tree
            # print (node.name,node.initial,node.elementCount)
    # print (tree.write(features=['name','initial','gene_block']))
    
###############################################################################
## Main function 
###############################################################################
'''@function: Go through each operon file in the target directory,
              write out newick treefile with reconstrution info, and a 
              mapping file for gane name and alphabet.
   @input   : tree in nwk format,and a mapping between alphabet and gene name
   @output  : tree in nwk format,gene g and a string of the info
''' 
if __name__ == "__main__":

    start = time.time()
    args = get_arguments()
    sessionID = uuid.uuid1()
    condition = chk_output_directory_path(args.OutputDirectory,sessionID)
    treeFile=args.TreeFile
    if condition:
        outputsession = args.OutputDirectory
        os.mkdir(outputsession)
        res = traverseAll(args.InputDataDirectory)
        for r in res:
            # check for DS.Store
            root,f = os.path.split(r)
            if "DS_Store" in f:
                continue
            mapping,genomes = parsing(r)
            mapping = mapping_write(mapping) # modify the string mapping to write out
            tree= Tree(treeFile) # using ete3 to read treeFile
            tree = reconstruct(genomes,tree)
            #if f == 'rplKAJL-rpoBC':
            #    for node in tree.iter_descendants("postorder"):
            #        if node.name == 'Node8' or node.name == 'Node15' or node.name =='Node18':
            #            print node.name,node.initial,node.elementCount,node.count
            outfile=open(outputsession+'/'+f+'_mapping','w')
            outfile.write(mapping)      
            outfile.close()
            tree.write(format=2, outfile=outputsession+'/'+f,features=['name',
            'initial','gene_block','deletion','duplication','split'])
    print (time.time() - start)
