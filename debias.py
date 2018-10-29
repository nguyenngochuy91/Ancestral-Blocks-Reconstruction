#!/usr/bin/env python
''' Author  : Huy Nguyen
    Program : debias
    Start   : 06/04/2017
    End     : 06/08/2017
'''
import os
import argparse
import shutil
from Bio import SeqIO
from ete3 import Tree
def parser_code():

    parser = argparse.ArgumentParser(description='The purpose of this script to debias tree based on parameter')
    
    parser.add_argument("-i", "--input_tree", help="Input tree that we want to debias")
    
    parser.add_argument("-o", "--pda_out", help="Output of pda to be store.")
                
    parser.add_argument("-s", "--tree_size", help="Reduce the size of the tree to this size, for example, you can reduce your number of species from 100 to 30 by input 30")
 
    parser.add_argument("-r", "--ref", help="Force to include the following species, here I force to include the reference species",default = None)                                       
    return parser.parse_args()
    
def parse_pda(handle):
    for line in handle.readlines():
        if len(line) > 100:
            return line.strip()

        
    
if __name__ == "__main__":
    args  = parser_code()
    input_tree = args.input_tree
    pda_out = args.pda_out
    difference = '/'.join(input_tree.split('/')[:-1])+"/not_included.txt"
    t = Tree(input_tree, format = 1)
    leaves_full = t.get_leaf_names()
    size = int(args.tree_size)
    ref = args.ref
    for item in leaves_full:
        if ref in item:
            ref = item
    outfile = open("keep.txt","w")
    outfile.write(ref+'\n')
    outfile.close()
    cmd1 = "pda -k {} {} {} -if keep.txt".format(size,input_tree,pda_out)
    os.system(cmd1)
#    tree = parse_pda(open(pda_out,"r"))
#    t = Tree(tree,format = 1)
#    leaves_partial = t.get_leaf_names()
#    outfile = open("phylo_order.txt","w")
#    
#    for name in leaves_partial:
#        name = '_'.join(name.split('_')[-2:])+'\n'
#        outfile.write(name)
#    outfile.close()
