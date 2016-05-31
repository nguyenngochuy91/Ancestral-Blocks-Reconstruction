#!/usr/bin/env python
''' Author  : Huy Nguyen
    Program : Providing visualization of the Reconstruction. Grouping genomes
              into group color
    Start   : 05/08/2016
    End     : 05/08/2016
'''
from ete3 import *
import argparse
import os

def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--Operon","-o", help="Operon file name")
    parser.add_argument("--Group","-g", help="Color Grouping")
    args = parser.parse_args()
    return args
    
# color group of genome
# create a dictionary, sadly this is done manually
def parse(file):
    color_dic={}
    group= open(file,'r')
    for line in group.readlines():
        key= line.split(':')[0]
        color=line.split(':')[1]
        value = color.split('\n')[0]
        color_dic[key]=value
    return color_dic


if __name__ == "__main__":
    start = time.time()
    args = get_arguments()
    tree= Tree(args.Operon)
    # using the color dic to color group
    color_dic = parse(args.Group)
    for node in tree.iter_descendants("postorder"):
        if not node.is_leaf():
            # create face contain initial set info
            initial = TextFace(node.initial)
            node.add_face(initial, column=0, position = "branch-top")
        else:
            name = node.name.split('_')
            # modify name to be normal, and append the gene block info to it
            node.name = name[0]+' '+ name[1]+':' + node.gene_block
            # get accesion number
            short = name[2]+'_'+name[3]
            # get the color
            for key in color_dic:
                if short in key:
                    color = color_dic[key]
                    break
            nstyle = NodeStyle()
            nstyle["fgcolor"] = color
            nstyle["vt_line_color"]=color
            nstyle["hz_line_color"]=color
            node.set_style(nstyle)
    tree.show()



