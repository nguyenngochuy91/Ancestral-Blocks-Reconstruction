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
    parser.add_argument("--Operon","-i", help="Operon file name")
    parser.add_argument("--Group","-g", help="Color Grouping")
    parser.add_argument("--Image","-o", help="Output Image")
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
    far,distance = tree.get_farthest_leaf()
    for node in tree.iter_descendants("postorder"):
        if not node.is_leaf():
            child1,child2 = node.get_children()
            color1 = child1.node_color
            color2 = child2.node_color
            if color1 == color2 and color1 !='mixed':
                node.add_features(node_color=color1)
            else:
                node.add_features(node_color='mixed')
        else:
            name = node.name.split('_')
            # modify name to be normal, and append the gene block info to it
            node.name = name[0]+' '+ name[1]
            # get accesion number
            short = name[2]+'_'+name[3]
            # get the color
            if short in color_dic:
                color = color_dic[short]
            else:
                color = "blue"
            node.add_features(node_color=color)
            node.add_face(TextFace(node.name), column =0, position ="aligned")
            # node.dist = distance
        if node.node_color != 'mixed':
            nstyle = NodeStyle()
            nstyle["fgcolor"] = color
            # nstyle["vt_line_color"]=color
            # nstyle["hz_line_color"]=color
            node.set_style(nstyle)
    tree_style = TreeStyle()
    tree_style.show_leaf_name = False
    tree_style.min_leaf_separation = 5
    tree_style.extra_branch_line_type = 0
    tree_style.draw_guiding_lines=True
    tree_style.guiding_lines_type = 1                                        
    # render the image
    tree.render(args.Image+'.pdf',dpi=1000,tree_style=tree_style)
#    tree.show(tree_style=tree_style)




