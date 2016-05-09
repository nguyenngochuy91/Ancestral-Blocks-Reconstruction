#!/usr/bin/env python
''' Author  : Huy Nguyen
    Program : Providing visualization of the Reconstruction. Grouping genomes
              into group color
    Start   : 05/08/2016
    End     : 05/08/2016
'''
from ete3 import *

# color group of genome
# create a dictionary, sadly this is done manually
alpha=('Cau','Bra','Sin','Agr','Mes','Bru')
beta=('Hel','Cam','Wol')
gamma=('Bor','Ral','Nit','Chr','Nei')
delta=('Xan','Xyl','She','Pse','Vib','Can','Pas','Hae','Buc','Yer','Sal','Esc')
# assign color
color_dic={}
color_dic[alpha]='cyan'
color_dic[beta]='black'
color_dic[gamma]='hotpink'
color_dic[delta]='brown'


def main():
    f= input("Which operon you want to show the tree:")
    tree= Tree('./reconstruction/'+f)
    # using the color dic to color group
    for node in tree.iter_descendants("postorder"):
        if node.is_leaf():
            # get the first 2 letter of the node name
            # print (node.name)
            name = node.name.split('_')
            # it is the first 3 index of the index 0 of name
            short = name[0][:3]
            # get the color
            for key in color_dic:
                if short in key:
                    color = color_dic[key]
                    break
            nstyle = NodeStyle()
            nstyle["fgcolor"] = color
            node.set_style(nstyle)
        else:
            # create face contain initial set info
            initial = TextFace(node.initial)
            node.add_face(initial, column=0, position = "branch-top")
            
    tree.show()

main()
