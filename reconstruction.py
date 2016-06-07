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
from findParent import *
from ete3 import *
# traverse and get the file
def traverseAll(path):
    res=[]
    for root,dirs,files in os.walk(path):
        for f in files:
            res.append(root+f)
    return res

class readable_dir(argparse.Action):
    def __call__(self,parser, namespace, values, option_string=None):
        prospective_dir=values
        if not os.path.isdir(prospective_dir):
           try:
               os.mkdir(prospective_dir)
           except OSError:
               print (argparse.ArgumentTypeError("readable_dir:{0} is not a readable dir".format(prospective_dir)))
        if os.access(prospective_dir, os.R_OK):
            setattr(namespace,self.dest,prospective_dir)
        else:
            raise argparse.ArgumentTypeError("readable_dir:{0} is not a readable dir".format(prospective_dir))

def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--InputDataDirectory","-i",action=readable_dir,help="This contain the translation result in term of alphabet")
    parser.add_argument("--OutputDirectory","-o", help="Output of this program will be stored in the path supplied here. It will make a new directory if path given is valid or it will raise an error")
    parser.add_argument("--TreeFile","-t", help="Tree file name")
    args = parser.parse_args()
    return args


def chk_output_directory_path(OutputDirectory,sessionID):
    if not os.path.exists(OutputDirectory + "_" + str(sessionID)):
        try:
           #os.mkdir(OutputDirectory + "_" + str(sessionID))
           return True
        except OSError:
           print ("Unable to create directory:", OutputDirectory)
           sys.exit()

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
    

def reconstruct(myfile,tree):
    tree= Tree(tree)
    result=parsing(myfile)
    genomes=result[1]
    mapping=result[0]
    # traverse the first time to assign name, and the node_type
    count =0
    for node in tree.iter_descendants("postorder"):
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
    
    # this serve to find which node has which children
    info = ''
    for node in tree.iter_descendants("postorder"):
        if not node.is_leaf():
            string = str(node.name)+':'
            for children in node.get_children():
                string += children.name+'\t'
            info += string +'\n'
    out= open('./structure.txt','w')
    out.write(info)
    out.close()
    # use findParent function to find geneblock
    for node in tree.iter_descendants("postorder"):
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
                                   )
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
                                   genome_list[0].deletion)
                tuple2=(genome_list[1].initial,genome_list[1].elementCount,
                                   genome_list[1].count,genome_list[1].dup,
                                   genome_list[1].deletion)
                # run the SS function
                mytuple=findSetInitial_SS(tuple1,tuple2)
            # at the end, from mytuple, insert the info into the internal node
            node.add_features(initial=mytuple[0],elementCount=mytuple[1],
                                  count = mytuple[2],dup=  mytuple[3],
                                  deletion = mytuple[4])
    # set mapping back to a string
    myString= ''
    for key in mapping:
        myString+=mapping[key]+','+key+'\t'
    myString+='\n'
    return (tree,myString)
            # print (node.name,node.initial,node.elementCount)
    # print (tree.write(features=['name','initial','gene_block']))
        
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
            myTuple=reconstruct(r,treeFile)
            tree = myTuple[0]
            mapping = myTuple[1]
            #if f == 'rplKAJL-rpoBC':
            #    for node in tree.iter_descendants("postorder"):
            #        if node.name == 'Node8' or node.name == 'Node15' or node.name =='Node18':
            #            print node.name,node.initial,node.elementCount,node.count
            outfile=open(outputsession+'/'+f+'_mapping','w')
            outfile.write(mapping)      
            outfile.close()
            tree.write(format=2, outfile=outputsession+'/'+f,features=['name','initial','gene_block','deletion'])
    print (time.time() - start)
