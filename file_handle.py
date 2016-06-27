#!/usr/bin/env python
''' Author  : Huy Nguyen
    Program : Directory handling, file parsing and writing
    Start   : 05/04/2016
    End     : 05/05/2016
'''
import os
import argparse
import time
import uuid

###############################################################################
## Helper function to parse arguments, check directoring, ...
###############################################################################
# traverse and get the file
def traverseAll(path):
    res=[]
    for root,dirs,files in os.walk(path):
        for f in files:
            res.append(root+f)
    return res

# check if directory is readable
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

# get the arguments from command line
def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--InputDataDirectory","-i",action=readable_dir,help="This contain the translation result in term of alphabet")
    parser.add_argument("--OutputDirectory","-o", help="Output of this program will be stored in the path supplied here. It will make a new directory if path given is valid or it will raise an error")
    parser.add_argument("--TreeFile","-t", help="Tree file name")
    parser.add_argument("--Method","-m",help="Choose method to reconstruct (global or local)")
    args = parser.parse_args()
    return args

# check if the directory already exists, or cant create one

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
    
# function to write out the mapping between an alphabet letter and the gene of the operon
def mapping_write(mapping):
    myString= ''
    for key in mapping:
        myString+=mapping[key]+','+key+'\t'
    myString+='\n'
    return myString
