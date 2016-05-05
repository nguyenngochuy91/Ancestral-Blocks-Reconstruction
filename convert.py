#!/usr/bin/env python
''' Author  : Huy Nguyen, and David C.Ream
    Program : Given the operon directory, for each operons file, get the info about the gene in each genomes, map it
              into alphabet letter , get the gap, map gap to '|' and write to ouputfile.
    Start   : 05/04/2016
    End     : /2016
'''

import os
import argparse
import time
import uuid
# traverse and get the file
def traverseAll(path):
    res=[]
    for root,dirs,files in os.walk(path):
        for f in files:
            res.append(root+f)
    return res


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--OperonDataDirectory","-n",action=readable_dir,help="This directory should contain files with gene name, start, stop, strand direction information for each genome.")
    parser.add_argument("--splitDistance","-d", type=int ,default=500,help="Splitting distance")
    parser.add_argument("--OutputDirectory","-o", help="Output of this program will be stored in the path supplied here. It will make a new directory if path given is valid or it will raise an error")
    args = parser.parse_args()
    return args


def chk_output_directory_path(OutputDirectory,sessionID):
    if not os.path.exists(OutputDirectory + "_" + str(sessionID)):
        try:
           #os.mkdir(OutputDirectory + "_" + str(sessionID))
           return True
        except OSError:
           print "Unable to create directory:", OutputDirectory
           sys.exit()
# given a dic of info of genome genes (info has position and strand), provide an string of input
# with gap as '|'

def toString(gene_list):
    
def convert(file):
    infile = open(file,'r')
    mapping =''
    body=''
    dic_map ={}
    dic_distance={}
    dic_distance['+1']={}
    for line in infile.readlines():
        if line[0] == '(':
            mapping = line.split('\t')[:-1] # 
            for item in mapping:
                item_split = item.split(',')
                key = item_split[0]
                value = item_split[1]
                dic_map[key]=value # {'astA': 'a'}
        else:
            string =''
            genome = line.split(':')[0] # (line.split(':')= ['NC_002696', '(astA,634700,635744,1)\t(astD,635730,637149,1)\t(astB,637145,638426,1)\t(astC,638435,639614,1)\t']
            genes_string = line.split(':')[1]
            # to deal with each genes, they are in tuple, where first is the name of the gene, follow by the position, and the strand it is on
            # should consider 2 type of strand (so i can put a gap correctly
            genes_string = genes_string .split('\t')[:-1] # ['(astA,634700,635744,1)', '(astD,635730,637149,1)', '(astB,637145,638426,1)', '(astC,638435,639614,1)']
            genes_string = list(set(genes_string))
            
            

if __name__ == "__main__":

    start = time.time()
    args = get_arguments()
    sessionID = uuid.uuid1()
    condition = chk_output_directory_path(args.OutputDirectory,sessionID)
    if condition:
        outputsession = args.OutputDirectory
        os.mkdir(outputsession)
        res = traverseAll(args.OperonDataDirectory)
        for r in res:
            root,f = os.path.split(r)
            outfile = open(outputsession+'/'+f,'w')
            
