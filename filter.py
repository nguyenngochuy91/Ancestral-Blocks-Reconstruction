#!/usr/bin/env python
''' Author  : Huy Nguyen
    Program : filter the all genbank file , only keep those that have 1 copy of chromosomes 
                filter our plasmid files
    Start   : 09/08/2016
    End     : 09/1/2016
'''
import os
import argparse
import shutil

## Traverses the genome information directory
def traverseAll(path):
    res=[]
    for root,dirs,files in os.walk(path):
        for f in files:
            res.append(root+'/'+f)
    return res   

## parse argument
def parse_args():
    parser = argparse.ArgumentParser()
                     
    parser.add_argument("--genomes_directory","-g", help="The directory that store all the genomes file (genomes_folder)")                   
    return parser.parse_args()
    
if __name__ == '__main__':
    args = parse_args()
    genomes_directory = args.genomes_directory
    res = traverseAll(genomes_directory)
    dic = {}
    for file in res:
        species_name = file.split('/')[-2]
        if species_name in dic and dic[species_name] == 2:
            continue
        infile = open(file,'r')
        infile.readline()
        # get the info in the secondline
        if "plasmid" in infile.readline():
            os.remove(file)
        else:
            if species_name not in dic:
                dic[species_name] = 1
            else:
                dic[species_name] = 2
                shutil.rmtree('/'.join(file.split('/')[:-1]))