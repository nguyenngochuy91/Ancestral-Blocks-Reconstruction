#!/usr/bin/python

from homolog4 import Homolog

import os

import argparse

def parser_code():

    parser = argparse.ArgumentParser()
                     
    
    parser.add_argument("--input","-i", default="./optimized_gene_block/",help="optimized gene block ")  
    parser.add_argument("--gene_name", "-g", default='gene_block_names_and_genes.txt',
                help="the gene_block_names_and_genes that stores the name of the operon and its genes")            
    parser.add_argument("--output","-o", default="result/",
                help="where the result be stored (result/)")
    parser.add_argument("-a", "--accession",  default='tree/accession_to_common.csv',
                help="Filter file, default as potential file, if NONE then not doing the parse filter")  
                            
    return parser.parse_args()


def get_accession(accession):
    dict = {}
    for line in open(accession,'r').readlines():
        line = line.strip().split(',')
        dict[line[0]]= line[1]
    return dict
        
## parse the gene block names and genes txt file
def parse(operon_genes_dict):
    result = {}
    infile = open(operon_genes_dict,'r')
    for line in infile.readlines():
        line = line.strip().split()
        result[line[0]] = line[1:]
    return result
## Traverses the genome information directory
def traverseAll(path):
    res=[]
    for root,dirs,files in os.walk(path):
        for f in files:
            res.append(root+'/'+f)
    return res    
    
## given an operon file (astCADBE.txt), format the info into format easier to read 
def formatOperon(operon,output,operon_genes_dict,accession_dict):
    alphabet = 'abcdefghijklmnop'
    operon_name = operon.split('/')[-1].split('.')[0]
    genes = sorted(operon_genes_dict[operon_name])
    outfile = open(output+operon_name,'w')
    for i in range(len(genes)):
        outfile.write(genes[i]+','+alphabet[i]+'\t')
    outfile.write('\n')
    result = {}
    for line in [i.strip() for i in open(operon).readlines()]:

            hlog = Homolog.from_blast(line)
            accession = hlog.accession()[:-2]
            start = str(hlog.start())
            end = str(hlog.stop())
            strand = str(hlog.strand())
            gene_name = hlog.blast_annotation()
            if accession in result:
                result[accession].append([gene_name, start, end, strand])
            else:
                result[accession]=[[gene_name, start, end, strand]]

    for species in accession_dict:
        outfile.write(species+':')
        if species in result:
            for item in result[species]:
                outfile.write(','.join(item)+'\t')
        outfile.write('\n')
    outfile.close()

if __name__ == "__main__":
    args                  = parser_code()
    input                 = args.input
    result                = args.output
    accession             = args.accession
    operon_genes_dict     = args.gene_name
    operon_genes_dict     = parse(operon_genes_dict)
    accession_dict        = get_accession(accession)
    try:
        os.makedirs(result)
    except:
        print ("Result dic already exists")
    # goes through all the file name in the optimized_gene_block dic
    res = traverseAll(input)
    for file in res:
        formatOperon(file,result,operon_genes_dict,accession_dict)
    
    