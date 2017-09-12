#!/usr/bin/python

from multiprocessing import Pool
import time
import os
import sys
import argparse
from homolog4 import *
from collections import defaultdict

# Copyright(C) 2015 David Ream
# Released under GPL version 3 licence. http://www.gnu.org/licenses/lgpl.html
# Do not remove this comment

# This exists to  make the main function easier to read. It contains code to run the argument parser, and does nothing else.
def parser_code():

    parser = argparse.ArgumentParser(description="Parse the results of a BLAST -m8 search and organize the results by specific gene blocks. The program will save the results in a directory designated by the user, or the default './blast_parse/'.")
                
    parser.add_argument("-i", "--infolder", dest="infolder", default='./blast_result/', metavar="DIRECTORY",
                help="A file that contains the path to every organism database that you are interested in.")
    
    parser.add_argument("-o", "--outfolder", dest="outfolder", metavar="DIRECTORY", default='./blast_parse/',
                help="Folder where the BLAST results will be stored. Default is the folder './blast_result/'.")
    
    # rename this option to gene block file... or something more accurately conveying the idea that this is the file
    # that is composed: gene_block_name: gene1, gene2, gene3... etc all fields tab delenated.
    parser.add_argument("-b", "--gene_block_query", dest="gene_block_query", default='gene_block_names_and_genes.txt', metavar="FILE",
                help="A file that contains the names and genes comprising the gene blocks that are under investigation.")

    parser.add_argument("-f", "--filter", dest="filter", default='', metavar="FILE",
                help="A file that contains the accession numbers of the organisms that are under investigation.")            
    
    parser.add_argument("-n", "--num_proc", dest="num_proc", metavar="INT", default = os.sysconf("SC_NPROCESSORS_CONF"), type=int,
                help="Number of processors that you want this script to run on. The default is every CPU that the system has.")
                
    parser.add_argument("-q", "--quiet", dest="quiet", action="store_true", default=False,
                help="Suppresses most program text outputs.")
                            
    return parser.parse_args()


def check_options(parsed_args):
    if os.path.isdir(parsed_args.infolder):
        infolder = parsed_args.infolder
    else:
        print "The infolder directory %s does not exist." % parsed_args.infolder
        sys.exit()
    
    # if the directory that the user specifies does not exist, then the program makes it for them. 
    if not os.path.isdir(parsed_args.outfolder):
        os.makedirs(parsed_args.outfolder)
    outfolder = parsed_args.outfolder
    if outfolder[-1] != '/':
        outfolder = outfolder + '/'
        
    if os.path.exists(parsed_args.gene_block_query):
        gene_block_query = parsed_args.gene_block_query
    else:
        print "The gene block query file %s does not exist." % parsed_args.gene_block_query
        sys.exit()
    
    
    if os.path.exists(parsed_args.filter):
        filter_file = parsed_args.filter
    elif parsed_args.filter == '':
        filter_file = parsed_args.filter
    else:
        print "The filter file %s does not exist." % parsed_args.filter
        sys.exit()
        
    # section of code that deals determining the number of CPU cores that will be used by the program
    if parsed_args.num_proc > os.sysconf("SC_NPROCESSORS_CONF"):
        num_proc = os.sysconf("SC_NPROCESSORS_CONF")
    elif parsed_args.num_proc < 1:
        num_proc = 1
    else:
        num_proc = int(parsed_args.num_proc)
    
    quiet = parsed_args.quiet    

    return infolder, outfolder, gene_block_query, filter_file, num_proc, quiet


#this function will return all of the files that are in a directory. os.walk is recursive traversal.
def returnRecursiveDirFiles(root_dir):
    result = []
    for path, dir_name, flist in os.walk(root_dir):
        for f in flist:
            fname = os.path.join(path, f)
            if os.path.isfile(fname):
                result.append(fname)
    return result


# The filter file here i think should be either a vacant value (such as '') or a user defined
# value.  I do not think by default it should be given, like I have made the default behavior.    
def parallel_blast_parse_dict(in_folder, out_folder, num_proc, filter_file, gene_block_dict):
    result = {}
    gene_block_out_folder = out_folder
    if not os.path.isdir(gene_block_out_folder):
        os.makedirs(gene_block_out_folder)
    if filter_file != '':
        tmp = returnRecursiveDirFiles(in_folder)
        nc_list = [i.strip() for i in open(filter_file).readlines()]
        fname_list = [i for i in tmp if i.split('/')[-1].split('.')[0] in nc_list]
    else:
        fname_list = returnRecursiveDirFiles(in_folder)
        
    # this is working correctly...
    #print "blast_parse fname_list", fname_list
    for fname in fname_list:
        for line in [i.strip() for i in open(fname).readlines()]:
            try:
                hlog = Homolog.from_blast(line)
            except:
                print "Error in function parallel_blast_parse_dict from script blast_parse.py, conversion from result to Homolog class failed.", line
            #hlog.Print()
            
            # this might have to be changed.... 
            #accession = hlog.accession()
            try:
                accession = hlog.accession()
            except:
                print "There was an error in the Homolog class, found in line function parallel_blast_parse_dict from script blast_parse.py"
            
            
            #predicted_gene = hlog.blast_annatation()
            predicted_gene = hlog.blast_annotation()
            
            '''
            try: # faster implementation than "if predicted_gene in gene_block_dict.keys():"
                gene_block = gene_block_dict[predicted_gene]
                # Debugging the missing gene block casABCDE12... no idea right now.
                if gene_block == 'casABCDE12':
                    print 'AFDFAFDSF'
                if gene_block in result.keys():    # check if the gene block is in the result dict already, if not make a new entry in the else clause
                    if accession in result[gene_block].keys(): # Check if the organism has been added to the gene block
                        result[gene_block][accession].append(hlog.ret_str())
                    else: # if the organim has not been added to the rest yet, add it
                        result[gene_block].update({accession:[hlog.ret_str()]})
                else: # add the gene block to the result
                    result.update({gene_block: {accession: [hlog.ret_str()]}})
            except:
                pass
                
            '''
            try: # faster implementation than "if predicted_gene in gene_block_dict.keys():"
                gene_block = gene_block_dict[predicted_gene]
                # Debugging the missing gene block casABCDE12... no idea right now.
                # ok, it's there, so omitting this, leaving the comment in for the time being though
                #if gene_block == 'casABCDE12':
                #    print 'AFDFAFDSF'
                if gene_block in result.keys():    # check if the gene block is in the result dict already, if not make a new entry in the else clause
                    if accession in result[gene_block].keys(): # Check if the organism has been added to the gene block
                        result[gene_block][accession].append(hlog.to_file())
                    else: # if the organims is not part of the gene block, add it
                        result[gene_block].update({accession:[hlog.to_file()]})
                else: # add the gene_block to the result
                    result.update({gene_block: {accession: [hlog.to_file()]}})
            except:
                pass
            
    #print sorted(result.keys()), len(result)
    
    '''
    # For the time being, i am going to cause each intermediate step in this pipeline to save in a folder called intermediate_for_debug
    intermediate_folder = './intermediate_for_debug/'
    if not os.path.isdir(intermediate_folder):
        os.makedirs(intermediate_folder)
    
    # For this step i will save the result in 'unfiltered_gene_block/'
    
    unfilter_folder = 'unfiltered_gene_block/'
    if not os.path.isdir(intermediate_folder + unfilter_folder):
        os.makedirs(intermediate_folder + unfilter_folder)
    '''    
    for gene_block in result.keys():
        
        # this code is omitted because it used to debugging purposes, and is currently unneeded
        '''
    	outfile = intermediate_folder + unfilter_folder + gene_block + '.txt'
        #print "outfile", outfile
        handle = open(outfile, 'w')
        
        for accession in result[gene_block].keys():
        	handle.write('\n'.join(result[gene_block][accession]) + '\n')
        handle.close()
        '''
        
        # save results where i actually want them to go:
        #print "plast_parse.py outfile", out_folder + gene_block + '.txt'
        handle = open(out_folder + gene_block + '.txt', 'w')
        for accession in result[gene_block].keys():
        	handle.write('\n'.join(result[gene_block][accession]) + '\n')
        handle.close()
        
        
# I have to figure out a better name for this function. The jist of what I am doing here is as follows:
# First, I will provide a file name that contains all of the hits for every organism that we are interested in.
# Then it sorts this homolog list first by organism, then by locus. (By the required input, the files already have 
# been screened for both eval cutoff and gene block membership. The function will then return a dictionary for that
# gene block. The intention is that this data structure will then be used to find the best hit for the locus out of the
# many redundant hits, however this functionality will be handled another function that i have yet to write/test.
def return_gene_block_list(fname):
    gene_block = fname.split('/')[-1].split('.')[0]
    hlog_list = [Homolog.from_file(i.strip()) for i in open(fname).readlines()]
    result_dict = {}
    for hlog in hlog_list:
        accession = hlog.accession()
        locus = hlog.locus()
        if accession not in result_dict.keys():
            result_dict.update({accession: {}})
        if locus not in result_dict[accession].keys():
            result_dict[accession].update({locus: [hlog]})
        else:
            result_dict[accession][locus].append(hlog)
    #print result_dict[accession]
    
    return result_dict
        
# currently not used
'''
# might  not use this: will see
def parallel_return_gene_block_list(infolder, outfolder, num_proc):
    pool = Pool(processes = num_proc)
    organism_dict_for_recovery = dict(pool.map(parallel_gene_block_fasta, genome_of_interest_list))
'''

# This function will take the organism-locus dict (per gene block file) and determine the best homolog.
def best_homolog_list(gene_block_dict, outfile):
    result = []
    
    for org in sorted(gene_block_dict.keys()):
        for locus in sorted(gene_block_dict[org].keys()):
            hlog_list = gene_block_dict[org][locus][1:]
            best_hit = gene_block_dict[org][locus][0]
            gene_count = defaultdict(int) # i am goign to use this, to see if a locus has more than one predicted gene, and the count ratio
            gene_count[best_hit.predicted_gene()] +=1
            for hlog in hlog_list:
                gene_count[hlog.predicted_gene()] +=1
                if best_hit.e_val() > hlog.e_val():
                    best_hit = hlog
            #print gene_count.keys()
            result.append(best_hit)
    handle = open(outfile, 'w')
    handle.write('\n'.join([i.ret_str() for i in result]))
    handle.close()


def return_gene_to_gene_block_dict(fname):
    gene_block_dict = {}
    for line in [i.strip().split('\t') for i in open(fname).readlines()]:
        gene_block = line[0]
        for gene in line[1:]:
            gene_block_dict.update({gene: gene_block})
    return gene_block_dict
 
def main():
    start = time.time()
    
    parsed_args = parser_code()
    
    infolder, outfolder, gene_block_query, filter_file, num_proc, quiet = check_options(parsed_args)
    
    if not quiet:
        print infolder, outfolder, gene_block_query, filter_file, num_proc, quiet
    
    # This code makes a dictionary mapping gene annotation to the gene block that it belong to
    gene_block_dict = return_gene_to_gene_block_dict(gene_block_query)
    #print "gene_block_dict", gene_block_dict
    
    #parallel_blast_parse_dict('./blast_parse/organism_raw_info/', './blast_parse/filtered_homologs/', num_proc, './genbank_pathway_lists/nc_filter_file.txt', gene_block_dict)
    parallel_blast_parse_dict(infolder, outfolder, num_proc, filter_file, gene_block_dict)
    
    
    #gene_block_dict = return_gene_block_list('./blast_parse/filtered_homologs/atpIBEFHAGDC.txt')
    
    #best_homolog_list(gene_block_dict, './blast_parse/processed_gene_block_files/atpIBEFHAGDC.txt')
    
    if not quiet:
        print time.time() - start

    # ./blast_parse.py -f phylo_order.txt
if __name__ == '__main__':
    main()
