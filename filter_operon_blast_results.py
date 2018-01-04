#!/usr/bin/env python

from multiprocessing import Pool
import time
import os
import sys
import argparse
import math
from homolog4 import *
from collections import defaultdict
import itertools

# Copyright(C) 2015 David Ream
# Released under GPL version 3 licence. http://www.gnu.org/licenses/lgpl.html
# Do not remove this comment

# This exists to  make the main function easier to read. It contains code to run the argument parser, and does nothing else.
def parser_code():

    parser = argparse.ArgumentParser(description="This program will be used to remove spurious results from a BLAST search organized by gene block.")
    
    parser.add_argument("-i", "--infolder", dest="infolder", default='./blast_parse/', metavar="DIRECTORY",
                help="A folder that contains the gene block BLAST results.")
    
    parser.add_argument("-o", "--outfolder", dest="outfolder", metavar="DIRECTORY", default='./optimized_gene_block/',
                help="Folder where the filtered results will be stored. Default is the folder './optimized_gene_block/'.")

    parser.add_argument("-f", "--filter", dest="filter", default='', metavar="FILE",
                help="A file that contains the gene blocks that are under investigation.  All others will be omitted from analysis an results.")            
    
    parser.add_argument("-n", "--num_proc", dest="num_proc", metavar="INT", default = os.sysconf("SC_NPROCESSORS_CONF"), type=int,
                help="Number of processors that you want this script to run on. The default is every CPU that the system has.")
                
    parser.add_argument("-e", "--eval", dest="eval", default='1e-10', metavar="FLOAT", type=float,
                help="Use this option to change the eval for the BLAST search that is permitted. Useful if you would like to investigate what altering the eval threshold will do to your results.")
                
    parser.add_argument("-g", "--max_gap", dest="max_gap", metavar="INT", default = 500, type=int,
                help="Largest allowable gap to be considered a gene block by the analysis.")
                
    parser.add_argument("-q", "--quiet", dest="quiet", action="store_true", default=False,
                help="Suppresses most program text outputs.")
    
    return parser.parse_args()


def check_options(parsed_args):
    if os.path.isdir(parsed_args.infolder):
        infolder = parsed_args.infolder
    else:
        print("The infolder directory %s does not exist." % parsed_args.infolder)
        sys.exit()
    
    # if the directory that the user specifies does not exist, then the program makes it for them. 
    if not os.path.isdir(parsed_args.outfolder):
        os.makedirs(parsed_args.outfolder)
    outfolder = parsed_args.outfolder
    if outfolder[-1] != '/':
        outfolder = outfolder + '/'
     
    if os.path.exists(parsed_args.filter):
        filter_file = parsed_args.filter
    elif parsed_args.filter == '':
        filter_file = parsed_args.filter
    else:
        print("The filter file %s does not exist." % parsed_args.filter)
        sys.exit()
        
    # section of code that deals determining the number of CPU cores that will be used by the program
    if parsed_args.num_proc > os.sysconf("SC_NPROCESSORS_CONF"):
        num_proc = os.sysconf("SC_NPROCESSORS_CONF")
    elif parsed_args.num_proc < 1:
        num_proc = 1
    else:
        num_proc = int(parsed_args.num_proc)
    
    # validate the input for the eval
    try:    
        e_val = float(parsed_args.eval)
    except:
        print("The e-value you entered is not a floating point number, please enter a floating point number, ex. '1e-3', or '12'.")
        sys.exit()
        
    # validate the input for the maximum allowed gap
    try:    
        max_gap = int(parsed_args.max_gap)
        if max_gap <= 0:
           print("The gap that you entered %s is a negative number, please enter a positive integer." % parsed_args.max_gap)
           sys.exit()
        else:
           pass
    except:
        print("The gap that you entered %s is not an integer, please enter a positive integer." % parsed_args.max_gap)
        sys.exit()
    
    quiet = parsed_args.quiet
    
    return infolder, outfolder, filter_file, num_proc, e_val, max_gap, quiet


#this function will return all of the files that are in a directory. os.walk is recursive traversal.
def return_recursive_dir_files(root_dir):
    result = []
    for path, dir_name, flist in os.walk(root_dir):
        for f in flist:
            fname = os.path.join(path, f)
            if os.path.isfile(fname):
                result.append(fname)
    return result

def return_file_list(infolder, filter_file):
    if filter_file == '':
        return return_recursive_dir_files(infolder)   
    else:
        filter_list = [i.strip() for i in open(filter_file)]
        return [i for i in return_recursive_dir_files(infolder) if os.path.basename(i).split('.')[0] in filter_list]

# This function will take a BLAST tabular result, and remove any hits that are worse than the eval threshold provided
# The return will be a list of those hits as homolog objects.
def filter_eval(fname, e_val):
    # make a list of homolog class objects
    h_list = [Homolog.from_blast(i) for i in open(fname).readlines()]
    
    result = list(filter((lambda x: x.e_val() <= e_val), h_list))

    return result

def resolve_multiple_ORF_hits(hlist):
    result = [hlist[0]]
    
    curr_homolog = hlist[0]
    for next_homolog in hlist[1:]:
        # we have multiple BLAST hits that share a ORF, resolve using e_val
        if curr_homolog.start() == next_homolog.start():
            # If current homolog has better eval 
            if curr_homolog.e_val() <= next_homolog.e_val():
                # print curr_homolog.organism(), curr_homolog.locus(), "is duplicated"
                pass
            # The current homolog has a worse eval, remove for the better example
            else:
                result.pop(-1)
                result.append(next_homolog)
                #print "This totally worked", next_homolog.organism()
        else:
            result.append(next_homolog)
        # Now that we are done testing the current and next homolog against    
        curr_homolog = next_homolog

    return result

# The purpose of this function is to take a list of homologs, that have been e_val (or potenitally another means) filtered.
# The return is all homologs from organisms that contain at least one neighborhood defined by max_gap.
def return_valid_organism_homologs(hlog_list, max_gap):
    org_dict = {}

    # Stage 1:  read the list of homologs in, and organize based on accession.  Each accession will have a list of homologs for a given gene block.
    # Prior to this, the program does not sort the output.
    # This section has been tested and validated to the best of my abilities.
    #print len(hlog_list)
    for item in hlog_list:
        accession = item.accession()
        #print accession
        if accession in list(org_dict.keys()):
            org_dict[accession].append(item)
        else:
            org_dict.update({accession:[item]})
    
    
    # Stage 2: Sort the list of homologs for each organism.  Determine gene blocks based on the max_gap criterion, 
    # and reject organisms without a gene block.  
    # This section has been tested, but not extensively.  I have added resolve_multiple_ORF_hits which is untested.
    for accession in sorted(org_dict.keys()):
        h_list = org_dict.pop(accession)
        h_list.sort(key=lambda x: x.start())
        
        # Here is where the code dealing explicitly with multiple hits to a single ORF goes:
        # currently, we only use best hit. Other resolution schemes can be envisioned.
        ORF_filetered_hlist = resolve_multiple_ORF_hits(h_list)
        org_dict.update({accession:ORF_filetered_hlist})
    
    # Stage 3: Organize the homologs into neighborhoods.  We remove any organisms that lack neighboring genes.
    # The return from this stage is a list of lists. Where each sub-list is a gene block, as defined by max_gap.
    # This version is not completely tested, but appears to be working when tested against known cases.
    
    neighborhood_dict = {}
    for accession in sorted(org_dict.keys()):
        hlist = org_dict.pop(accession)
        gene_block_list, neighborhood_found = group_homologs(hlist, max_gap)
        geneblock =''
        if accession =="NC_002516" or accession =="NC_004463":
            for blocks in hlist:
                if "cai" in blocks.blast_annotation():
                    
                    geneblock += blocks.blast_annotation()+ ',' + str(blocks.query_start())+ ',' + str(blocks.query_stop()) +'\t'
        if len(geneblock) != 0:
            print(accession)
            print(geneblock)
        
        if neighborhood_found:
            neighborhood_dict.update({accession:gene_block_list})
            org_dict.update({accession:hlist})

        # Enable the print organism bit to see what we filter out initially...    
        else: # do nothing, there are no neighborhoods that have been recovered
            #print "accession", accession, "is missing."
            #print "Organism ", hlist[0].organism(), "is missing."
            #print hlist
            pass
    
    # An explanation on what each of these returned dictionaries contains:
    # Each has organisms that only contain gene neighborhoods.
    # neighborhood_dict: accession keyed dict whose data is a list of neighborhods. organism1:[[h1, h2], [h3, h4], [h5]]
    # org_dict: accession keyed dict whose data is a list of ungrouped homologs. organism1:[h1, h2, h3, h4, h5]
    # org_dict differs from the inpupt (besides being a dict and not a list) by 
    return neighborhood_dict, org_dict

# I think this version is more readable than those i have made in the past. 
# It can take either a sorted, or unsorted list of homologs.   
def group_homologs(lst_homologs, max_gap):
    # bypassing pass by reference in python that leads to potentially undesirable behavior downstream
    list_homologs = [i for i in lst_homologs]
    
    # Step 1: Sort the input list by the start position of the ORF
    list_homologs.sort(key=lambda x: x.start())
    
    # Step 2: Group homologs into gene blocks as defined my max_gap, and report these groups.
    result, neighborhood_found = homolog_list_grouping_function(list_homologs, max_gap)
    return result, neighborhood_found

# This function will take a list of ordered homologs, which have had their redundant BLAST thits removed, and group them by a max_gap constraint.
# The return is a list of lists. Single genes and gene blocks will both be represented as groups.
def homolog_list_grouping_function(list_homologs, max_gap):
    result = []
    neighborhood = [list_homologs[0]]
    neighborhood_found = False
    
    for i in list_homologs[1:]:
        #look at current
        start = neighborhood[-1].start() #start = list_homologs[i].start()
        stop = neighborhood[-1].stop() #stop = list_homologs[i].stop()
        # look at next
        start_n = i.start() #start_n = list_homologs[i+1].start()
        stop_n = i.stop() #stop_n = list_homologs[i+1].stop()
        
        # We have found neighboring genes, as defined by max_gap
        if math.fabs(start - stop_n) < max_gap or math.fabs(stop - start_n) < max_gap:
            neighborhood_found = True
            neighborhood.append(i)
        # These genes do not neighbor eachother
        else: 
            result.append(neighborhood)
            neighborhood = [i]
    result.append(neighborhood)
    #print list_homologs[0].organism(), "result", result, "neighborhood_found ", neighborhood_found  
    return result, neighborhood_found 

# fast implementation of an order preserving make unique function
def make_unique(lst, function):#lambda a: a.start): 
    seen = {}
    result = []
    for item in lst:
        marker = function(item)
        if marker not in seen:
            seen.update({marker:1})
            result.append(item)
    return result

# This function will take the grouped lists, and determine which gene block genes do not occur as part of a neighborhood.
# These singletons genes will further be filtered for the best e-val as reported by BLAST. (Though other filtering schemes could be developed later)
def return_best_singleton_genes(grouped_lists):
    # the result will be a list of lists.
    result = []

    # Step 1: determine all recovered genes (for a given gene block) in this organism
    
    # return a list of evey homolog found by BLAST as a single list
    organism_genes = list(itertools.chain(*grouped_lists))
    # make list above unique based on the blast_annotation
    unique_in_organism_by_annotation = make_unique(organism_genes, lambda x: x.blast_annotation())
    # resulting in a list of annotations that are unique to the organism
    organism_annotation_list = [i.blast_annotation() for i in unique_in_organism_by_annotation]
    
    #print "organism_annotation_list", organism_annotation_list
    
    # Step 2: determine genes that are found in neighborhoods as defined by the max_gap criterion
    
    # find genes that are grouped
    neighborhoods_only = [i for i in grouped_lists if len(i) > 1]
    # make into a single list of homologs
    neighborhood_hlog_list = list(itertools.chain(*neighborhoods_only))
    # make list above unique based on the blast_annotation
    unique_neighborhood_by_annotation = make_unique(neighborhood_hlog_list, lambda x: x.blast_annotation())
    # resulting in a list of annotations that are unique to neighboring genes
    neighborhood_annotation_list = [i.blast_annotation() for i in unique_neighborhood_by_annotation]
    
    #print "neighborhood_annotation_list", neighborhood_annotation_list
    

    # Step 3: return a list of singletons genes, that are only found as singletons in the current genome
    single_annotation_list = list(set(organism_annotation_list) - set(neighborhood_annotation_list))
    singleton_genes = [i for i in grouped_lists if len(i) == 1]
    filtered_singlgeton_genes = [i for i in singleton_genes if i[0].blast_annotation() in single_annotation_list]
    
    # Step 4: Find the singleton gene that has the most significant e-value per BLAST annotation
    best_singleton_dict = {}
    for tmp in filtered_singlgeton_genes:
        gene = tmp[0]
        # a singleton gene with this annotation has already been found, use e-val to determine wich is the more significant hit
        if gene.blast_annotation() in list(best_singleton_dict.keys()):
            old_gene = best_singleton_dict.pop(gene.blast_annotation())
            if old_gene.e_val() <= gene.e_val(): # the existing homolog is a more significant hit
                best_singleton_dict.update({gene.blast_annotation(): old_gene})
            else: # the new homolog is a more significant hit
                best_singleton_dict.update({gene.blast_annotation(): gene})
        else: # we have not seen this annotation yet, store it in the dictionary
            best_singleton_dict.update({gene.blast_annotation(): gene})
            
    # Step 5: return a list of lists, [[s1], [s2], [s3]] for each entry in the best_singleton_dict.
    for i in list(best_singleton_dict.keys()):
        result.append([best_singleton_dict[i]])
    return result
            
    

# This version is not rigorously tested, but seems to work correctly.
# Return a list of lists of homologs = [[],[]], number of splits, and number of duplications. 
# unique_genes_in_organism, len_gene_block are integers grouped_list is a list of lists, and only contains groups 2 or more.
def optimize_neighborhoods(grouped_lists):

    # Step 1: make the grouped lists into a single list of homologs
    org_hlog_list = list(itertools.chain(*grouped_lists))
    
    neighborhoods_only = [i for i in grouped_lists if len(i) > 1]
    grouped_hlog_list = list(itertools.chain(*neighborhoods_only))
    
    # Step 2: determine the number of unique genes in both neighborhoods, and the organism.
    # To better explain: I need to know the number of unique genes the organism contains for the gene block.
    # I also need to know the number of unique genes found in neighborhoods. 
    number_unique_genes_in_organism = len(make_unique(org_hlog_list, lambda x: x.blast_annotation()))
    number_unique_in_neighborhoods = len(make_unique(grouped_hlog_list, lambda x: x.blast_annotation()))
    
    '''
    # Debugging.  This does check out.  kinda interesting stuff though here, there are some inline tandem repeats where gene name is different.
    # Everything looks like it works correctly though
    print org_hlog_list[0].organism(), "number_unique_genes_in_organism", number_unique_genes_in_organism, "number_unique_in_neighborhoods", number_unique_in_neighborhoods
    if number_unique_in_neighborhoods == 1:
        print "neighborhoods_only", neighborhoods_only
        for group in neighborhoods_only:
            for gene in group:
                print gene.blast_annotation(), gene.genbank_annotation(), gene.locus(), gene.start(), gene.stop()
        #print "grouped_lists",grouped_lists
    '''
    # Step 3: greedy algorithm to determine the best neighborhoods to report as a final result
    
    optimal = False
    num_in_list = 1 # this is the number of elements per list reurned
    best_duplicates = 0 
    splits = number_unique_genes_in_organism - number_unique_in_neighborhoods
    while not optimal:
        for group in itertools.combinations(grouped_lists, num_in_list):
            #all_homologs_in_grouping = [item for sublist in group for item in sublist]
            all_homologs_in_grouping = list(itertools.chain(*group))
            
            #print all_homologs_in_grouping
            #unique_in_set = len(MakeUnique(all_homologs_in_grouping, lambda a: a.predicted_gene))
            unique_in_set = len(make_unique(all_homologs_in_grouping, lambda x: x.blast_annotation()))
            #if unique_in_set == len_unique_grouped: # we have an optimal solution, perhaps not global optima
            if unique_in_set == number_unique_in_neighborhoods: # we have an optimal solution, perhaps not global optima
                duplicates =  int(math.fabs(len(all_homologs_in_grouping) - number_unique_in_neighborhoods))
                if not optimal:
                    optimal = True
                    best_grouping = list(group)
                    best_duplicates = duplicates
                    best_split = splits
                elif duplicates < best_duplicates:
                    best_grouping = list(group)
                    best_duplicates = duplicates
        splits+=1
        num_in_list+=1
    #print "splits " , splits, ": best_split ", best_split
    #print "Best grouping as found by the program\n", best_grouping
    
    # Step 4: determine the best (if necessary) singleton genes to complete
    if number_unique_genes_in_organism != number_unique_in_neighborhoods:
        # This step takes time, so only perform it when you have to
        best_singletons = return_best_singleton_genes(grouped_lists)
        #print "Difference", number_unique_genes_in_organism - number_unique_in_neighborhoods, len(best_singletons)
        #print "singletons", best_singletons , ' '.join([i.blast_annotation() for i in list(itertools.chain(*best_singletons))])
        best_grouping = best_grouping + best_singletons
    
    
    return best_grouping, best_split, best_duplicates#, len_unique_grouped

# does not seem to be in use
'''
def parallel_filter_operons(arg_tuple):
    fname, outfolder, max_gap, e_val = arg_tuple
'''    
   
def main():
    start = time.time()
    
    parsed_args = parser_code()
    
    infolder, outfolder, filter_file, num_proc, e_val, max_gap, quiet = check_options(parsed_args)
    
    if not quiet:
        print(infolder, outfolder, filter_file, num_proc, e_val, max_gap)
    
    file_list = return_file_list(infolder, filter_file)
    parallel_list_param = [(i, outfolder, max_gap, e_val) for i in file_list]

    for fname in file_list:
        #print fname
        hlog_list = filter_eval(fname, e_val)
        # here we return two dictionaries that are keyed by org, and contain at least one gene block. Orgs without gene blocks are omitted.
        neighborhood_dict, org_dict = return_valid_organism_homologs(hlog_list, max_gap)
        '''if "caiTABCDE" in fname:
            for accession in neighborhood_dict:
                if accession =="NC_002516" or accession =="NC_004463":
                    print accession
                    for blocks in neighborhood_dict[accession]:
                        geneblock =''
                        for homolog in blocks:
                            geneblock += homolog.blast_annotation()+ ',' + str(homolog.query_start())+ ',' + str(homolog.query_stop()) +'\t'
                        print geneblock '''
        # open a file handle to the output 
        head, tail = os.path.split(fname)
        outfile = outfolder + tail
        handle = open(outfile, 'w')
        
        # Save the result, in the result folder,  just like the good little program you are.
        for org in sorted(neighborhood_dict.keys()):
            best_grouping, best_split, best_duplicates = optimize_neighborhoods(neighborhood_dict[org])
            grouping_list = sorted(list(itertools.chain(*best_grouping)), key=lambda x: x.start())
            handle.write('\n'.join([i.to_file() for i in grouping_list])+'\n')
        handle.close()

    if not quiet:
             print(time.time() - start)

    # ./filter_operon_blast_results.py -f phylo_order.txt

if __name__ == '__main__':
    main()
