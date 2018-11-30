#!/usr/bin/env python

# Copyright(C) 2015 David Ream
# Released under GPL version 3 licence. http://www.gnu.org/licenses/lgpl.html
# Do not remove this comment

# This script generates a phylogenetic tree from a list of genbank accession numbers and a marker gene.

# I do not know how many of these i need, i will edit later
import os
import sys
import time
import argparse
import shutil
from Bio import SeqIO,SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio import Application
from Bio.Application import _Option
from Bio.Align.Applications import MuscleCommandline
from Bio import Phylo, AlignIO
import subprocess
import traceback
from ete3 import Tree

# Globals
# The three output file names will be stored in the output directory that is supplied by the user
#newick_tree = 'out_tree.nwk'
#outfile_for_asma = 'accession_to_common.txt'
#accession_to_common = 'accession_to_common.csv'
#phylo_order_new = 'phylo_order_new.txt'

# This exists to  make the main function easier to read. It contains code to run the argument parser, and does nothing else.
def parser_code():

    parser = argparse.ArgumentParser(description='The purpose of this script is to build a newick format phylogenetic tree from a list of genomes and a marker gene.')
    
    parser.add_argument("-G", "--genbank_directory", dest="genbank_directory", metavar="DIRECTORY", default='./genomes/',
                help="Folder containing all genbank files for use by the program.")
    
    parser.add_argument("-o", "--outfolder", dest="outfolder", metavar="DIRECTORY", default='./tree/',
                help="Directory where the results of this program will be stored.")
                
    parser.add_argument("-f", "--filter", dest="filter", metavar="FILE", default='None',
                help="File restrictiong which accession numbers this script will process. If no file is provided, filtering is not performed.")
                
    parser.add_argument("-m", "--marker_gene", dest="marker_gene", metavar="STRING", default='rpob',
                help="This is a single marker gene that will be used to construct phylogenetic trees.")
                
    parser.add_argument("-t", "--tree", dest="tree_file", metavar="FILE", default='None',
                help="Newick format tree file which will be used to bypass tree creation.")
    parser.add_argument("-r","--ref",default='NC_000913',help = 'The reference genome')
    parser.add_argument("-q", "--quiet", dest="quiet", action="store_true", default=False,
                help="Suppresses most program text outputs.")
                                        
    return parser.parse_args()


def check_options(parsed_args):
    
    # check the genbank folder
    if os.path.isdir(parsed_args.genbank_directory):
        genbank_directory = parsed_args.genbank_directory
    else:
        print("The folder %s does not exist." % parsed_args.genbank_directory)
        sys.exit()
    
    if parsed_args.filter == 'None' or os.path.exists(parsed_args.filter):
        filter_file = parsed_args.filter
    else:   
        print("The file %s does not exist." % parsed_args.filter)
        sys.exit()
    
    outfolder = parsed_args.outfolder
    
    
    # check the output directory
    if os.path.isdir(parsed_args.outfolder):
        outfolder = parsed_args.outfolder
    else:
        os.mkdir(parsed_args.outfolder)
    outfolder = parsed_args.outfolder
    
    marker_gene = parsed_args.marker_gene
    
    # check some information about the tree file
    if parsed_args.tree_file == 'None' or os.path.exists(parsed_args.tree_file):
        tree_file = parsed_args.tree_file
    else:
        print("The file %s does not exist." % parsed_args.tree_file)
        sys.exit()
    
    quiet = parsed_args.quiet
    
    #return genbank_directory, outfolder, filter_file, marker_gene
    return genbank_directory, outfolder, filter_file, marker_gene, tree_file, quiet
    
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
    if filter_file == '' or filter_file == 'None':
        return return_recursive_dir_files(infolder)   
    else:
        filter_list = [i.strip() for i in open(filter_file)]

        #print "filter_list", filter_list
        #print "return",[i for i in return_recursive_dir_files(infolder) if os.path.basename(i).split('.')[0] in filter_list]
        return [i for i in return_recursive_dir_files(infolder) if os.path.basename(i).split('.')[0] in filter_list]
        

    
def make_common_to_accession_dict(infolder, filter_file):
    org_paths = return_file_list(infolder, filter_file)
    common_to_accession_dict = {}
    for org in org_paths:
        seq_record = next(SeqIO.parse(open(org), "genbank"))
        accession = seq_record.annotations['accessions'][0]
        organism_tmp = seq_record.annotations['organism'].replace(' ', '_')
        organism_tmp = organism_tmp.replace('[','')
        organism_tmp = organism_tmp.replace(']','')
        organism = '_'.join(organism_tmp.split('_')[:2])
        common_to_accession_dict.update({organism:accession})
    return common_to_accession_dict
        
# this function will take a marker gene and a list of genbank files
# and construct a fasta file for their sequences.  For the time being at least
# this function will only do protein marker genes. Later i will expand this function
# to take RNA genes like 16s as well.
# TODO : expand the functionality of this function to make use of RNA genes and give the user the ability to list more than one marker that they are interested in
def make_target_fasta(marker, infolder, filter_file, marker_fasta):
    org_paths = return_file_list(infolder, filter_file)
    #print org_paths
    #print "infolder", infolder
    #print "filter_file", filter_file
    print("org_paths",org_paths)
    result = [] # this is not great
    common_to_accession_dict = {}
    orgs_with_marker = []
    marker = marker.lower()
    #print len(org_paths)
    for org in org_paths:
        genes_found = []
        seq_record = next(SeqIO.parse(open(org), "genbank"))
        accession = seq_record.annotations['accessions'][0]
        organism_tmp = seq_record.annotations['organism'].replace(' ', '_')
        organism_tmp = organism_tmp.replace('[','')
        organism_tmp = organism_tmp.replace(']','')
        # Put code here to determine the format of the organisms' english name. currently i am using genus species, but strain can also be used
        organism = '_'.join(organism_tmp.split('_')[:2])#+"_"+accession
        #if(organism == "Natranaerobius_thermophilus") :
                #print accession
        # Here we store the {organism:accession} information so that we can build a new list that is needed for the visualization pipeline
        common_to_accession_dict.update({organism:accession})
        found = False
        for fnum, feature in enumerate(seq_record.features):
            #if((organism == "Natranaerobius_thermophilus") and feature.type == 'CDS') :
                #print feature           
            if feature.type == 'CDS':
                #start = feature.location._start.position
                start = feature.location.start
                stop = feature.location.end
                try: 
                    gene = feature.qualifiers['gene'][0]
                    gene = gene.lower()
                except:
                    gene = 'unknown'
                if gene == marker:
                    genes_found.append(gene)
                    seq = feature.qualifiers['translation'] # this returns the protein product, not suitable for RNA products like 16s
                    #result.append(">%s|%s|%s" % (accession, organism, gene))
                    #print organism+str(seq)
                    result.append(">%s,%s" % (organism, accession))
                    result.append(''.join(seq))
                    orgs_with_marker.append(accession)
                    found = True
                    break
        #found = True
        #print marker
        #This method is to get the rpoB gene protein sequence in those genomes which doesn't have it in their CDS record, but have its coordinates in their misc_feature.
        if(not found):
            #file = open(organism+".txt","w")
            for fnum, feature in enumerate(seq_record.features):
            #if((organism == "Natranaerobius_thermophilus") and feature.type == 'CDS') :
                #print feature           
                if feature.type == 'misc_feature':
                #start = feature.location._start.position
                    start = feature.location.start
                    stop = feature.location.end
                    try: 
                        note = feature.qualifiers['note'][0]
                        noteSplit = note.split(";")
                        region = 'unknown'
                        for i in noteSplit:
                            if('Region' in i):
                                region = i.split(":")[1].strip()
                                break
                    except:
                        region = 'unknown'
                    #file.write(region+"\n")
                    if region.lower() == marker.lower():
                        #print str(start)+" "+(str(stop))
                        genes_found.append(gene)
                        seq = seq_record.seq
                        seq = str(seq[start:stop].translate()) # this returns the protein product, not suitable for RNA products like 16s
                        #result.append(">%s|%s|%s" % (accession, organism, gene))
                        #print organism+seq
                        result.append(">%s,%s" % (organism, accession)) #make sure they are different
                        result.append(''.join(seq))
                        orgs_with_marker.append(accession)
                        found = True
                        break
        #print found   #file.close()
    #outfile = tmp_directory + "distmat_marker.fa"
    #print "outfile", outfile
    #handle = open(outfile, 'w')
    handle = open(marker_fasta, 'w')
    handle.write('\n'.join(result))
    handle.close()
    
    
    #print "orgs_with_marker", len(orgs_with_marker)
    #create_distmat(outfile)
    
    #print "common_to_accession_dict", common_to_accession_dict.keys()
    #return "./msa_tmp/distmat_marker.ph", common_to_accession_dict
    #print len(common_to_accession_dict)
    return common_to_accession_dict

##Not Used in the main Code
def return_tree_order_list(newick_tree_file, common_to_accession_dict):
    result = []
    tree = Phylo.read("out_tree.nwk", "newick")
    for clade in tree.find_clades():
        if clade.name != None:
            #accession, common_tmp, marker = clade.name.split('|')
            #common = '_'.join(common_tmp.split('_')[:2])
            common = clade.name
            accession = common_to_accession_dict[common]
            result.append(','.join([accession, common]))
            clade = common
    #print '\n'.join(result)
    
    '''
    outfile_for_asma = './asma_outlist.txt'
    handle = open(outfile_for_asma, 'w')
    handle.write('\n'.join(result))
    handle.close()
    '''
    
    outfile_for_asma = './accession_to_common.txt'
    handle = open(outfile_for_asma, 'w')
    handle.write('\n'.join(result))
    handle.close()
    
    
    phylo_order_new = './phylo_order_new.txt'
    handle = handle = open(phylo_order_new, 'w')
    handle.write('\n'.join([i.split(',')[0] for i in result]))
    handle.close()
    
    
    #Phylo.write(tree, "./test_tree.nwk", "newick")


def make_newick_tree(marker_fasta, tree_outfile,reference):
    ## using muscle
    # make the alignment file
    temp_align = '.'.join(marker_fasta.split('.')[:-1]) + '.aln' 
    cm1 ="muscle -in "+marker_fasta+ " -out "+temp_align
    os.system(cm1)
    #make the tree using clustal
    cm2 ="clustalw -infile="+temp_align+" -tree=1"
    # have to wait for few second for the aln file actually comes out lol
    os.system(cm2)
    temp_tree = '.'.join(marker_fasta.split('.')[:-1]) + '.ph' # that's what this file gets named by default, and i'm sick of looking for the cmd line arg to fix.
    print(temp_tree)
    print("modifying")
    #modify for negative branch
    modify_tree = '.'.join(marker_fasta.split('.')[:-1]) + '.new'
    cm3 = "sed -e 's,:-[0-9\.]\+,:0.0,g' "+temp_tree+" > "+modify_tree   
    os.system(cm3)
    
    if reference == 'NC_000913':
        t= Tree(modify_tree)
        ancestor = t.get_common_ancestor("Campylobacter_jejuni_NC_002163","Nitrosomonas_europaea_NC_004757")
        t.set_outgroup(ancestor)
        t.write(format = 1, outfile = modify_tree)
    # dealing with negative branch length
    #print "marker_fasta",marker_fasta
    #print "temp_tree", temp_tree
    # move the created tree file to the location i say its going
    shutil.copy(modify_tree, tree_outfile)
    

def return_tree_order_list_2(newick_tree_file, common_to_accession_dict, accession_to_common_outfile, phylo_order_new_outfile):
    result = []
    tree = Phylo.read(newick_tree_file, "newick")
    for clade in tree.find_clades():
        if clade.name != None:
            #accession, common_tmp, marker = clade.name.split('|')
            #common = '_'.join(common_tmp.split('_')[:2])
            common_tmp = clade.name
            common_tmp = common_tmp.split('_')
            common = common_tmp[0]+'_'+common_tmp[1]
            accession = common_to_accession_dict[common]
            result.append(','.join([accession, common])) #make sure they are different
            clade = common
    #print '\n'.join(result)
    
    '''
    outfile_for_asma = './asma_outlist.txt'
    handle = open(outfile_for_asma, 'w')
    handle.write('\n'.join(result))
    handle.close()
    '''
    
    #outfile_for_asma = './accession_to_common.txt'
    #handle = open(outfile_for_asma, 'w')
    #print "accession_to_common_outfile", accession_to_common_outfile
    handle = open(accession_to_common_outfile, 'w')
    handle.write('\n'.join(result))
    handle.close()
    
    
    #phylo_order_new = './phylo_order_new.txt'
    #handle = handle = open(phylo_order_new, 'w')
    #print "phylo_order_new_outfile", phylo_order_new_outfile
    handle = handle = open(phylo_order_new_outfile, 'w')
    handle.write('\n'.join([i.split(',')[0] for i in result]))
    handle.close()
    
    
    #Phylo.write(tree, "./test_tree.nwk", "newick")


def main():

    #newick_tree = 'out_tree.nwk'
    #accession_to_common = 'accession_to_common.txt'
    ## outfile_for_asma = 'accession_to_common.txt'
    #phylo_order_new = 'phylo_order_new.txt'
    
    start = time.time()    

    parsed_args = parser_code()
    reference   = parsed_args.ref
    #genbank_directory, outfolder, filter_file, marker_gene = check_options(parsed_args)
    genbank_directory, outfolder, filter_file, marker_gene, tree_file, quiet = check_options(parsed_args)
    
    # This folder and its contents will be removed after this program concludes its run. i just need a place to store a bunch of intermediate files.
    tmp_directory = outfolder + "msa_tmp/"
    #print "tmp_directory from create newick", tmp_directory
    
    try:
        os.mkdir(tmp_directory)
    except:
        pass
        #traceback.print_exc()
        #print "directory Not made"
    #marker_gene = "rpob"
    
    if not quiet:
        print(genbank_directory, outfolder, filter_file, marker_gene)
    
    # set the paths for the three output files
    newick_tree_outfile = outfolder + 'out_tree.nwk'
    accession_to_common_outfile = outfolder + 'accession_to_common.csv' # formally: outfile_for_asma
    phylo_order_new_outfile = outfolder + 'phylo_order_new.txt'
    
    # Execute code fork that builds a tree from scratch
    if tree_file == 'None':
        #print "got here"
        # set the distmat file that is created in "make_target_fasta()" for i have no idea why
        marker_fasta = tmp_directory + "distmat_marker.fa"
        #print "got here 1"
        common_to_accession_dict = make_target_fasta(marker_gene, genbank_directory, filter_file, marker_fasta)
        #print "got here 2"
        make_newick_tree(marker_fasta, newick_tree_outfile,reference)
        # cm = "sed -e 's,:-[0-9\.]\+,:0.0,g' "+newick_tree_outfile+" > "+newick_tree_outfile
        # os.system(cm)t.write(format=1, outfile="new_tree.nw")
        #print "got here 3"
        return_tree_order_list_2(newick_tree_outfile, common_to_accession_dict, accession_to_common_outfile, phylo_order_new_outfile)
    
    
        #tree_file, common_to_accession_dict = make_target_fasta(marker_gene, genbank_directory, filter_file, distmat_outfile)
    
    
    
        #shutil.copy(tree_file, outfile)
    else:
        #print "Got here 2"
        common_to_accession_dict = make_common_to_accession_dict(genbank_directory, filter_file)
        return_tree_order_list_2(tree_file, common_to_accession_dict, accession_to_common_outfile, phylo_order_new_outfile)
    #return_tree_order_list(outfile, common_to_accession_dict)
    
    # get rid of the temp files
    # shutil.rmtree(tmp_directory)
    
    if not quiet:
        print(time.time() - start)
# ./create_newick_tree.py -i ./gram_positive_test_organisms/  -f gram_positive_phylo_order.txt    
if __name__ == '__main__':
    main()    

