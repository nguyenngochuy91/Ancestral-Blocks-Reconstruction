#!/usr/bin/env python

from multiprocessing import Pool
import time
import os
import sys
import argparse

# Copyright(C) 2015 David Ream
# Released under GPL version 3 licence. http://www.gnu.org/licenses/lgpl.html
# Do not remove this comment

# This exists to  make the main function easier to read. It contains code to run the argument parser, and does nothing else.
def parser_code():

    parser = argparse.ArgumentParser(description='Run BLAST on a set of searchable database using a query that is specified by the user. Results will be stored by the accession number of the database. Currently this script will only accept protein queries, but I will update to automatically run on all types of genes, as most of the information needed for this behavior exists.')

    parser.add_argument("-d", "--database_folder", dest="database_folder", metavar="DIRECTORY", default='./db/',
                help="Folder containing all BLAST searchable databases to be used by the program.")
                 
    parser.add_argument("-o", "--outfolder", dest="outfolder", metavar="DIRECTORY", default='./blast_result/',
                help="Folder where the BLAST results will be stored.")
    
    parser.add_argument("-f", "--filter", dest="filter", metavar="FILE", default='NONE',
                help="File restrictiong which accession numbers this script will process. If no file is provided, filtering is not performed.")
                          
    parser.add_argument("-n", "--num_proc", dest="num_proc", metavar="INT", default = os.sysconf("SC_NPROCESSORS_CONF"), type=int,
                help="Number of processors that you want this script to run on. The default is every CPU that the system has.")
    
    # Fix this option, ultimately it should be a folder, that reads in a gene block file with the name of gene block(s) with a list of gene names and types.
    # The program will then take two files, protein and nucleic acid queries and run them. (This would require that there are two seperate, and complementary
    # blast databases within this folder. My desire is that it would be two subfolders, 'protein/' and 'rna/' which house these sets of data.
    parser.add_argument("-u", "--query", dest="query", default='gene_block_query.fa', metavar="FILE",
                help="A file that contains the BLAST query for every gene of interest in the dataset.")
                
    parser.add_argument("-e", "--eval", dest="eval", default='1e-10', metavar="FLOAT", type=float,
                help="eval for the BLAST search.")
    
    parser.add_argument("-q", "--quiet", dest="quiet", action="store_true", default=False,
                help="Suppresses most program text outputs.")
                
    return parser.parse_args()


def check_options(parsed_args):
    # section of code that checks the database entry    
    if os.path.isdir(parsed_args.database_folder):
        database_folder = parsed_args.database_folder
    else:
        print("The database directory %s does not exist." % parsed_args.database_folder)
        sys.exit()

    # if the directory that the user specifies does not exist, then the program makes it for them. 
    if not os.path.isdir(parsed_args.outfolder):
        os.makedirs(parsed_args.outfolder)
    outfolder = parsed_args.outfolder
    
    # Check the filter file
    if parsed_args.filter == 'NONE' or os.path.exists(parsed_args.filter):
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
    
    # Check the query file
    if os.path.exists(parsed_args.query):
        query_file = parsed_args.query
    else:
        print("The query file %s does not exist." % parsed_args.query)
        sys.exit()
    
    e_val = parsed_args.eval
        
    quiet = parsed_args.quiet
    
    return database_folder, outfolder, filter_file, num_proc, query_file, e_val, quiet
 
        
#this function will return all of the files that are in a directory. os.walk is recursive traversal.
def returnRecursiveDirFiles(root_dir):
    result = []
    for path, dir_name, flist in os.walk(root_dir):
        for f in flist:
            fname = os.path.join(path, f)
            if os.path.isfile(fname):
                result.append(fname)
    return result


# This code right now only deals with protein, but I will add functionality later for nucleotides. 
# Just moving the project along here, but this is a critical flaw moving forward.
def do_parallel_blast(arg_tuple):
    db, query_file, blast_result_folder, num_processors, eval_threshold = arg_tuple
    out_file = "%s%s.txt" % (blast_result_folder, db.split('/')[-1].split('.')[0])
    #print "db", db
    #print "got here"
    #cmd = "blastall -p tblastn -a %i -i %s -d %s -e %s -o %s -m 9" % (num_processors, query_file, db, eval_threshold, out_file)
    #cmd = "blastall -p tblastn -a %i -i %s -d %s -e %s -o %s -m 9" % (1, query_file, db, eval_threshold, out_file)
    #cmd = "blastall -p blastp -a %i -i %s -d %s -e %s -o %s -m 9" % (1, query_file, db, eval_threshold, out_file)
    cmd ="blastp -num_threads %i -query %s -db %s -evalue %s -out %s -outfmt 6 -seg yes" % (1, query_file, db, eval_threshold, out_file)
    #print cmd
    os.system( cmd )


#def parallel_blast(infile, query, folder, num_proc, e_val = '1e-10'):
def parallel_blast(database_folder, outfolder, filter_file, num_proc, query_file, e_val):
    # you kinda have to trust me here, but having blast run on as many threads per CPU as you have total processors is fastest
    # I have no idea why this is... ugh.
    
    unfiltered_db_list = [i for i in returnRecursiveDirFiles(database_folder) if i.split('/')[-1].split('.')[-1] == 'ffc']
    # print 'filter_file',filter_file
    if filter_file == '' or filter_file == 'NONE' or len(filter_file)<5:
        db_list = unfiltered_db_list
    else:
        filter_list = [i.strip() for i in open(filter_file).readlines()]
        db_list = [i for i in unfiltered_db_list if i.split('/')[-1].split('.')[0] in filter_list]
    
    #print len(unfiltered_db_list), len(db_list)
    
    #blast_arg_list = [(i, query_file, outfolder, num_proc, e_val) for i in db_list]
    blast_arg_list = [(i, query_file, outfolder, 1, e_val) for i in db_list]
    pool = Pool(processes = num_proc)
    pool.map(do_parallel_blast, blast_arg_list)


def main():
    
    start = time.time()
    
    parsed_args = parser_code()
    
    database_folder, outfolder, filter_file, num_proc, query_file, e_val, quiet = check_options(parsed_args)
    
    if not quiet:
        print(database_folder, outfolder, filter_file, num_proc, query_file, e_val, quiet)
    
    #parallel_blast(infile, query, folder, num_proc, e_val)
    parallel_blast(database_folder, outfolder, filter_file, num_proc, query_file, e_val)

    if not quiet:
        print(time.time() - start)
    
    # ./blast_script.py -d ./db/ -o ./blast_result/ -q ./gene_block_query.fa 
    
if __name__ == '__main__':
    main()

