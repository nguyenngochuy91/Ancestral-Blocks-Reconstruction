#!/usr/bin/env python
''' Author  : Huy Nguyen
    Program : Create operon tree
    Start   : 01/01/2018
    End     : 
'''
from Bio import SeqIO
import argparse
import os
import homolog4
import shutil
from ete3 import Tree
## Traverses the genome information directory
def traverseAll(path):
    res=[]
    for root,dirs,files in os.walk(path):
        for f in files:
            res.append(root+"/"+f)
    return res    
    
### parsing the argument from user input

def parser_code():

    parser = argparse.ArgumentParser()
                     
    parser.add_argument("--genomes_directory","-g", help="The directory that store all the genomes file (E_Coli/genomes)")  
    parser.add_argument("--gene_blocks","-b", help="The gene_block_names_and_genes.txt file, this file stores the operon name and its set of genes") 
    parser.add_argument("--reference","-r", help="The ncbi accession number for the reference genome (NC_000913 for E_Coli and NC_000964 for B_Sub)")  
    parser.add_argument("--filter","-f", help="The filter file for creating the tree (E_Coli/phylo_order.txt for E_Coli or B_Sub/phylo_order.txt for B-Sub)")         
    parser.add_argument("--output","-o", help="Output directory to store the result",default = "result")   
    
    return parser.parse_args()

### Given the optimized gene block data, and the genes in each species, provide a dictionary
### that has key as an operon, value is a dictionary where key is the species, and key is the gene sequence, start,stop and strand
def generate_operon(optimized_gene_block,db):
    db_dict = {}
    operon_dict ={}
    db_list = [i for i in traverseAll(db) if i.split('/')[-1].split('.')[-1] == 'ffc']
    for organism in db_list:
        accession_num = organism.split('/')[-1].split('.')[0]
        for record in SeqIO.parse(organism,"fasta"):
            info = record.id.split('|')
            name      = '_'.join(info[1].split('_')[:2])
            together  = name+"_"+accession_num
            start     = int(info[4])
            stop      = int(info[5])
            strand    = int(info[6])
            if together not in db_dict:
                db_dict[together] = {}
            db_dict[together][(start,stop,strand)] = record.seq.tostring()
    operon_list = traverseAll(optimized_gene_block)
    for operon in operon_list:
        operon_name = operon.split("/")[-1].split(".")[0]
        operon_dict[operon_name]= {}
        handle = open(operon,"r")
        for line in handle.readlines():
            h = homolog4.Homolog.from_blast(line)
            name   = '_'.join(h.organism().split('_')[:2])
            accession_num = h.accession().split(".")[0]
            together  = name+"_"+accession_num
            start  = h.start()
            stop   = h.stop()
            strand = h.strand()
            seq    = db_dict[together][(start,stop,strand)]
            if together not in operon_dict[operon_name]:
                operon_dict[operon_name][together] = [(seq,start,stop,strand)]
            else:
                operon_dict[operon_name][together].append((seq,start,stop,strand))
            
    return operon_dict

### given the gene position, determine whether there is a gene block, and then concatenate
### and assign split as a string of 500 letter "N"
def concatenate(potential):
    d = {}
    for gene in potential:
        seq    = gene[0]
        start  = gene[1]
        stop   = gene[2]
        strand = gene[3]
        if strand in d:
            d[strand].append((start,stop,seq))
        else:
            d[strand] = [(start,stop,seq)]
    has_block = False # flag to check whether there is a block
    wholeString = ""
    for strand in d:
        d[strand].sort()
    for strand in d:

        substring = d[strand][0][2] # always store the first gene sequence
        for i in range(len(d[strand])-1):
            gene1 = d[strand][i]
            gene2 = d[strand][i+1]
            if abs(gene2[0]-gene1[1])>500:
                substring+="N"*500
            else:
                has_block = True
            substring+=gene2[2]
        wholeString+= substring

    if has_block:
        return wholeString[:-1]+"\n"
    else:
        return ""

### given the operon dict above, generate fasta file to do alignment
def generate_fasta(operon_dict,fasta):
    for operon in operon_dict:
        outfile = open(fasta+operon,'w')
        for species in operon_dict[operon]:
            potential = operon_dict[operon][species]
            combine   = concatenate(potential)
            if combine:
                outfile.write(">"+species+"\n")
                outfile.write(combine)
        outfile.close()
    return None

### given the fasta, generate a multiple alignment and tree 
def generate_tree(fasta,tree):
    files = traverseAll(fasta)
    ## using muscle
    # make the alignment file
    for file in files:
        name = file.split('/')[-1]
        tree_name = tree+name+"/"
        try:
            os.mkdir(tree_name)
        except:
            print (tree_name+" directory already created")
        temp_align = tree_name +name+ '.aln' 
        cm1 ="muscle -in "+file+ " -out "+temp_align
        os.system(cm1)
        #make the tree using clustal
        cm2 ="clustalw -infile="+temp_align+" -tree=1"
        # have to wait for few second for the aln file actually comes out lol
        os.system(cm2)
        temp_tree = tree_name + name+ '.ph' # that's what this file gets named by default, and i'm sick of looking for the cmd line arg to fix.
        print(temp_tree)
        print("modifying")
        #modify for negative branch
        modify_tree = tree_name + name+ '.new'
        cm3 = "sed -e 's,:-[0-9\.]\+,:0.0,g' "+temp_tree+" > "+modify_tree   
        os.system(cm3)
        os.remove(temp_tree)
        # dealing with negative branch length
        #print "marker_fasta",marker_fasta
        #print "temp_tree", temp_tree
        # move the created tree file to the location i say its going

if __name__ == '__main__':
    args                       = parser_code()
    reference                  = args.reference
    genomes_directory          = args.genomes_directory
    reference                  = args.reference
    filter_file                = args.filter
    gene_block_names_and_genes                = args.gene_blocks
    outdir                     = args.output
    # check if we are going to output results in the current directory
    dirs = genomes_directory.split('/')
    outdir+='/'
    try:
        os.mkdir(outdir+'/')
    except:
        print ("output directory has already been created")
    if len(dirs)>=3: # means that we have to go to subdirectory
        parent_dir = outdir+dirs[0]+"/"
    else:
        parent_dir = outdir
    ##########################################################################
    # finding gene blocks           
    db = parent_dir+ 'db'
    gene_block_names_and_genes = dirs[0]+"/"+'gene_block_names_and_genes.txt'
    gene_block_query           = parent_dir +'gene_block_query.fa'
    blast_result = parent_dir+'blast_result/'
    blast_parse = parent_dir+'blast_parse/'
    optimized_gene_block = parent_dir+'optimized_gene_block/'
    tree = parent_dir +"operon_tree/"    
    try:
        os.mkdir(tree)
    except:
        print ("directory tree has already been created")
#    ### format a database for faster blasting. output in db
#    cmd1 ='./format_db.py -i {} -o {}'.format(genomes_directory,db)
#    os.system(cmd1)
#    print ('cmd1:',cmd1)
#    
#    ### Given the gene_block_names_and_genes.txt, create a gene_block_query.fa using the reference gene bank file. output in file gene_block_query.fa
#    cmd2 ='./make_operon_query.py -i {} -b {} -r {} -o {}'.format(genomes_directory,gene_block_names_and_genes,reference,gene_block_query)
#    os.system(cmd2)
#    print ('cmd2:',cmd2)
#    
#    ### blasting using db vs the gene_block_query.fa above. output in blast_result
#    cmd3 ='./blast_script.py -u {} -d {} -o {}'.format(gene_block_query,db,blast_result)
#    os.system(cmd3)
#    print ('cmd3:',cmd3)
#    
#    ### parsing the blast result directory into files that group by operon names, output in blast_parse
#    cmd4 ='./blast_parse.py -b {} -i {} -o {}'.format(gene_block_names_and_genes,blast_result,blast_parse)
#    os.system(cmd4)
#    print ('cmd4:',cmd4)    
#    
#    ### filtering the gene blocks so that we have the most optimal gene blocks given the blast parse directory, ouput to optimized_gene_block
#    cmd5 ='./filter_operon_blast_results.py -i {} -o {}'.format(blast_parse,optimized_gene_block)
#    os.system(cmd5)
#    print ('cmd5:',cmd5)
    
    ### from the filter_operon_blast_results, create a tree directory for the alignment and tree

    operon_dict = generate_operon(optimized_gene_block,db)
    fasta = parent_dir+"fasta/"
    try:
        os.mkdir(fasta)
    except:
        print ("directory fasta has already been created")
    generate_fasta(operon_dict,fasta)
    ### from the fasta file, generate a alignment and a tree
    generate_tree(fasta,tree)
    
    