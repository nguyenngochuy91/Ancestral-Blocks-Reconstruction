#!/usr/bin/env python
''' Author  : Huy Nguyen
    Program : ROAGUE pipeline, which include the finding of orthoblocks in set of given taxa
                and reconstruct the ancestral gene blocks.
                This is basically a program that calls several other program
    Start   : 09/08/2016
    End     : 09/1/2016
'''

import argparse
import os


## Traverses the genome information directory
def traverseAll(path):
    res=[]
    for root,dirs,files in os.walk(path):
        for f in files:
            res.append(root+'/'+f)
    return res    
    
### parsing the argument from user input

def parser_code():

    parser = argparse.ArgumentParser()
                     
    parser.add_argument("--genomes_directory","-g", help="The directory that store all the genomes file (E_Coli/genomes)")  
    parser.add_argument("--gene_blocks","-b", help="The gene_block_names_and_genes.txt file, this file stores the operon name and its set of genes") 
    parser.add_argument("--reference","-r", help="The ncbi accession number for the reference genome (NC_000913 for E_Coli and NC_000964 for B_Sub)")  
    parser.add_argument("--filter","-f", help="The filter file for creating the tree (E_Coli/phylo_order.txt for E_Coli or B_Sub/phylo_order.txt for B-Sub)")  
    parser.add_argument("--method","-m", help="The method to reconstruc ancestral gene block, we support either global or local")       
    parser.add_argument("--output","-o", help="Output directory to store the result",default = "result")                  
    return parser.parse_args()

if __name__ == '__main__':
    args                       = parser_code()
    genomes_directory          = args.genomes_directory
    reference                  = args.reference
    filter_file                = args.filter
    method                     = args.method
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
        
    ### format a database for faster blasting. output in db
    db = parent_dir+ 'db'
    cmd1 ='./format_db.py -i {} -o {}'.format(genomes_directory,db)
    os.system(cmd1)
    print ('cmd1:',cmd1)
    
    ### Given the gene_block_names_and_genes.txt, create a gene_block_query.fa using the reference gene bank file. output in file gene_block_query.fa
    gene_block_names_and_genes = dirs[0]+"/"+'gene_block_names_and_genes.txt'
    gene_block_query           = parent_dir +'gene_block_query.fa'
    cmd2 ='./make_operon_query.py -i {} -b {} -r {} -o {}'.format(genomes_directory,gene_block_names_and_genes,reference,gene_block_query)
    os.system(cmd2)
    print ('cmd2:',cmd2)
    
    ### blasting using db vs the gene_block_query.fa above. output in blast_result
    blast_result = parent_dir+'blast_result/'
    cmd3 ='./blast_script.py -u {} -d {} -o {}'.format(gene_block_query,db,blast_result)
    os.system(cmd3)
    print ('cmd3:',cmd3)
    
    ### parsing the blast result directory into files that group by operon names, output in blast_parse
    blast_parse = parent_dir+'blast_parse/'
    cmd4 ='./blast_parse.py -b {} -i {} -o {}'.format(gene_block_names_and_genes,blast_result,blast_parse)
    os.system(cmd4)
    print ('cmd4:',cmd4)    
    
    ### filtering the gene blocks so that we have the most optimal gene blocks given the blast parse directory, ouput to optimized_gene_block
    optimized_gene_block = parent_dir+'optimized_gene_block/'
    cmd5 ='./filter_operon_blast_results.py -i {} -o {}'.format(blast_parse,optimized_gene_block)
    os.system(cmd5)
    print ('cmd5:',cmd5)

    ### create newick tree file. output in tree directory
    tree_dir = parent_dir+'tree/'
    cmd6 ='./create_newick_tree.py -G {} -f {} -r {} -o {}'.format(genomes_directory,filter_file,reference,tree_dir)
    os.system(cmd6)
    print ('cmd6:',cmd6)    
    
    ### from the filter_operon_blast_results, create a result directory. output in result directory
    accession = tree_dir + 'accession_to_common.csv'
    result    = parent_dir + 'result/'
    cmd7      = './get_result.py -g {} -a {} -o {} -i {}'.format(gene_block_names_and_genes,accession,result,optimized_gene_block)
    os.system(cmd7)
    print ('cmd7:',cmd7)
    
    ### generate group.txt file to color taxa based on class
    group    = parent_dir+'group.txt'
    cmd8     ='./group.py -i {} -o {} -a {}'.format(genomes_directory,group,filter_file) 
    os.system(cmd8)
    print ('cmd8:',cmd8)    
    ###########################################################################
    ## reconstruction of ancestral gene blocks
    
    ### convert the file in directory result and store it in new_result
    new_result = parent_dir+'new_result/'
    cmd9       = './convert.py -i {} -o {}'.format(result,new_result)
    os.system(cmd9)
    print ('cmd9:',cmd9)    
    
    ### reconstructing the gene block using global method
    reconstruction = parent_dir+'reconstruction_'+method
    tree = tree_dir+'out_tree.nwk'
    cmd10 ='./reconstruction.py -i {} -t {} -o {} -m {}'.format(new_result,tree,reconstruction,method)
    os.system(cmd10)
    print ('cmd10:',cmd10)
    
    ### create a visualization of the reconstruction and store it in directory visualization
    visualization = parent_dir+'/visualization/'
    try:
        os.mkdir(visualization)
    except:
        print ("visualization directory is already created" )
    res = traverseAll(reconstruction)
    # go through all the reconstruction ancestral file and create a pdf visualization of them, using the group text to color group
    for file in res:
        operon = file.split('/')[-1]
        if 'mapping' in operon:
            continue
        outfile = visualization+operon
        cmd11 = './show_tree.py -i {} -g {} -o {}'.format(file,group,outfile)
        os.system(cmd11)
        
        