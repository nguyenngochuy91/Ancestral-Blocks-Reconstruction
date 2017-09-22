# ROAGUE: **R**econstruction **o**f **A**ncestral **G**ene Blocks **U**sing **E**vents
## Purpose

ROAGUE is a tool to reconstruct ancestors of gene blocks in prokaryotic genomes. Gene blocks are genes co-located on the chromosome. In many cases, gene blocks are
conserved between bacterial species, sometimes as operons, when genes are co-transcribed. The conservation is rarely absolute: gene loss, gain, duplication, block
splitting and block fusion are frequently observed. 

ROAGUE accepts a set of species and a gene block in a reference species. It then finds all gene blocks, orhtologous to the reference gene blocks, and reconsructs their
ancestral states.

## Requirements
* [Python 3+](https://www.python.org/download/releases/3.0/)
* [Biopython 1.63+](http://biopython.org/wiki/Download)
* [Muscle Alignment](https://www.drive5.com/muscle/downloads.htm)
* [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* [ETE3](http://etetoolkit.org/download/) (python framework for tree)

## Installation
User can either use github interface Download or type the following command in command line:
```bash
git clone https://github.com/nguyenngochuy91/Ancestral-Blocks-Reconstruction
```
For the requirements, everything but ete3 can be installed using the following command line:
```bash
sudo apt-get install python-biopython blast2 ncbi-blast+ muscle
```

For ete3, check installation instructions on this website: http://etetoolkit.org/download/

## Usage

The easiest way to run the project is to execute the script [roague](https://github.com/nguyenngochuy91/Ancestral-Blocks-Reconstruction/blob/master/roague.py). 

1. The users can run this script on the example data sets provided in directory [E_Coli](https://github.com/nguyenngochuy91/Ancestral-Blocks-Reconstruction/tree/master/E_Coli) and [B_Sub](https://github.com/nguyenngochuy91/Ancestral-Blocks-Reconstruction/tree/master/B_Sub). The two following command line will run [roague](https://github.com/nguyenngochuy91/Ancestral-Blocks-Reconstruction/blob/master/roague.py) on our 2 directories. The final results (pdf files of our ancestral reconstructions) are stored in `result/E_Coli/visualization` and `result/B_Sub/visualization` directory.
### E_Coli
```bash
./roague.py -g E_Coli/genomes/ -b E_Coli/gene_block_names_and_genes.txt -r NC_000913 -f E_Coli/phylo_order.txt -m global
```

### B_Sub
```bash
./roague.py -g B_Sub/genomes/ -b B_Sub/gene_block_names_and_genes.txt -r NC_000964 -f B_Sub/phylo_order.txt -m global
```
2. If the users wants to run the program on their own datasets, then they have to provide the following inputs:
  1. Directory that store all the genomes file to study in genbank format 
  2. Gene block text file that stores gene blocks in a reference species (this reference has to be in the genomes directory). The gene block format is tab delimited. The first column is the gene block name, then followed by the genes' name. For example, here is an example `gene_block_names_and_genes.txt` from Escheria coli K-12 MG1655.
```bash
astCADBE	astA	astB	astC	astD	astE
atpIBEFHAGDC	atpI	atpH	atpC	atpB	atpA	atpG	atpF	atpE	atpD
caiTABCDE	caiA	caiE	caiD	caiC	caiB	caiT
casABCDE12	casE	casD	casA	casC	casB	cas1	cas2
chbBCARFG	chbG	chbF	chbC	chbB	chbA	chbR
``` 
Besides, the user can also provide a filter text file. This filter file specifies the species to be included in the reconstruction analysis. The reason is that there might be families of species that are over representative. This will reduce phylogenetic diversity and cause bias in our ancestral reconstruction. Hence, it is recomended to run [PDA](http://www.cibiv.at/software/pda/#download) on generated tree before proceeding further steps in our analysis. In order to achieve this, the user can follow the following instructions:
   1. Generate a phylogenetic trees from the genomes directory
   ```bash
   ./create_newick_tree.py -G genomes_directory -o tree_directory -f NONE -r ref_accession
   ```
   usage: create_newick_tree.py [-h] [-G DIRECTORY] [-o DIRECTORY] [-f FILE]
                             [-m STRING] [-t FILE] [-r REF] [-q]

  optional arguments:
    -h, --help            show this help message and exit
    
    -G DIRECTORY, --genbank_directory DIRECTORY
                          Folder containing all genbank files for use by the
                          program.
                          
    -o DIRECTORY, --outfolder DIRECTORY
                          Directory where the results of this program will be
                          stored.
                          
    -f FILE, --filter FILE
                          File restrictiong which accession numbers this script
                          will process. If no file is provided, filtering is not
                          performed.

    -r REF, --ref REF     The reference genome number such as NC_000913 for E_Coli 


   2. Download and install [PDA](http://www.cibiv.at/software/pda/#download). Debias the phylogenetic tree using `PDA` program:
   ```bash
   ./debias.py -i tree_directory/out_tree.nwk -o pda_result.txt -s num -r ref_accession
   ```
   usage: debias.py [-h] [-i INPUT_TREE] [-o PDA_OUT] [-s TREE_SIZE] [-r REF]



  optional arguments:
  -h, --help            show this help message and exit
  
  -i INPUT_TREE, --input_tree INPUT_TREE
                        Input tree that we want to debias
                        
  -o PDA_OUT, --pda_out PDA_OUT
                        Output of pda to be store.
                        
  -s TREE_SIZE, --tree_size TREE_SIZE
                        Reduce the size of the tree to this size, for example,
                        you can reduce your number of species from 100 to 30
                        by input 30
                        
  -r REF, --ref REF     Force to include the following species, here I force
                        to include the reference species

   3. Run [roague](https://github.com/nguyenngochuy91/Ancestral-Blocks-Reconstruction/blob/master/roague.py)
  ```bash
  ./roague.py -g genomes_directory -b gene_block_names_and_genes.txt -r NC_000964 -f phylo_order.txt -m global
  ```

roague.py [-h] [--genomes_directory GENOMES_DIRECTORY]
                 [--gene_blocks GENE_BLOCKS] [--reference REFERENCE]
                 [--filter FILTER] [--method METHOD]

Optional arguments:

  -h, --help            show this help message and exit

  --genomes_directory GENOMES_DIRECTORY, -g GENOMES_DIRECTORY
                        The directory that store all the genomes file in genbank format.
                        (`E_Coli/genomes`)

  --gene_blocks GENE_BLOCKS, -b GENE_BLOCKS
                        The gene_block_names_and_genes.txt file, this file
                        stores the operon name and its set of genes.

  --reference REFERENCE, -r REFERENCE
                        The ncbi accession number for the reference genome.
                        (NC_000913 for E_Coli and NC_000964 for B_Sub)

  --filter FILTER, -f FILTER
                        The filter file for creating the tree.
                        (`E_Coli/phylo_order.tx` for E_Coli or
                        `B_Sub/phylo_order.txt` for B_Sub)

  --method METHOD, -m METHOD
                        The method to reconstruct ancestral gene block, we
                        support either global or local approach.

## Examples

Here are two gene blocks that were generated through our program. 
1. Gene block paaABCDEFGHIJK:

This gene block codes for genes involved in the catabolism of phenylacetate and it is not conserved between the group of studied bacteria.

![paaABCDEFGHIJK](https://github.com/nguyenngochuy91/Ancestral-Blocks-Reconstruction/blob/master/image/paa_global.jpg "Gene block paaABCDEFGHIJK")
2. Gene block atpIBEFHAGDC:

This gene block catalyzes the synthesis of ATP from ADP and inorganic phosphate and it is very conserved between the group of studied bacteria.

![atpIBEFHAGDC](https://github.com/nguyenngochuy91/Ancestral-Blocks-Reconstruction/blob/master/image/atp_global.jpg "Gene block atpIBEFHAGDC")
## Credits
1. http://bioinformatics.oxfordjournals.org/content/early/2015/04/13/bioinformatics.btv128.full 



