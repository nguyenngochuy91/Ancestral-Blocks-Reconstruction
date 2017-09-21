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

a. The user can run this script on the example data sets provided in directory [E_Coli](https://github.com/nguyenngochuy91/Ancestral-Blocks-Reconstruction/tree/master/E_Coli) and [B_Sub](https://github.com/nguyenngochuy91/Ancestral-Blocks-Reconstruction/tree/master/B_Sub). The two following command line will run [roague](https://github.com/nguyenngochuy91/Ancestral-Blocks-Reconstruction/blob/master/roague.py) on our 2 directories. The final results (pdf files of our ancestral reconstructions) are stored in `result/E_Coli/visualization` and `result/B_Sub/visualization` directory.
### E_Coli
```bash
./roague.py -g E_Coli/genomes/ -b E_Coli/gene_block_names_and_genes.txt -r NC_000913 -f E_Coli/phylo_order.txt -m global
```

### B_Sub
```bash
./roague.py -g B_Sub/genomes/ -b B_Sub/gene_block_names_and_genes.txt -r NC_000964 -f B_Sub/phylo_order.txt -m global
```
b. The user can provide his/her own set of bacteria species in genbank format. Preferably, the user should download all the available bacterial genomes from the NCBI website. The ftp link for that is [here](ftp://ftp.ncbi.nih.gov/genomes/archive/old_refseq/Bacteria/all.gbk.tar.gz). (about 30GB) 


Usage: 

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



