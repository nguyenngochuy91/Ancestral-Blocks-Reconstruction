# ROAGUE: **R**econstruction **o**f **A**ncestral **G**ene Blocks **U**sing **E**vents
## Purpose

ROAGUE is a tool to reconstruct ancestors of gene blocks in prokaryotic genomes. Gene blocks are genes co-located on the chromosome. In many cases, gene blocks are
conserved between bacterial species, sometimes as operons, when genes are co-transcribed. The conservation is rarely absolute: gene loss, gain, duplication, block
splitting and block fusion are frequently observed. 

ROAGUE accepts a set of species and a gene block in a reference species. It then finds all gene blocks, orhtologous to the reference gene blocks, and reconsructs their
ancestral states.

## Requirements
* Python 2.7.6+
* Biopython 1.63+ (python-biopython)
* Muscle Alignment
* ncbi-tools (ncbi-tools-bin)
* BLAST2 (blast2)
* BLAST+ (ncbi-blast+)
* ete3 (python framework for trees)

## Installation
User can either use github interface Download or type the following command in command line:
```bash
git clone https://github.com/nguyenngochuy91/Ancestral-Blocks-Reconstruction
```
For the requirements, everything but ete3 can be installed using the following command line:
```bash
sudo apt-get install python-biopython ncbi-tools-bin blast2 ncbi-blast+ muscle
```

For ete3, user can check installation instructions on this website: http://etetoolkit.org/download/

## Usage

Caution:
There are several different between running python2 or python3 (such as divide function). My program is written in python3, if you want to run it correctly in python2, please uncooment the second line in the program "findParent_local.py".

The easiest way to run the project is to execute the script [roague](https://github.com/nguyenngochuy91/Ancestral-Blocks-Reconstruction/blob/master/roague.py). The user can run this script on the two data sets provided in directory [E.Coli](https://github.com/nguyenngochuy91/Ancestral-Blocks-Reconstruction/tree/master/E.Coli) and [B.Sub](https://github.com/nguyenngochuy91/Ancestral-Blocks-Reconstruction/tree/master/B.Sub). The two following command line will run [roague](https://github.com/nguyenngochuy91/Ancestral-Blocks-Reconstruction/blob/master/roague.py) on our 2 directories.
1. E.Coli: The final results (pdf file of our ancestral reconstruction) are stored in 'E.Coli/visualization' directory.
```bash
./roague.py -g E.Coli/genomes/ -b E.Coli/gene_block_names_and_genes.txt -r NC_000913 -f E.Coli/phylo_order.txt -m global
```
2. B.Sub: The final results (pdf file of our ancestral reconstruction) are stored in 'B.Sub/visualization' directory.
```bash
./roague.py -g B.Sub/genomes/ -b B.Sub/gene_block_names_and_genes.txt -r NC_000964 -f B.Sub/phylo_order.txt -m global
```

Each accompanying script can be run on its own as well, and each help for each script can be found by
using the -h or --help option.

```bash
./roague.py -h
```

usage: 
roague.py [-h] [--genomes_directory GENOMES_DIRECTORY]
                 [--gene_blocks GENE_BLOCKS] [--reference REFERENCE]
                 [--filter FILTER] [--method METHOD]

optional arguments:

  -h, --help            show this help message and exit

  --genomes_directory GENOMES_DIRECTORY, -g GENOMES_DIRECTORY
                        The directory that store all the genomes file
                        (E.Coli/genomes)

  --gene_blocks GENE_BLOCKS, -b GENE_BLOCKS
                        The gene_block_names_and_genes.txt file, this file
                        stores the operon name and its set of genes

  --reference REFERENCE, -r REFERENCE
                        The ncbi accession number for the reference genome
                        (NC_000913 for E.Coli and NC_000964 for B.Sub)

  --filter FILTER, -f FILTER
                        The filter file for creating the tree
                        (E.Coli/phylo_order.txt for E.Coli or
                        B.Sub/phylo_order.txt for B.Sub)

  --method METHOD, -m METHOD
                        The method to reconstruc ancestral gene block, we
                        support either global or local

## Credits
1. http://bioinformatics.oxfordjournals.org/content/early/2015/04/13/bioinformatics.btv128.full 



