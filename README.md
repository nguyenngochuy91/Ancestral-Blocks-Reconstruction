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
1. E.Coli
```bash
./roague.py -g E.Coli/genomes/ -b E.Coli/gene_block_names_and_genes.txt -r NC_000913 -f E.Coli/phylo_order.txt -m global
```
2. B.Sub 
```bash
./roague.py -g B.Sub/genomes/ -b B.Sub/gene_block_names_and_genes.txt -r NC_000964 -f B.Sub/phylo_order.txt -m global
```

Each accompanying script can be run on its own as well, and each help for each script can be found by
using the -h or --help option.


## Credits
1. http://bioinformatics.oxfordjournals.org/content/early/2015/04/13/bioinformatics.btv128.full 



