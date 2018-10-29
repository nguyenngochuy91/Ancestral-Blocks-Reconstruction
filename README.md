# ROAGUE: **R**econstruction **o**f **A**ncestral **G**ene Blocks **U**sing **E**vents
## Purpose

ROAGUE is a tool to reconstruct ancestors of gene blocks in prokaryotic genomes. Gene blocks are genes co-located on the chromosome. In many cases, gene blocks are conserved between bacterial species, sometimes as operons, when genes are co-transcribed. The conservation is rarely absolute: gene loss, gain, duplication, block splitting and block fusion are frequently observed. 

ROAGUE accepts a set of species and a gene block in a reference species. It then finds all gene blocks, orhtologous to the reference gene blocks, and reconsructs their ancestral states.

## Requirements
*[Wget](https://www.gnu.org/software/wget/) 
* [Conda](https://conda.io/miniconda.html) (package manager so we don't have to use sudo)
* [Python 3+](https://www.python.org/download/releases/3.0/)
* [Biopython 1.63+](http://biopython.org/wiki/Download)
* [Clustalw](http://www.clustal.org/clustal2/#Download)
* [Muscle Alignment](https://www.drive5.com/muscle/downloads.htm)
* [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* [ETE3](http://etetoolkit.org/download/) (python framework for tree)
* [PDA](http://www.cibiv.at/software/pda/#download) (optional if you want to debias your tree base on Phylogenetic Diversity)

## Installation
Users can either use github interface Download button or type the following command in command line:
```bash
git clone https://github.com/nguyenngochuy91/Ancestral-Blocks-Reconstruction
```
Install Miniconda (you can either export the path everytime you use ROAGUE, or add it to the .bashrc file). Before using
the following command line, users will need to install [Wget](https://www.gnu.org/software/wget/)
```bash
wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O Miniconda-latest-Linux-x86_64.sh
bash Miniconda-latest-Linux-x86_64.sh -b -p ~/anaconda_ete/
export PATH=~/anaconda_ete/bin:$PATH;

```

Install Biopython and ete3 using conda (highly recommended install biopython with conda)
```bash
conda install -c bioconda biopython ete3
```
Install ete_toolchain for visualization
```bash
conda install -c etetoolkit ete_toolchain
```

Install BLAST, ClustalW, MUSCLE 
```bash
conda install -c bioconda blast clustalw muscle
```

For PDA, check installation instructions on this website: [PDA](http://www.cibiv.at/software/pda/#download)

## Usage

The easiest way to run the project is to execute the script [ROAGUE](https://github.com/nguyenngochuy91/Ancestral-Blocks-Reconstruction/blob/master/roague.py) which is inside the directory [Ancestral-Blocks-Reconstruction]. 

### Run on example datasets
The users can run this script on the example data sets provided in directory [E_Coli](https://github.com/nguyenngochuy91/Ancestral-Blocks-Reconstruction/tree/master/E_Coli) and [B_Sub](https://github.com/nguyenngochuy91/Ancestral-Blocks-Reconstruction/tree/master/B_Sub). The two following command lines will run [roague](https://github.com/nguyenngochuy91/Ancestral-Blocks-Reconstruction/blob/master/roague.py) on our 2 directories. The final results (pdf files of our ancestral reconstructions) are stored in `result/E_Coli/visualization` and `result/B_Sub/visualization` directory by default.
#### E_Coli
```bash
./roague.py -g E_Coli/genomes/ -b E_Coli/gene_block_names_and_genes.txt -r NC_000913 -f E_Coli/phylo_order.txt -m global
```

#### B_Sub
```bash
./roague.py -g B_Sub/genomes/ -b B_Sub/gene_block_names_and_genes.txt -r NC_000964 -f B_Sub/phylo_order.txt -m global
```
### Run on users' specific datasets
If the users wants to run the program on their own datasets, then they have to provide the following inputs:
  1. Directory that stores all the genomes file to study in genbank format 
  2. Gene block text file that stores gene blocks in a reference species (this reference has to be in the genomes directory). The gene block format is tab delimited. The first column is the gene block name, then followed by the genes' name. For example, here is the `gene_block_names_and_genes.txt` file from Escheria coli K-12 MG1655.
```bash
astCADBE	astA	astB	astC	astD	astE
atpIBEFHAGDC	atpI	atpH	atpC	atpB	atpA	atpG	atpF	atpE	atpD
caiTABCDE	caiA	caiE	caiD	caiC	caiB	caiT
casABCDE12	casE	casD	casA	casC	casB	cas1	cas2
chbBCARFG	chbG	chbF	chbC	chbB	chbA	chbR
``` 
   3. Run ROAGUE, the output is stored in directory `result`.
  ```bash
  ./roague.py -g genomes_directory -b gene_block_names_and_genes.txt -r ref_accession -m global -o result
  ```
  ```
  usage: roague.py [-h] [--genomes_directory GENOMES_DIRECTORY]
                 [--gene_blocks GENE_BLOCKS] [--reference REFERENCE]
                 [--filter FILTER] [--method METHOD] [--output OUTPUT]

optional arguments:
  -h, --help            show this help message and exit
  --genomes_directory GENOMES_DIRECTORY, -g GENOMES_DIRECTORY
                        The directory that store all the genomes file
                        (E_Coli/genomes)
  --gene_blocks GENE_BLOCKS, -b GENE_BLOCKS
                        The gene_block_names_and_genes.txt file, this file
                        stores the operon name and its set of genes
  --reference REFERENCE, -r REFERENCE
                        The ncbi accession number for the reference genome
                        (NC_000913 for E_Coli and NC_000964 for B_Sub)
  --filter FILTER, -f FILTER
                        The filter file for creating the tree
                        (E_Coli/phylo_order.txt for E_Coli or
                        B_Sub/phylo_order.txt for B-Sub)
  --method METHOD, -m METHOD
                        The method to reconstruc ancestral gene block, we
                        support either global or local
  --output OUTPUT, -o OUTPUT
                        Output directory to store the result

  ```
   
Besides, the users can also provide a filter text file. This filter file specifies the species to be included in the reconstruction analysis. The reason is that there might be families of species that are over representative in our genomes directory. This will reduce phylogenetic diversity and cause bias in our ancestral reconstruction. Hence, it is recomended to run [PDA](http://www.cibiv.at/software/pda/#download) on generated tree before proceeding further steps in our analysis. In order to achieve this, the user can follow the following instructions:
   1. Generate a phylogenetic tree from the genomes directory
   ```bash
   ./create_newick_tree.py -G genomes_directory -o tree_directory -f NONE -r ref_accession
   ```
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
  -r REF, --ref REF     The reference genome number, such as NC_000913 for E_Coli
  -q, --quiet           Suppresses most program text outputs.

   ```
   2. Download and install [PDA](http://www.cibiv.at/software/pda/#download). Debias the phylogenetic tree using `PDA` program:
   ```bash
   ./debias.py -i tree_directory/out_tree.nwk -o pda_result.txt -s num -r ref_accession
   ```
   ```
   usage: debias.py [-h] [-i INPUT_TREE] [-o PDA_OUT] [-s TREE_SIZE] [-r REF]


optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_TREE, --input_tree INPUT_TREE
                        Input tree that we want to debias
  -o PDA_OUT, --pda_out PDA_OUT
                        Output of pda to be store.
  -s TREE_SIZE, --tree_size TREE_SIZE
                        Reduce the size of the tree to this size
  -r REF, --ref REF     Force to include the following species, here I force
                        to include the reference species

   ```
   3. Run ROAGUE,  the output is stored in directory `result`. 
  ```bash
  ./roague.py -g genomes_directory -b gene_block_names_and_genes.txt -r ref_accession -f phylo_order.txt -m global -o result
  ```




## Examples

Here are two gene blocks that were generated through our program. 
1. Gene block paaABCDEFGHIJK:

This gene block codes for genes involved in the catabolism of phenylacetate and it is not conserved between the group of studied bacteria.

![paaABCDEFGHIJK](https://github.com/nguyenngochuy91/Ancestral-Blocks-Reconstruction/blob/master/images/paa_global_edit.png "Gene block paaABCDEFGHIJK")
2. Gene block atpIBEFHAGDC:

This gene block catalyzes the synthesis of ATP from ADP and inorganic phosphate and it is very conserved between the group of studied bacteria.

![atpIBEFHAGDC](https://github.com/nguyenngochuy91/Ancestral-Blocks-Reconstruction/blob/master/images/atp_global_edit.png "Gene block atpIBEFHAGDC")
## Credits
1. http://bioinformatics.oxfordjournals.org/content/early/2015/04/13/bioinformatics.btv128.full 



