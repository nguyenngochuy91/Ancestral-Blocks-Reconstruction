# ROAGUE: **R**econstruction **o**f **A**ncestral **G**ene Blocks **U**sing **E**vents
## Purpose

ROAGUE is a tool to reconstruct ancestors of gene blocks in prokaryotic genomes. Gene blocks are genes co-located on the chromosome. In many cases, gene blocks are
conserved between bacterial species, sometimes as operons, when genes are co-transcribed. The conservation is rarely absolute: gene loss, gain, duplication, block
splitting and block fusion are frequently observed. 

ROAGUE accepts a set of species and a gene block in a reference species. It then finds all gene blocks, orhtologous to the reference gene blocks, and reconsructs their
ancestral states.

## Prerequisites
* [http://etetoolkit.org/download/] ete3 

## Installation
User can either use github interface Download or type the following command in command line:
```bash
git clone https://github.com/nguyenngochuy91/Ancestral-Blocks-Reconstruction
```
User has to install python, and ete3 (would recommend use anaconda package). The instruction can be followed from this site:
http://etetoolkit.org/download/

## Usage

Caution:
There are several different between running python2 or python3 (such as divide function
). My program is written in python3, if you want to run it correctly in python2, please uncooment the second line in the program "findParent_local.py".
Here is the step by step usage:
* Step 1: Parse each operon file in the result dic, and convert it into appropriate form:
  1. For each gene in an operon file, map it into an alphabet letter
  2. For each genomes, use the start, stop position and strand to find neighbor gene pairs, only display those genomes info that actually have at least orthoblocks
  3. Use the command line below and the output will be stored in directory new_result
```bash
./convert.py  -i result_dic/ -o new_result/ 
```

* Step 2: Using each operon file from new result, and the tree input (this tree was built using muscle alignment of the 33 taxa on the rpOb gene marker), then provide ancestral reconstruction method depends on user choice (global or local)
 1. Use the command line below and the output will be stored in directory reconstruction. Method used is global
```bash
./reconstruction.py -i new_result/ -t muscle.ph -o reconstruction_global/ -m global 
```
* Step 3: Provide a visualization of the ancestral reconstruction using ete3. ete3 also provides a grouping theme depends on the class of the taxa. You can uncomment the line 103 to render the file into image, however, you need to provide the where to output the render file
  * Use the command line below for each operon that you like, here I use not so conserved operon paaABCDEFGHIJK:
```bash
./show_tree.py -g group.txt -i reconstruction_global/paaABCDEFGHIJK 
```

  ![Image of paaABCDEFGHIJK](https://github.com/nguyenngochuy91/Ancestral-Blocks-Reconstruction/blob/master/image/paa_global.jpg)
  * Use the command line below for each operon that you like, here I use a highly conserved operon rplKAJL-rpoBC:
```bash
./show_tree.py -g group.txt -i reconstruction_global/rplKAJL-rpoBC 
```
  ![Image of paaABCDEFGHIJK](https://github.com/nguyenngochuy91/Ancestral-Blocks-Reconstruction/blob/master/image/rpl_global.jpg)

  * Render the file:

```bash
./show_tree.py -g group.txt -i reconstruction_global/caiTABCDE -o caiTABCDE_image
```

## Credits
1. http://bioinformatics.oxfordjournals.org/content/early/2015/04/13/bioinformatics.btv128.full 



