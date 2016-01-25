# Ancestral-BLocks-Reconstruction project
## Sypnosis
This project provide ancestral reconstruction tools dedicated to bacteria genomes.
## Model Assumption
* Phylogenetic Tree
  1. Our phylogenetic tree is given, it is binary and rooted.
  2. Leaves are populated by orthoblocks. At least one leaf has a reference operon.
  3. The model is agnostic to gene order.
* Relationship between parent nodes and their children
  1. Given a parent gene blocks, its children gene blocks can't have any gene that is not in the parent gene blocks. (Hard assumption)
  2. There are 3 types of events that can happen from a parent to a child:
     * Split      : If two genes in one taxon are neighboring and their homologs in the other taxon are not, then that is defined as a single split event. The distance is the minimal number of split events identified between the compared genomes.
     * Deletion    : A gene exists in the operon in the one taxon, but its homolog cannot be found in an orthoblock in another taxon. Note that the definition of homolog, e-value 10−10 is strict, and may result in false negatives. The deletion distance is the number of deletion events identified between the compared target genomes.
     * Duplication : A duplication event is defined as having gene j in a gene block in the source genome, and homologous genes (j′,j″)(j′,j″) in the homologous block in the target genome. The duplication distance is the number of duplication events counted between the source and target genomes. The duplication has to occur in a gene block to be tallied.
  3. Multiple events from parent to children are possible.
  4. Events are treated as independent.
## Installation
TODO: Describe the installation process
## Usage
TODO: Write usage instructions
## History
TODO: Write history
## Credits
TODO: Write credits
## License



