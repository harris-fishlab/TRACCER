# TRACCER
TRACCER detects convergent evolutionary rate shifts while correcting for phylogenetic relatedness, sampling biases, and common tree-making artifacts.

It is generally used on a genome-wide set of gene trees or protein trees, but it equivalently works on arbitrary genomic regions.

TRACCER calculates the background rate of convergence for your specific phylogeny by generating synthetic trees in which each branch is replaced by the same branch from a random tree in the dataset.
Millions of these synthetic trees are scored for convergent rate shifts to determine the background score distribution, then the actual tree scores are compared to this distribution to determine their individual significance. This permutation strategy can be computationally demanding, but it ensures the final scores are representative of the actual odds of these convergent patterns appearing in your specific  biological data and phylogeny.

An in-depth description of these approaches are at
https://doi.org/10.1093/molbev/msab226

Requirements:
1) Species tree. A single newick tree on the top line, with branch lengths representing relatedness.
2) File of multiple newick trees (thousands, genome-wide) all fixed to the species phylogeny, but individual lineages can be missing from individual trees. The format should be (\t means tab).
   
          unique_gene_name_1\t(newick_tree)
          unique_gene_name_2\t(newick_tree)
          
          or
          
          >unique_gene_name_1
          (newicktree)
          >unique_gene_name_2
          (newicktree)
          
3) List of species that have a trait of interest
4) Python 3, with numpy (https://numpy.org/install/) and ete3 (http://etetoolkit.org/download/) installed.

To use: python3 TRACCER.py --mastertree=path/to/speciestree --gtrees=path/to/multipletreefile --hastrait=speciesA,speciesB,speciesC --outgroup=speciesG,speciesH --cpus=numberofcpus

For more details and other flags: python3 TRACCER.py --help

Reach out to stephen_treaster (at) hms (dot) harvard (dot) edu for help. Still optimizing the user experience and would love feedback!
