# TRACCER
TRACCER detects convergent rate shifts in conjunction with a trait of interest while correcting for phylogentic distance and common tree-making artifacts.

It is generally used on gene or protein trees, but it equivalently works on arbitrary non-coding regions.

The significance of each gene tree is based on the distribution of scores from randomly permuting gene trees.
Specifically, branches defined by the extant lineages that share it as an ancestor are shuffled into synthetic trees.
Millions of these are scored, and actual gene scores are compared to this distribution to determine significance.

An in-depth description of these approaches are at
https://academic.oup.com/mbe/advance-article/doi/10.1093/molbev/msab226/6330627 [TODO update link after advace-access is updated to full]

Requirements:
1) Master tree. A single newick tree on the top line.
2) Fasta-like file of multiple trees, all fixed to the species phylogeny, but individual lineages can be missing.
3) List of species that have a trait of interest
4) Python 3, with numpy (https://numpy.org/install/) and ete3 (http://etetoolkit.org/download/) installed.

To use: python3 TRACCER.py --mastertree=pathtospeciestree --gtrees=pathtomultipletreefile --hastrait=speciesA,speciesB,speciesC --outgroup=speciesG,speciesH --cpus=numberofcpus
 
Species names must match those on the trees.

Outgroup species define one side of the species tree.

For more details: python3 TRACCER.py -help
