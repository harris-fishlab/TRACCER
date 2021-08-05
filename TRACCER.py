#TRACCER detects convergent rate shifts in conjunction with a trait of interest while correcting for phylogentic distance, sampling biases, and common tree-making artifacts.
#It is generally used on gene or protein trees, but it equivalently works on arbitrary non-coding regions.
#The significance of each tree is based on the distribution of scores from randomly permuting branches across all trees.
#Specifically, branches defined by the extant lineages that share it as an ancestor are shuffled into a synthetic tree. Millions of these are scored, and actual tree scores are compared to this distribution to determine significance.
#An in depth description of these approaches are at https://academic.oup.com/mbe/advance-article/doi/10.1093/molbev/msab226/6330627?login=true [TODO update link after advace-access is updated to full]
#Requirements:
#1. Species tree. A single newick tree on the top line
#2. Fasta-like file of multiple newick trees, all fixed to the species phylogeny, but individual lineages can be missing.
#3. List of species bearing a trait of interest.
#4. Python 3, with numpy (https://numpy.org/install/) and ete3 (http://etetoolkit.org/download/)
#Soon to require pandas (https://pandas.pydata.org/)
#TRACCER.py --help for descriptions of the settings
#Reach out to stephen_treaster (at) hms (dot) harvard (dot) edu for help. Still optimizing the user experience and would love feedback!

import sys
import numpy as np
import random
import math
import os
import bisect
import timeit
import psutil
from ete3 import Tree
from multiprocessing import Pool
from optparse import OptionParser
from collections import defaultdict
from itertools import product,chain,cycle,combinations

parser = OptionParser()
parser.add_option("--mastertree", dest="mastertree_path", action='store',type='string',
        help="File with the species phylogeny in newick format on the first line. --mastertree=pathtomastertree")
parser.add_option("--gtrees", dest="genetrees_path", action='store',type='string',
        help="path to genetrees in fasta-like-newick format. --gtrees=pathtogtrees")
parser.add_option("--outgroup", dest="outgroup", action='store',type='string',default=False,
        help="comma separated list representing one side of the species tree, eg. --outgroup=speciesX,speciesY,speciesZ")
parser.add_option("--hastrait", dest="groupA", action='store',type='string',
        help="comma separated list of species that have the trait of interest, eg. --hastrait=speciesA,speciesB,etc")
parser.add_option("--cpus", dest="cpus", action='store',type='int', default=1,
        help="Number of cpus to use. More is faster, should use as many as you have, but requires more memory. --cpus=12")
parser.add_option("--maxpermutations", dest="maxpermutations", action='store',type='int', default=2000000,
        help="Permutations are computationally intensive. This is dynamically calculated based on time allotment, but you can restrict the maximum per cycle to save memory.")
parser.add_option("--time", dest="time", action='store',type='float', default=1.5,
        help="Time in hours to finish. Recommend ~10 with 12 cpus. High numbers can inflate memory requirements, leading to invisible crashing and timeouts. Start low, increase if permutations are too low for sensitivity. --time=H")
parser.add_option("--outname", dest="outname", action='store',type='string', default='Unnamed',
        help="name of the output result files eg. --outname=experiment1")
parser.add_option("--scalar_transform", dest="scalar_transform", action='store',type='string',default='rank',
        help="The Nth root transformation applied to the pairwise scalars, weighing each comparison by distance on the species phylogeny. Default is a ranking transformation. Accepts integers. --scalar_transform=2 means squareroot, etc.")
parser.add_option("--min_species", dest="min_species", action='store',type='float', default=.5,
        help="Discard gene trees with fewer extant lineages than this. Can be a fraction (eg, 0.5) of all species on the master tree, or a flat number (integer) cutoff.")
parser.add_option("--min_Tspecies", dest="min_Tspecies", action='store',type='float', default=.5,
        help="Discard gene trees with fewer TRAIT lineages than this. Can be a fraction (eg, 0.5) of all species on the master tree, or a flat number (integer) cutoff.")
parser.add_option("--min_NTspecies", dest="min_NTspecies", action='store',type='float', default=.5,
        help="Discard gene trees with fewer NONTRAIT lineages than this. Can be a fraction (eg, 0.5) of all species on the master tree, or a flat number (integer) cutoff.")
parser.add_option("--lackstrait", dest="groupB", action='store',type='string', default=None,
        help="Comma separated list of species that do not have the trait of interest. eg. --lackstrait=speciesC,speciesD,etc. This is NOT REQUIRED unless you are trying to prune species off the trees by leaving them out of both --lackstrait and --hastrait. --skipped can do the same thing.")
parser.add_option("--skipped", dest="groupC", action='store',type='string', default=None,
        help="Comma separated list of species to prune off of all trees and not be included in the analysis")
parser.add_option("--random", dest="random", action='store_true',
        help="Use random species selections of trait and non-trait; should still provide --hastrait=species to prevent randomly creating actual trait list. Can be used as a rough control run.")
parser.add_option("--sensitivity", dest="sensitivity", action='store_true',
        help="Run in sensitivity diagnostic mode, checking the robustness of the species phylogeny. Genetrees not required. STILL UNDER DEVELOPMENT FOR PUBLIC USE.")
options, args = parser.parse_args()

def help_text(*args):
    if args: print('\n'.join(args))
    print("\nTRACCER correlates relative evolutionary rates with a phenotype of interest and corrects for phylogenetic relatedness.\n")
    print("To use:\npython3 TRACCER.py --gtrees=allgenetrees_path --mastertree=mastertree_path --hastrait=speciesA,speciesC,speciesF --outname=myresultsname --outgroup=SpeciesG,SpeciesH\n")
    print("Use --help for more details on the flags available.\n")
    print("Genetrees must be in a fasta-like-newick format, like so:\n>genename1\n(newickstring);\n>genename2\n(newickstring);\n")
    print("The mastertree file is a single newick tree on the top line with all lineages included.")
    print("Outgroup species are required to parse ancestry when some lineages are missing on individual gene trees. List all the species, comma separated, no spaces, from one side of the root.")
    print("TRACCER requires the numpy and ete3 libraries installed and accessible to Python3.\n")
    print("""TRACCER is RAM intensive, and currently will stall indefinitely during a memory shortage.
If it does not update its progress at regular intervals, it likely ran out of memory.
To avoid these problems, use lower durations, eg. --time=3 , then slowly increase after successful runs.
~28GB should cover roughly 30,000 trees with ~60 species and 12 cpus for 10 hours, providing 10M permutations.
We are working to improve the resource efficiency and to set maximum memory usage directly.""")
    sys.exit()

if not options.groupA: help_text("\nMISSING REQUIRED FLAG '--hasttrait=speciesA,speciesB,speciesC'")
if not options.mastertree_path: help_text("\nMISSING REQUIRED FLAG '--mastertree=path_to_speciestree.nh'")
if not options.genetrees_path and not options.sensitivity: help_text("\n MISSING REQUIRE FLAGS, either '--gtrees=path_to_genetrees.fa' or '--sensitivity' (for a diagnostic on phylogeny convergence robustness)")
if not options.outgroup: help_text("\nMISSING REQUIRED FLAG '--outgroup=speciesX,speciesY,speciesZ' (list all terminal species on one side of the root)")
availmem = psutil.virtual_memory().available / (1024**3)
if availmem < 30: print("Available memory is less than 30GB, which is roughly what is required for a dataset of 30,000 trees, with 62 species, running for 10 hours with 12 cpus. If memory runs out, TRACCER will stall without any feedback.")
if options.cpus < 6: print("Highly recommend increasing --cpus to 6+. To get high sensitivity (10M+ permutations) requires 12 cpus for 10 hours on 30,000 tree with 62 mammals and 7000 unique topologies. These numbers are dependant on the number of unique topologies.")

#ADAPED from scipy.stats.rankdata to minimize importing non-standard libraries.
def rankdata(a, method='average', *, axis=None):
    if method not in ('average'):
        raise ValueError('unknown method "{0}": note this is a truncated scipy function'.format(method))
    if axis is not None:
        a = np.asarray(a)
        if a.size == 0:
            # The return values of `normalize_axis_index` are ignored.  The
            # call validates `axis`, even though we won't use it.
            # use scipy._lib._util._normalize_axis_index when available
            np.core.multiarray.normalize_axis_index(axis, a.ndim)
            dt = np.float64 if method == 'average' else np.int_
            return np.empty(a.shape, dtype=dt)
        return np.apply_along_axis(rankdata, axis, a, method)
    arr = np.ravel(np.asarray(a))
    sorter = np.argsort(arr, kind='quicksort')
    inv = np.empty(sorter.size, dtype=np.intp)
    inv[sorter] = np.arange(sorter.size, dtype=np.intp)
    arr = arr[sorter]
    obs = np.r_[True, arr[1:] != arr[:-1]]
    dense = obs.cumsum()[inv]
    # cumulative counts of each unique value
    count = np.r_[np.nonzero(obs)[0], len(obs)]
    return .5 * (count[dense] + count[dense - 1] + 1)

def initialize_species_groups():
    with open(options.mastertree_path) as infile:
        for row in infile:
            if row.startswith('>'): continue
            if '\t' in row: row = row.split('\t')[-1]
            try:
                mastertree = Tree(row)
            except:
                help_text("MASTER SPECIES TREE MISFORMATTED. Confirm it fulfills newick requirements.")
            break
    allspecies = frozenset([x.name for x in mastertree.get_leaves()])
    groupA = set(options.groupA.split(','))
    if options.random:
        halfA = set(random.sample(allspecies,int(len(groupA)/2)))
        halfB = set(random.sample(allspecies,len(groupA) - len(halfA)))
        groupA = halfA|halfB
    #Option to define groupB (lacks trait), or groupC (skipped), such that some species can be left out.
    if options.groupB and options.groupC:
        help_text("USE EITHER --lackstrait or --skipped. NOT BOTH")
    if options.groupB:
        groupB = set(options.groupB.split(','))
        skippedspecies = allspecies - (groupA | groupB)
        allspecies = frozenset(groupA | groupB)
        mastertree.prune(list(allspecies))
    elif options.groupC:
        groupC = set(options.groupC.split(','))
        allspecies = frozenset(allspecies - groupC)
        skippedspecies = list(groupC)
        groupB = allspecies - groupA
        mastertree.prune(list(allspecies))
    else:
        groupB = allspecies - groupA
        skippedspecies = ['none']
    print("Has Trait: %s" %','.join(groupA))
    print("Lacks Trait: %s" %','.join(groupB))
    print("Skipped species: %s" %','.join(skippedspecies))
    return groupA,groupB,allspecies,mastertree

def process_outgroup():
    outgroup = set(options.outgroup.split(','))
    outgroup = outgroup&allspecies
    dummy = next(iter(allspecies-outgroup)) #because of weird ete3 bug 'cannot set self as outgroup'
    mastertree.set_outgroup(mastertree&dummy) #so pick a random not-outgroup terminal branch and set it.
    if len(outgroup) > 1: #then reset it to the outgroup ancestor
        rootanc = mastertree.get_common_ancestor(*outgroup)
        mastertree.set_outgroup(rootanc)
    else: #or set it to just the outgroup if only one lineage is left
        outgroupStr = next(iter(outgroup))
        mastertree.set_outgroup(mastertree&outgroupStr)
    outgroup_order = [outgroup] #user chosen outgroup is first; then walk through rooted tree
    #each possible outgroup is needed, in traverse order, in case lineages are missing from gene trees
    #so they can be correctly rooted
    for node in mastertree.traverse():
        nextoutgroup = set([x.name for x in node.get_descendants() if x.name != ''])
        if len(nextoutgroup) == 0 or len(nextoutgroup) == len(allspecies): continue
        outgroup_order.append(nextoutgroup)
    return mastertree,outgroup_order

def ddset():
    return defaultdict(set)
def determine_mastersets():
    #branches can be identified by the extant lineages that share it, called species sets
    #setpaths are the path of such branches between pairs of extant species
    master_setpaths = defaultdict(ddset)
    master_setdists = {}
    for A,B in chain(product(groupA,groupB),product(groupB,groupA)):
        ancestor = mastertree.get_common_ancestor(A,B)
        tempnode = mastertree.get_leaves_by_name(name=A)[0]
        while tempnode != ancestor:
            nodeleaves = frozenset([x.name for x in tempnode.get_leaves()])
            master_setpaths[A][B].add(nodeleaves)
            master_setdists[nodeleaves] = tempnode.dist
            tempnode = tempnode.up
    return master_setpaths,master_setdists

def rescale_min_species():
    if options.min_species <= 1:   options.min_species = int(len(allspecies)*options.min_species)
    else:                         options.min_species = int(options.min_species)
    if options.min_Tspecies <= 1:  options.min_Tspecies = int(len(groupA)*options.min_Tspecies)
    else:                         options.min_Tspecies = int(options.min_Tspecies)
    if options.min_NTspecies <= 1: options.min_NTspecies = int(len(groupB)*options.min_NTspecies)
    else:                         options.min_NTspecies = int(options.min_NTspecies)
    if options.genetrees_path:
        print("minimum species in trees:  %d / %d" %(options.min_species,len(allspecies)))
        print("minimum Trait species:     %d / %d" %(options.min_Tspecies,len(groupA)))
        print("minimum Non-Trait species: %d / %d" %(options.min_NTspecies,len(groupB)))
    return options

def calculate_scalars(gtree,groupA,groupB,master_dists):
    longestdist = max([master_dists[A][B]+master_dists[B][A] for A,B in product(groupA,groupB)])
    scalars = defaultdict(dict)
    scalarlist = []
    for A,B in product(groupA,groupB):
        #from mastertree, calculate distance based scalar for each pair, with each branch scaled by number of uses
        #start at leaf, go UP repeatedly until you hit common ancestor. Add branch length * species sharing it at each step.
        #multiplying by "species use" stretches the branch to scale for sampling biases
        tempscalar = 0
        anc = gtree.get_common_ancestor(A,B)
        node = gtree.get_leaves_by_name(A)[0]
        while node != anc: #walk up from A to anc
            tempscalar += node.dist*len(node.get_leaves())#*branch use
            node = node.up
        node = gtree.get_leaves_by_name(B)[0]
        while node != anc: #walk up from B to anc
            tempscalar += node.dist*len(node.get_leaves()) #A->B will be same as B->A, they are added together.
            node = node.up
        if options.scalar_transform != 'rank':
            tempscalar = longestdist/tempscalar
            tempscalar = tempscalar**(1/float(options.scalar_transform))
            scalars[A][B] = tempscalar
        else:
            scalars[A][B] = (longestdist/tempscalar)
            scalarlist.append(scalars[A][B])
        scalars[B][A] = scalars[A][B]
        #final scalar minimum is ~1, for the furthest comparison. Maximum is the longestcomparison/shortest
    if options.scalar_transform == 'rank':
        ranks = rankdata(scalarlist,method='average')
        for i,(A,B) in enumerate(product(groupA,groupB)):
            scalars[A][B] = ranks[i]
            scalars[B][A] = ranks[i]
    else:
        print("Using %s root scalar transformation" %(options.scalar_transform))
    #Check most important comparisons, which is used when doing sensitivity
    #TODO: actually enable this for users
    most_important_sisters = []
    for A in groupA:
        most_important_sisters.append(max(groupB, key=lambda x: scalars[A][x]))
    return scalars,most_important_sisters

def calc_signal_distribution(magchangesorted):
    #This transforms randomized fixed-tree with random-lengths into a name describing the ranked pattern.
    signaldistribution = []
    for i,(A,B) in enumerate(product(groupA,groupB)):
        Ai = magchangesorted.index(A)
        Bi = magchangesorted.index(B)
        signaldistribution.append((Ai-Bi)*master_scalar[A][B])
    return sum(signaldistribution)/len(signaldistribution)

def ancestor_distances(gtree,ggroupA,ggroupB):
    #TODO Refactor to pandas
    MRCAdistsSingle = defaultdict(dict)
    distlist = []
    for A,B in product(ggroupA,ggroupB):
        #Do NOT use get_common_ancestor OR get_distance very often. They are SLOW
        ancestor = gtree.get_common_ancestor(A,B)
        distA = gtree.get_distance(ancestor,A)
        distB = gtree.get_distance(ancestor,B)
        MRCAdistsSingle[A][B] = distA
        MRCAdistsSingle[B][A] = distB
    return MRCAdistsSingle

def calcGeneSetDists(gtree):
    setdists = {}
    for node in gtree.iter_descendants():
        nodeleaves = frozenset([x.name for x in node.get_leaves()])
        if nodeleaves in master_setdists:
            setdists[nodeleaves] = node.dist
    return setdists

def roottree(tree,gspecies):
    for outgroup in outgroup_order:
        if len(gspecies) > len(gspecies&outgroup) >= 1:#outgroup can't be everything and can't be nothing.
            outgroup = gspecies&outgroup
            if len(outgroup) == 1: #different syntax for setting a single vs common ancestor as root
                tree.set_outgroup(tree&next(iter(outgroup)))
            else:
                dummy = next(iter(gspecies-outgroup)) #because of weird ete3 bug 'cannot set self as outgroup'
                tree.set_outgroup(tree&dummy) #so pick a random not-outgroup terminal branch and set it.
                rootanc = tree.get_common_ancestor(*outgroup) #then rotate back to actual outgroup
                tree.set_outgroup(rootanc)
            break #Only use the first outgroup that works. They are ordered from the master tree
    return tree

def prep_tree_data(name_and_gtree):
    name,gtree = name_and_gtree #annoying, but multiple arguments don't map well through pool
    try: #in case tree error, misformated, duplicated species
        gtree = Tree(gtree)
        gspecies = set([x.name for x in gtree.get_leaves()])
        gtree.prune(list(allspecies&gspecies)) #remove extra species not in master tree
        gspecies = allspecies&gspecies
        gtree = roottree(gtree,gspecies)
        ggroupA = gspecies & groupA
        ggroupB = gspecies & groupB
        #TODO Refactor these error flags to be more descriptive and user friendly
        if len(ggroupB) < options.min_NTspecies: return [name,'lowNTspecies',len(ggroupB)]
        if len(ggroupA) < options.min_Tspecies: return [name,'lowTspecies',len(ggroupA)]
        if len(gspecies) < options.min_species: return [name,'lowSpecies',len(gspecies)]
        #Remove trees without many branch differences. ie, many zero length branches
        #these numbers were determined by the distribution of average lengths on the UCSC 100-Way Alignment
        alldists = [x.dist for x in gtree.get_descendants()]
        uniquebranchlengthratio = len(set(alldists)) / len(alldists)
        if uniquebranchlengthratio < .2: return [name,'branchestoosimilar',str(uniquebranchlengthratio)]
        shortbranchratio = len([x for x in alldists if x < .0001]) / len(alldists)
        if shortbranchratio > .75: return [name,'shortbranchratio',str(shortbranchratio)]
    except Exception as e:
        print("bad tree: %s: %s" %(name, e))
        return [name,'Preperror',str(e)]
    try:
        dists = ancestor_distances(gtree,ggroupA,ggroupB)
    except Exception as e:
        return [name,'ancDisterror',str(e)]

    geneSetDists = calcGeneSetDists(gtree)
    #save lengths now, because juicing missing lineages, once implemented, will disrupt them
    #score is calculated later, after median rescaling dists to RERs
    result =  {'groupAlen':len(ggroupA),
              'groupBlen':len(ggroupB),
              'uniquebranchlengthratio':uniquebranchlengthratio,
              'shortbranchratio':shortbranchratio}
    return {'name':name,'dists':dists,'result':result,'geneSetDists':geneSetDists,'ggroupA':ggroupA,'ggroupB':ggroupB}

def calcscores(bundle):
    score = calculate_tree_score(bundle['ggroupA'],bundle['ggroupB'],bundle['RERs'])
    bundle['result']['score'] = score
    return bundle

def processtrees(genetrees_path):
    names2trees = []
    usednames = set()
    skippablesymbols = {'#','\n',"''"}
    with open(genetrees_path) as infile:
        #Parse fasta-like newick trees into names and trees.
        for row in infile:
            if row[0] == '>':
                name = row[1:].strip()
                if name in usednames:
                    print("%s name occurs at least twice" %name)
                    raise Exception("Naming overlap. Unique-ify your tree names.")
            elif row[0] == '(':
                if name in usednames:
                    print("%s has multiple trees" %name)
                    raise Exception("Fasta-like misformmating. One name to multiple trees.")
                names2trees.append((name,row.strip()))
                usednames.add(name)
            elif row[0] in skippablesymbols: continue
            else:
                print(row)
                help_text("Skipping this row. Fasta-like format may have errors.")
    with Pool(processes=options.cpus) as p:
        print("Processing %d input trees... (~10 minutes for 30k trees with 12 cpus)" %len(names2trees))
        processedtrees = list(p.imap_unordered(prep_tree_data,names2trees))
    badbundles,goodbundles = [],[]
    for bundle in processedtrees:
        #TODO this structure is not intuitive. Refactor to remove the nested results dict.
        if type(bundle) is not list: goodbundles.append(bundle)
        else:                        badbundles.append(bundle)
    goodbundles,allmedians = calculate_RERs(goodbundles)
    if len(goodbundles) < 1000:
        print("Less than 1000 gene trees passed quality checks. There is simply not enough variation in this dataset to identify meaningful convergence or generate a meaningful score-space. TRACCER will try anyway, but it is only validated on proteome wide datasets (20,000+ trees).")

    print("scoring...")
    goodbundles = list(map(calcscores,goodbundles))
    print("Trees that passed quality checks: %d" %len(goodbundles))
    print("Excluded trees are noted with the reason in the notes file.")
    #flatten results
    results = {}
    treetypes = defaultdict(list)
    all_species_setdists = defaultdict(list)
    for bundle in goodbundles:
        results[bundle['name']] = bundle['result']
        treetypes[frozenset(list(bundle['dists'].keys()))].append(bundle['name'])
        for nodeleaves,dist in bundle['geneSetDists'].items():
            all_species_setdists[nodeleaves].append(dist) #to be used for perms

    #Sometimes no gene trees contain all species. Adding mastertree dists for each branch gives something for those branches to work with, but they will always be the same.
    for x in master_setdists: all_species_setdists[x].append(master_setdists[x])
    print("Clade representation is the number of trees with all members of a given clade")
    print("Poorly represented clades will have insufficient sampling to impart significance")
    print("Least Represented Clade: %d" %min([len(all_species_setdists[x]) for x in master_setdists]))
    print("Most Represented Clade:  %d" %max([len(all_species_setdists[x]) for x in master_setdists]))
    print("Average Representation:  %d" %np.average([len(all_species_setdists[x]) for x in master_setdists]))
    return results,badbundles,all_species_setdists,len(goodbundles),allmedians,treetypes

def calculate_RERs(goodtreebundles):
    #TODO Refactor to pandas
    #Due to tree making artifacts median RERs are often much lower than 1. This forces the median to be 1.
    #update all pair comparison by skipping short artifacts, calculating medians, dividing by medians, and return medians, as they are needed for calculating RERs in permutations later
    print('calculating RERs...')
    allABdists = defaultdict(lambda: defaultdict(list))
    for bundle in goodtreebundles: #compile ABdists across all trees to get medians
        for A,Bdict in bundle['dists'].items():
            for B,value in Bdict.items():
                if value < .0001:
                    continue #values below this are zeroes or artifacts of fixed trees and throw off median calcs
                allABdists[A][B].append(value)
    #check completeness, alert user if any comparisons are poorly represented in the dataset.
    for A in groupA:
        for B in groupB:
            if len(allABdists[A][B]) < 10:
                allABdists[A][B].append(master_dists[A][B]) #add masterdist to guarantee it isn't empty
                print("Poor quality comparison, few trees have it with good branches: %s %s"%(A,B))
            if len(allABdists[B][A]) < 10: #add masterdists to guarantee it isn't empty
                allABdists[B][A].append(master_dists[B][A])
                print("Poor quality comparison, few trees have it with good branches: %s %s"%(B,A))

    medians = defaultdict(dict)
    for A,Bdict in allABdists.items():
        for B,values in Bdict.items():
            medians[A][B] = np.median(values)
    for i,bundle in enumerate(goodtreebundles):
        for A,Bdict in bundle['dists'].items():
            for B,value in Bdict.items():
                newvalue = value/medians[A][B]
                #minimum RER is 1/100; anything slower should be guttered, as it is likely due to short-sequence/artifacts
                if   newvalue >= .01:  bundle['dists'][A][B] = math.log(newvalue)
                #else:                 bundle['dists'][A][B] = math.log(.0099)
                else:                  bundle['dists'][A][B] = -4.61522052184
        goodtreebundles[i]['RERs'] = bundle['dists']
    return goodtreebundles,medians

def calc_permutation_scores(permcount):
    with Pool(processes=options.cpus) as p:
        rawperms = p.map(permutation_worker,list(range(permcount)))
    #flatten, because rawperms is list of dicts with one score per treetype.
    flat_permutations = defaultdict(list)
    for permdict in rawperms:
        for treetype,score in permdict.items():
            flat_permutations[treetype].append(score)
    return flat_permutations

def permutation_worker(permcount):
    #Branches are defined by the set of lineages that share it.
    #For each branch of master tree, grab from random tree
    #Sample is faster than repeated choice, but the sets are not a constant size. If it chooses too large an i, use choice as a backup, which is slower, but relatively uncommon.
    permutation_setdict = {}
    randIs = random.sample(list(range(goodtreecount)),len(master_setdists))
    for i,(nodeleaves,dist) in enumerate(master_setdists.items()):
        try: #To allow incomplete trees to contibute, the number of trees with each nodeset is not constant.
            permutation_setdict[nodeleaves] = all_species_setdists[nodeleaves][randIs][i]
        except: #has to work, since it's specific for that nodeset
            permutation_setdict[nodeleaves] = random.choice(all_species_setdists[nodeleaves])
    permutation_dists = convert_setdists(permutation_setdict)
    permutation_score = {}
    for treetype in treetypes:
        ggroupA,ggroupB = groupA&treetype,groupB&treetype
        permutation_score[treetype] = calculate_tree_score(ggroupA,ggroupB,permutation_dists)
    return permutation_score

def generate_and_score_background_brownian_trees(permcount):
    with Pool(processes=options.cpus) as p:
        rawperms = p.map(background_worker,list(range(permcount)))
    rawperms = sorted(rawperms)
    permmedian = np.median(rawperms)
    permlow = [x for x in rawperms if x <= permmedian]
    permhigh = rawperms[len(permlow):]
    return rawperms,permlow,permhigh,permmedian

def background_worker(permcount):
    #Brownian motion all branches of master tree. Identical to the non-background function, but does not have foundation for juicing.
    permutation_setdict = dict()
    for nodeleaves,dist in master_setdists.items():
        permutation_setdict[nodeleaves] = dist * np.random.lognormal(mean=0,sigma=.1)
    #convert specsetdists to normal AB pairwise comparison, RER, and median correction
    permutation_dists = convert_setdists(permutation_setdict)
    permutation_score = calculate_tree_score(groupA,groupB,permutation_dists)
    return permutation_score

def convert_setdists(genesetdists):
    permutation_dists = defaultdict(dict)
    for A,B in product(groupA,groupB):
        newvalueAB = sum([genesetdists[x] for x in master_setpaths[A][B]])/allmedians[A][B]
        newvalueBA = sum([genesetdists[x] for x in master_setpaths[B][A]])/allmedians[B][A]
        if     newvalueAB >= .01: permutation_dists[A][B] = math.log(newvalueAB)
        else:                     permutation_dists[A][B] = math.log(.0099)
        if     newvalueBA >= .01: permutation_dists[B][A] = math.log(newvalueBA)
        else:                     permutation_dists[B][A] = math.log(.0099)
    return permutation_dists

def generate_and_score_brownian_trees_pool():
    #TODO: Fully implement juicing for users
    hastraitcombos = [frozenset(x) for x in hastraitcombos] #Foundation for juicing
    senscount = options.maxpermutations
    with Pool(processes=options.cpus) as p:
        results = p.starmap(generate_and_score_tree_worker,zip(range(senscount),cycle(hastraitcombos)))
    return results

def generate_and_score_tree_worker(poolID,juicedspecies):
    changemap = {}
    tempdists = {}
    for nodeleaves,dist in master_setdists.items():
        jitter = np.random.lognormal(mean=0,sigma=.1)
        if len(nodeleaves) == 1: #Terminal branches
            (singlespec,) = nodeleaves #unpacks single member sets
            jitter = np.random.lognormal(mean=0,sigma=.1) #currently terminal and internal branches are treated the same
            ''' #TODO: Add juiced trait-bearing species
            if singlespec in juicedspecies:
                jitter = np.random.uniform(2000000,2000001) #arbitrary large number,but not constant
            '''
            ''' #TODO: Add juiced sisters (non-trait bearing)
            if singlespec in ransisters:
                #jitter = 0 #jitter doesn't matter, will be treated as less than arts, and forced to 0
                jitter = np.random.uniform(2000000,2000001)
            '''
            changemap[singlespec] = jitter #set jitter for the terminal branches
        tempdists[nodeleaves] = dist * jitter #Change each branch, internal and terminal, by jitter

    RERs = convert_setdists(tempdists)
    score = calculate_tree_score(groupA,groupB,RERs)

    #Generate a name for this tree that roughly describes the pattern.
    #Sort jittered species by terminal rate of change
    magchangesorted = [x[0] for x in sorted(list(changemap.items()), key=lambda x: x[1])]
    avgAi = sum([magchangesorted.index(species) for species in groupA])/len(groupA)
    avgBi = sum([magchangesorted.index(species) for species in groupB])/len(groupB)
    avgrankdiff = avgAi-avgBi
    signaldistribution = calc_signal_distribution(magchangesorted)
    #Get the name right for later parsing
    juicedspecies = ','.join(juicedspecies)
    juiceID = [avgrankdiff,signaldistribution,poolID]
    juiceID = '|'.join([str(x) for x in juiceID])
    return (juiceID,score)

def calculate_tree_score(ggroupA,ggroupB,RERs):
    usedscalars = []
    comparisonscores = []
    test = defaultdict(list)
    #positives,negatives = [],[] #TODO: track discordant lineages
    for A,B in product(ggroupA,ggroupB):
        ABdiff = RERs[A][B]-RERs[B][A] #both present, regular
        #if ABdiff >= 0: positives.append(B) # track discordant lineages
        #elif ABdiff < 0: negatives.append(B)# track discoradant lineages
        comparisonscores.append(ABdiff)
        usedscalars.append(master_scalar[A][B])
    absscores = [abs(x) for x in comparisonscores]
    ranks = rankdata(absscores,method='average')
    polarity = [-1 if x < 0 else 0 if x == 0 else 1 for x in comparisonscores]
    ranksPolScal = [a*b*c for a,b,c in zip(ranks,polarity,usedscalars)]
    rankGeneScore = sum(ranksPolScal)
    return rankGeneScore

def determine_p(score,permlow,permhigh,permmedian):
    if score < permmedian:
        rankposition = bisect.bisect_left(permlow,score)
        rankpvalue = rankposition/len(permlow)
        if rankpvalue == 0: rankpvalue = 1/len(permlow)/2
    else:
        rankposition = bisect.bisect_right(permhigh,score)
        rankpvalue = 1-(rankposition/len(permhigh))
        if rankpvalue == 0: rankpvalue = 1/len(permhigh)/2
    return rankpvalue

def permutations_wrapper(results):
    #TODO Refactor to support pandas, which should drop memory costs as well.
    #TODO Scale settings by maximum memory.
    #Store all permutation scores for each treetype. Keep extending scores through each loop.
    #Check significance of every tree in that treetype each loop.
    #Repeat until a treetype no longer contains significant trees.
    #Record the significance of those trees.
    #Remove that treetype to free up memory; no longer need to check.
    #Make minimum 5000 permutations for each treetype.
    global treetypes #Treetypes are removed as they are exhausted of significant hits. Used in the pool workers.
    p_thresholds,pthreshold,minp = [], 0.5, 0.00001 #roughly proteome size cutoff for last cycle
    while pthreshold > minp:
        p_thresholds.append(pthreshold)
        pthreshold = pthreshold/2
    permsize = 2000 #starting permutations for first significance threshold
    master_permutation_scores = defaultdict(list)
    finalresults = []
    print("Calculating permutations and score distribution...")
    print("Increase --cpus=N and/or --time=H to increase sensitivity, but both increase memory usage.")
    print("Will stall without an error message if it runs out of memory")
    for i,pthreshold in enumerate(p_thresholds):
        if len(treetypes) == 0: break
        starttime = timeit.default_timer() #Track each loop time usage to rescale later loops
        if permsize < 100: permsize = 100
        if permsize > options.maxpermutations: permsize = options.maxpermutations
        print("Cycle: %d out of %d" %(i+1,len(p_thresholds)))
        print("Significance Threshold:            %1.6f"%pthreshold)
        print("Unique tree topologies remaining:  %d"%len(treetypes))
        print("Permuations for each topology:     %d"%permsize)
        print("Working... should use ~1/15th of the --time=%d setting. If it goes consideraly longer, it stalled."%options.time)
        permutation_scores = calc_permutation_scores(permsize)
        #TODO replace permuation scores with numpy array to save memory
        for treetype,scores in permutation_scores.items():
            #keep adding permutations to treetypes that retain a significant hit
            master_permutation_scores[treetype].extend(scores)
        #in case the permutation score distribution isn't centered on zero
        #make sure p values are based on distance to the median.
        #otherwise they will be artificially suppressed/enhanced by the opposite side.
        #and allows for different slope shapes on each side
        for treetype,treenames in treetypes.items():
            #TODO replace permutation scores with numpy array to save memory
            sortpermutation_scores = sorted(master_permutation_scores[treetype])
            permmedian = np.median(sortpermutation_scores)
            permlow = [x for x in sortpermutation_scores if x <= permmedian]
            permhigh = sortpermutation_scores[len(permlow):]
            tempfinalresults = []
            for name in treenames:
                pvalue = determine_p(results[name]['score'],permlow,permhigh,permmedian)
                if pvalue > pthreshold or i+1 == len(p_thresholds): #last cycle grab everything
                    if results[name]['score'] >= 0: pressure = "accelerated"
                    else:                           pressure = "constrained"
                    outresult = [name,
                            results[name]['groupAlen'],
                            results[name]['groupBlen'],
                            len(sortpermutation_scores),
                            pressure,
                            pvalue
                            ]
                    tempfinalresults.append(outresult) #store trees less significant than threshold, but don't commit them unless all are less significant.
            #Don't empty treetypes unless they are all done or last cycle
            donetreestoremove = set([x[0] for x in tempfinalresults])
            if len(donetreestoremove) == len(treenames) or i+1 == len(p_thresholds):
                treetypes[treetype] = [x for x in treetypes[treetype] if x not in donetreestoremove]
                finalresults.extend(tempfinalresults) #commit
        #remove empty (fully processed) tree types and permutations
        treetypes_lastlength = len(treetypes)  #for calculating permutations per time
        treetypes = {k:v for k,v in treetypes.items() if len(v) > 0}
        master_permutation_scores = {k:v for k,v in master_permutation_scores.items() if k in treetypes}
        #estimate time for permutations, allocate time per pthreshold
        timechange = timeit.default_timer() - starttime
        print("Minutes used for cycle: %1.1f"%(timechange/60))
        if i > 1: print("Estimated minutes remaining:       %1.1f" %(timechange/60*(len(p_thresholds)-i)))
        #number of perms to do, based on number done per second, and how much time to do it in.
        #needs to be divided by treetypes back at the top
        permsize = (options.time)*3600/len(p_thresholds) / (timechange/(permsize*treetypes_lastlength))
        if len(treetypes) > 0:
            permsize = int(permsize/len(treetypes))
        else:
            permsize = 0
        pthreshold = pthreshold/2
    return finalresults

def savenotes(badbundles):
    with open(options.outname+'.TRACCER.notes.txt', 'w') as out:
        out.write('Has trait:'+'\n')
        out.write(','.join(groupA)+'\n')
        out.write('Lacks trait:'+'\n')
        out.write(','.join(groupB)+'\n')
        for result in badbundles:
            out.write('\t'.join([str(x) for x in result]))
            out.write('\n')

def saveresults(finalresults):
    finalresults = sorted(finalresults,key=lambda x: x[5]) #sort by p value
    if finalresults[0][3] < len(finalresults)*100:
        print("Permutations on the most significant hits may be insufficient. Re-run with more --time and/or --cpus to improve sensitivity. This will use more memory however. If TRACCER runs out of memory, it will stall indefinitely.")
    headers = ['Tree name','Trait species','Non-trait species','Permutations','Selective Pressure','p-value','FDR']
    with open(options.outname+'.TRACCER.txt', 'w') as out:
        out.write('\t'.join(headers))
        out.write('\n')
        for i,result in enumerate(finalresults):
            if options.sensitivity: FDR = 25000*result[2] #p*rough proteome size
            else: FDR = len(finalresults)*result[5]/(i+1) #number by chance / number actual
            result.append(FDR)
            out.write('\t'.join([str(x) for x in result]))
            out.write('\n')

def main():
    #TODO: properly scope these into pool workers to minimize memory bloat.
    #Constants:
    global allspecies
    global groupA #local changes when species are missing
    global groupB #local changes when species are missing
    global master_scalar
    global master_dists
    global mastertree
    global outgroup_order
    global most_important_sisters
    global all_species_setdists
    global master_setpaths
    global master_setdists
    global goodtreecount
    global permmedian
    global allmedians
    #treetypes changes as they are depleted of significant hits
    global treetypes

    starttime = timeit.default_timer() #Track each loop time usage to rescale later loops
    groupA, groupB, allspecies, mastertree = initialize_species_groups()
    #TODO: subtract prep time from overall options.time
    mastertree,outgroup_order = process_outgroup()
    master_dists = ancestor_distances(mastertree,groupA,groupB)
    master_scalar,most_important_sisters = calculate_scalars(mastertree,groupA,groupB,master_dists)
    master_setpaths,master_setdists = determine_mastersets()
    options = rescale_min_species() #overwrite user fraction minimums with actual values
    if not options.sensitivity:
        results,badbundles,all_species_setdists,goodtreecount,allmedians,treetypes = processtrees(options.genetrees_path)
        savenotes(badbundles)
        timechange = timeit.default_timer() - starttime
        options.time = options.time - timechange/3600
        finalresults = permutations_wrapper(results)
    elif options.sensitivity:
        #overwrite treetypes and results. Make and score brownian motion trees based on species tree
        print("Running sensitivity simulation...")
        allmedians = master_dists
        results = generate_and_score_brownian_trees_pool() #TODO: add options for userjuicing
        print("Running sensitivity stats...")
        #Generate many brownian trees and score them.
        perms,permlow,permhigh,permmedian = generate_and_score_background_brownian_trees(options.maxpermutations)
        finalresults = []
        for name,score in results:
            pvalue = determine_p(score,permlow,permhigh,permmedian)
            outresult = [name,score,pvalue,options.maxpermutations]
            finalresults.append(outresult)
    else: help_text()

    saveresults(finalresults)
    print("Done. Results recorded to files prefixed with --outname.")

if __name__ == "__main__":
    main()
