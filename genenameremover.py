#!/usr/bin/python3

import argparse, dendropy

parser = argparse.ArgumentParser()
parser.add_argument("treefile")
parser.add_argument("genename")

args = parser.parse_args()
gene = str(args.genename)
treefile = str(args.treefile)

# Go through each tip of the gene tree
tree = dendropy.Tree.get(path=treefile, schema="newick")
tips = str(tree.taxon_namespace)

print(tips)
class(tips)

#for name in tips:
    #print(name)
