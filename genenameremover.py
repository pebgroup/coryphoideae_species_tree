#!/usr/bin/python3

import argparse, dendropy

parser = argparse.ArgumentParser()
parser.add_argument("genename")
args = parser.parse_args()

gene = str(args.genename)

# Go through each tip of the gene tree
tree = dendropy.Tree.get(path=treefile, schema="newick")
tips = str(tree.taxon_namespace)

for name in tips:
    print(name)