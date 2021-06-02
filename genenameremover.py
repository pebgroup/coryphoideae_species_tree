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
tipslist = tips.split(", ")

print("tips",tips,"\n")
print("tipslist",tipslist,"\n")

#for j in range(len(tips)):
#   print(tips[j],j)
