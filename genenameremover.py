#!/usr/bin/python3

import argparse, dendropy

parser = argparse.ArgumentParser()
parser.add_argument("treefile")
parser.add_argument("genename")

args = parser.parse_args()
gene = str(args.genename)
treefile = str(args.treefile)

# Formatting the input files
tree = dendropy.Tree.get(path=treefile, schema="newick")
tips = str(tree.taxon_namespace)
tips = tips.strip("[") 
tips = tips.strip("]")
tips = tips.replace("'","")
tipslist = tips.split(", ")

#print(tree.taxon_namespace)
for name in tree.taxon_namespace:
   print(getattr(name,"label"))
   setattr(name,"label",getattr(name,"label").replace("-{}".format(gene)))

print(tree.taxon_namespace)

#Looping through each tip in the tree
#for j in range(len(tipslist)):
  # print(tipslist[j],j)

