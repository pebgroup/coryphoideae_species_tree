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

#Removing the gene name from all the taxon labels
for name in tree.taxon_namespace:
   setattr(name,"label", getattr(name,"label").replace("-{}".format(gene),""))

#Defining a new TaxonNamespace with the edited names
newnames =dendropy.TaxonNamespace(tree.taxon_namespace)

#Adding the new namespace to the old tree
tree.taxon_namespace = newnames
tree.reconstruct_taxon_namespace()

#Writing the new tree
tree.write(path=treefile, schema="newick")


