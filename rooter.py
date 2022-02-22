#!/usr/bin/python3

import argparse, dendropy, subprocess

parser = argparse.ArgumentParser()
parser.add_argument("--treefile")
parser.add_argument("--gene")

args = parser.parse_args()
treefile = str(args.treefile)
gene = args.gene

tree = dendropy.Tree.get(path=treefile, schema="newick")
tips = str(tree.taxon_namespace)

#Rooting tree on outgroup species. 
# 1079, Aphandra natalia
# 1080, Dypsis ambositrae
# 1081, Eugeissona tristis
# 1082, Nypa fructicans

if "1079" in tips and "1080" in tips and "1081" in tips and "1082" in tips:
    cmd = "pxrr -t "+treefile+" -g 1079,1080,1081,1082 -o {gene}_rooted.tre".format(gene=gene)
else:
    if  "1079" in tips and "1080" in tips and "1081" in tips: 
        cmd = "pxrr -t "+treefile+" -g 1079,1080,1081 -o {gene}_rooted.tre".format(gene=gene)
    elif "1079" in tips and "1080" in tips and "1082" in tips:
        cmd = "pxrr -t "+treefile+" -g 1079,1080,1082 -o {gene}_rooted.tre".format(gene=gene)
    elif "1079" in tips and "1081" in tips and "1082" in tips:
        cmd = "pxrr -t "+treefile+" -g 1079,1081,1082 -o {gene}_rooted.tre".format(gene=gene)
    elif "1080" in tips and "1081" in tips and "1082" in tips:
        cmd = "pxrr -t "+treefile+" -g 1080,1081,1082 -o {gene}_rooted.tre".format(gene=gene)
    elif "1079" in tips and "1080" in tips:
        cmd = "pxrr -t "+treefile+" -g 1079,1080 -o {gene}_rooted.tre".format(gene=gene)
    elif "1079" in tips and "1081" in tips:
        cmd = "pxrr -t "+treefile+" -g 1079,1081 -o {gene}_rooted.tre".format(gene=gene)
    elif "1079" in tips and "1082" in tips:
        cmd = "pxrr -t "+treefile+" -g 1079,1082 -o {gene}_rooted.tre".format(gene=gene)
    elif "1080" in tips and "1081" in tips:
        cmd = "pxrr -t "+treefile+" -g 1080,1081 -o {gene}_rooted.tre".format(gene=gene)
    elif "1080" in tips and "1082" in tips:
        cmd = "pxrr -t "+treefile+" -g 1080,1082 -o {gene}_rooted.tre".format(gene=gene)
    elif "1081" in tips and "1082" in tips:
        cmd = "pxrr -t "+treefile+" -g 1081,1082 -o {gene}_rooted.tre".format(gene=gene)
    elif "1079" in tips:
        cmd = "pxrr -t "+treefile+" -g 1079 -o {gene}_rooted.tre".format(gene=gene)
    elif "1080" in tips:
        cmd = "pxrr -t "+treefile+" -g 1080 -o {gene}_rooted.tre".format(gene=gene)
    elif "1081" in tips:
        cmd = "pxrr -t "+treefile+" -g 1081 -o {gene}_rooted.tre".format(gene=gene)
    elif "1082" in tips:
        cmd = "pxrr -t "+treefile+" -g 1082 -o {gene}_rooted.tre".format(gene=gene)
    else:
        cmd = "mv "+treefile+" /home/owrisberg/Coryphoideae/work_flow/11_tree_building/01_genetrees/genetrees_no_outgroup"
        raise ValueError("Outgroup not in ", treefile)

subprocess.call(cmd, shell=True)