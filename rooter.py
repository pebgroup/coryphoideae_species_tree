#!/usr/bin/python3

import argparse, dendropy, subprocess

parser = argparse.ArgumentParser()
parser.add_argument("treefile")
parser.add_argument("genename")

args = parser.parse_args()
treefile = str(args.treefile)
gene = str(args.genename)

tree = dendropy.Tree.get(path=treefile, schema="newick")
tips = str(tree.taxon_namespace)

if "1079-{}".format(gene) in tips and "1080-{}".format(gene) in tips and "1081-{}".format(gene) in tips and "1082-{}".format(gene) in tips:
    cmd = "pxrr -t "+treefile+" -g 1079-{},1080-{},1081-{},1082-{} -o temp.tre".format(gene)
else:
    if  "1079-{}".format(gene) in tips and "1080-{}".format(gene) in tips and "1081-{}".format(gene) in tips: 
        cmd = "pxrr -t "+treefile+" -g 1079-{},1080-{},1081-{} -o temp.tre".format(gene)
    elif "1079-{}".format(gene) in tips and "1080-{}".format(gene) in tips and "1082-{}".format(gene) in tips:
        cmd = "pxrr -t "+treefile+" -g 1079-{},1080-{},1082-{} -o temp.tre".format(gene)
    elif "1079-{}".format(gene) in tips and "1081-{}".format(gene) in tips and "1082-{}".format(gene) in tips:
        cmd = "pxrr -t "+treefile+" -g 1079-{},1081-{},1082-{} -o temp.tre".format(gene)
    elif "1080-{}".format(gene) in tips and "1081-{}".format(gene) in tips and "1082-{}".format(gene) in tips:
        cmd = "pxrr -t "+treefile+" -g 1080-{},1081-{},1082-{} -o temp.tre".format(gene)
    elif "1079-{}".format(gene) in tips and "1080-{}".format(gene) in tips:
        cmd = "pxrr -t "+treefile+" -g 1079-{},1080-{} -o temp.tre".format(gene)
    elif "1079-{}".format(gene) in tips and "1081-{}".format(gene) in tips:
        cmd = "pxrr -t "+treefile+" -g 1079-{},1081-{} -o temp.tre".format(gene)
    elif "1079-{}".format(gene) in tips and "1082-{}".format(gene) in tips:
        cmd = "pxrr -t "+treefile+" -g 1079-{},1082-{} -o temp.tre".format(gene)
    elif "1080-{}".format(gene) in tips and "1081-{}".format(gene) in tips:
        cmd = "pxrr -t "+treefile+" -g 1080-{},1081-{} -o temp.tre".format(gene)
    elif "1080-{}".format(gene) in tips and "1082-{}".format(gene) in tips:
        cmd = "pxrr -t "+treefile+" -g 1080-{},1082-{} -o temp.tre".format(gene)
    elif "1081-{}".format(gene) in tips and "1082-{}".format(gene) in tips:
        cmd = "pxrr -t "+treefile+" -g 1081-{},1082-{} -o temp.tre".format(gene)
    elif "1079-{}".format(gene) in tips:
        cmd = "pxrr -t "+treefile+" -g 1079-{} -o temp.tre".format(gene)
    elif "1080-{}".format(gene) in tips:
        cmd = "pxrr -t "+treefile+" -g 1080-{} -o temp.tre".format(gene)
    elif "1081-{}".format(gene) in tips:
        cmd = "pxrr -t "+treefile+" -g 1081-{} -o temp.tre".format(gene)
    elif "1082-{}".format(gene) in tips:
        cmd = "pxrr -t "+treefile+" -g 1082-{} -o temp.tre".format(gene)
    else:
        raise ValueError("Outgroup not in tree")

subprocess.call(cmd, shell=True)