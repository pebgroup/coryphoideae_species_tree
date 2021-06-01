#!/usr/bin/python3

import argparse, dendropy, subprocess

parser = argparse.ArgumentParser()
parser.add_argument("treefile")
args = parser.parse_args()
treefile = str(args.treefile)

tree = dendropy.Tree.get(path=treefile, schema="newick")
tips = str(tree.taxon_namespace)

if "1079" in tips and "1080" in tips and "1081" in tips and "1082" in tips:
    cmd = "pxrr -t "+treefile+" -g 1079,1080,1081,1082 -o temp.tre"
else:
    if  "1079" in tips and "1080" in tips and "1081" in tips: 
        cmd = "pxrr -t "+treefile+" -g 1079,1080,1081 -o temp.tre"
    elif "1079" in tips and "1080" in tips and "1082" in tips:
        cmd = "pxrr -t "+treefile+" -g 1079,1080,1082 -o temp.tre"
    elif "1079" in tips and "1081" in tips and "1082" in tips:
        cmd = "pxrr -t "+treefile+" -g 1079,1081,1082 -o temp.tre"
    elif "1080" in tips and "1081" in tips and "1082" in tips:
        cmd = "pxrr -t "+treefile+" -g 1080,1081,1082 -o temp.tre"
    elif "1079" in tips and "1080" in tips:
        cmd = "pxrr -t "+treefile+" -g 1079,1080 -o temp.tre"
    elif "1079" in tips and "1081" in tips:
        cmd = "pxrr -t "+treefile+" -g 1079,1081 -o temp.tre"
    elif "1079" in tips and "1082" in tips:
        cmd = "pxrr -t "+treefile+" -g 1079,1082 -o temp.tre"
    elif "1080" in tips and "1081" in tips:
        cmd = "pxrr -t "+treefile+" -g 1080,1081 -o temp.tre"
    elif "1080" in tips and "1082" in tips:
        cmd = "pxrr -t "+treefile+" -g 1080,1082 -o temp.tre"
    elif "1081" in tips and "1082" in tips:
        cmd = "pxrr -t "+treefile+" -g 1081,1082 -o temp.tre"
    elif "1079" in tips:
        cmd = "pxrr -t "+treefile+" -g 1079 -o temp.tre"
    elif "1080" in tips:
        cmd = "pxrr -t "+treefile+" -g 1080 -o temp.tre"
    elif "1081" in tips:
        cmd = "pxrr -t "+treefile+" -g 1081 -o temp.tre"
    elif "1082" in tips:
        cmd = "pxrr -t "+treefile+" -g 1082 -o temp.tre"
    else:
        raise ValueError("Outgroup not in tree")

subprocess.call(cmd, shell=True)