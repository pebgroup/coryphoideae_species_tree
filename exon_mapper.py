#!/usr/bin/python3

from dataclasses import replace
import os, argparse, subprocess
import numpy as np
import pandas as pd
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("--outdir")
parser.add_argument("--gene")
parser.add_argument("--file_ending")
args = parser.parse_args()
outdir = args.outdir
gene = args.gene
file_ending = args.file_ending

#Creating output file ending
output_file_ending = file_ending.replace(".fasta", "_mapped.fasta" )


# import recovery statistics
df = pd.read_csv("../../03_hybpiper/seq_lengths.txt", sep="\t")
df = df.drop(0)


# for each aligned gene, identify the two samples with the highest recovery stats
#gene = fn.split("_")[1]
if gene in df.columns:
	df_red = df[['Species',gene]].sort_values(gene, ascending=False)
	df_red = df_red.reset_index()

	# Find the correct gene in the fasta file and save it as exon1
	sp = df_red["Species"][0]
	with open(gene+"_aligned.fasta", "r") as fasta_file:
		for record in SeqIO.parse(fasta_file, "fasta"):
			if record.id.startswith("exon1"):
				exon1 = record
				break
	exon1.id = "exon1"

	# second sample
	sp = df_red["Species"][1]
	with open(gene+"_aligned.fasta", "r") as fasta_file:
		for record in SeqIO.parse(fasta_file, "fasta"):
			if record.id.startswith("exon2"):
				exon2 = record
				break
	exon2.id = "exon2"

	with open(gene+"_temp.fasta", "w") as output_handle:
		SeqIO.write([exon1, exon2], output_handle, "fasta")
	subprocess.call("mafft --keeplength --add "+gene+"_temp.fasta "+gene+file_ending+" > "+outdir+"/"+gene+output_file_ending,shell=True)
	subprocess.call("rm "+gene+"_temp.fasta", shell=True)
else:
	print(gene+" not there!!!!")