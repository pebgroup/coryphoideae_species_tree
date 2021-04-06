#!/usr/bin/python3

import os, argparse, subprocess
from Bio import SeqIO
from Bio.Seq import Seq

parser = argparse.ArgumentParser()
parser.add_argument("sample")
args = parser.parse_args()
sample = str(args.sample)

# depth required to KEEP (i.e. anything <trshld will be discarded)
trshld = 2

# Get all subdirectories in the current working directory. these are the loci recovered by hybpiper
loci = next(os.walk(sample))[1]
sequences = {}

for locus in loci: 

	#exon
	#pth = sample+'/'+locus+'/'+sample+'/sequences/FNA/'+locus+'.FNA'	

	#supercontig
	pth = sample+'/'+locus+'/'+sample+'/sequences/intron/'+locus+'_supercontig.fasta'	
	
	if os.path.isfile(pth):
		for record in SeqIO.parse(pth, "fasta"):
			record.id = record.id+'_'+locus
			#sequences.append(record)
			sequences[record.id] = record
			
with open('../5_coverage/'+sample+'.fasta', "w") as outfile:
 	SeqIO.write(list(sequences.values()), outfile, "fasta")

print(sample+'.fasta generated')
	
# BWA index targets
cmd = 'bwa index ../5_coverage/'+sample+'.fasta'
subprocess.call(cmd,shell=True)
print(sample+'.fasta indexed')

# BWA mem paired reads
cmd = 'bwa mem ../5_coverage/'+sample+'.fasta ../2_trimmed/'+sample+'_1P.fastq ../2_trimmed/'+sample+'_2P.fastq | samtools view -b -o ../5_coverage/'+sample+'.bam'
subprocess.call(cmd,shell=True)
print('paired reads mapped to '+sample+'.fasta')

# BWA mem unpaired reads
cmd = 'bwa mem ../5_coverage/'+sample+'.fasta ../2_trimmed/'+sample+'_UN.fastq | samtools view -b -o ../5_coverage/'+sample+'_up.bam'
subprocess.call(cmd,shell=True)
print('unpaired reads mapped to '+sample+'.fasta')

# merge BAM files
cmd = 'samtools merge ../5_coverage/'+sample+'_all.bam ../5_coverage/'+sample+'.bam ../5_coverage/'+sample+'_up.bam'
subprocess.call(cmd,shell=True)
print('BAMs merged')

# sort and index BAM files
cmd = 'samtools sort ../5_coverage/'+sample+'_all.bam -o ../5_coverage/'+sample+'_all_sorted.bam'
subprocess.call(cmd,shell=True)
cmd = 'samtools index ../5_coverage/'+sample+'_all_sorted.bam'
subprocess.call(cmd,shell=True)
print('BAM indexed and sorted')

# remove duplicates
cmd = 'java -jar ~/software/picard.jar MarkDuplicates I=../5_coverage/'+sample+'_all_sorted.bam O=../5_coverage/'+sample+'_all_sorted_deduplicated.bam M=../5_coverage/'+sample+'marked_dup_metrics.txt REMOVE_DUPLICATES=true'
subprocess.call(cmd,shell=True)
print('reads deduplicated for sample '+sample)

# calculate coverage
cmd = 'samtools depth ../5_coverage/'+sample+'_all_sorted_deduplicated.bam > ../5_coverage/'+sample+'.cov'
subprocess.call(cmd,shell=True)
print('coverage calculated for sample '+sample)

# define function to replace nth position of sequence with N
def n2N(sqnc, pstn):
	sqnc = list(sqnc)
	sqnc[int(pstn)-1] = "N"
	return "".join(sqnc)

# process coverage
with open('../5_coverage/'+sample+'.cov', "r") as covfile:
	for line in covfile:
		line = line.strip()
		LINE = line.split("\t")
		if int(LINE[2]) < trshld:
			sequences[LINE[0]].seq = n2N(sequences[LINE[0]].seq,LINE[1])

# remove unnecessary leading and trailing Ns
for nm in sequences.keys():
	sequences[nm].seq = sequences[nm].seq.strip("N")
	if isinstance(sequences[nm].seq, str):
		sequences[nm].seq = Seq(sequences[nm].seq)

print('coverage trimming completed, keeping only positions with coverage of '+str(trshld)+' or above')

# write outfile
with open('../5_coverage/'+sample+'_trimmed.fasta', "w") as outfile:
	SeqIO.write(list(sequences.values()), outfile, "fasta")
print('trimmed seqs written to '+sample+'_trimmed.fasta')

# remove unnecessary files
#cmd = "find ../coverage -type f ! -name '*.fasta' -delete"
#subprocess.call(cmd, shell=True)
#print('tidied up.')