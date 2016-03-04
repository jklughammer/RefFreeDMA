#! /usr/bin/env python

__author__ = 'JKlughammer'

import os
import math
import sys



infile_name = str(sys.argv[1])
tolerance = int(sys.argv[2]) #can be different from max_read_length but sofar set to the same value and therefore doesn't do anything
max_read_length = int(sys.argv[3])
working_dir=str(sys.argv[4])
cLimit=float(sys.argv[5])


file_name = infile_name.split(".")[0]
infile = open(infile_name, 'r')
refFasta = open(file_name + "_" + str(tolerance) + "_collapse.fa", 'w')
refFastq = open(file_name + "_" + str(tolerance) + "_collapse.fq", 'w')
outRef = open(file_name + "_" + str(tolerance) + "_collapse.ref", 'w')



count = 0
name = []
orig = []
conv = []
count = 0
start = 1
include = 1
last_include = 0

#reads need to be sorted by converted sequence and from smalles to largest

print("Now collapsing reads...")

for line in infile:
	count = count + 1
#reastart matching after encountering the last non TGG read (actually not needed, because non TGG reads are filtered out before (prepareReads.sh))
#	if (line.split("\t")[2].strip()[0:3] != "TGG") & (last_include == 0):
#		print("no TGG")
#		name = [line.split("\t")[0].strip()]
#		orig = [line.split("\t")[1].strip()]
#		conv = [line.split("\t")[2].strip()]
#		include = 0
#		continue
	if count == 1:
		name = [line.split("\t")[0].strip()]
		orig = [line.split("\t")[1].strip()]
		conv = [line.split("\t")[2].strip()]

#collect all reads that are substrings of the next (longer) read in the converted form.
#Only require the first [tolerance] bases to be a substring of the next read 		
	if conv[len(conv) - 1][0:tolerance] in line.split("\t")[
		2].strip():  
		name.append(line.split("\t")[0].strip())
		orig.append(line.split("\t")[1].strip())
		conv.append(line.split("\t")[2].strip())
		last_include = 1
		continue
	end = count - 1
	cMerge = []

#built consensus for group of matching reads on the basis of the real read sequence
	cCount = dict()
	tCount = dict()
	tChange = dict()
	cChange = dict()
	for read in orig:
		for i in range(0, len(read)):
			if read[i] == "C":
				if i in cCount.keys():
					cCount[i] = cCount[i] + 1
				else:
					cCount[i] = 1
					tCount[i] = 0
			if read[i] == "T":
				if i in tCount.keys():
					tCount[i] = tCount[i] + 1
				else:
					tCount[i] = 1

	for key in cCount.keys():
		quot=float(cCount[key])/(float(cCount[key])+float(tCount[key]))
		if quot < cLimit:
			tChange[key] = quot
		else:
			cChange[key] = quot

	ref = list(max(orig, key=len))
	for idx in cChange.keys():
		ref[idx] = "C"
	for idx in tChange.keys():
		ref[idx] = "T"

	ref = "".join(ref)

#output consensus read if it is a TGG read	
	fastaName = "dedRef_" + str(start) + "-" + str(end) + "_" + str(len(cChange.keys()))
	if include <> 0:
		###for fasta
		refFasta.writelines(">" + fastaName + "\n")
		refFasta.writelines(max(conv, key=len) + "\n")
		###for fastq
		refFastq.writelines("@" + fastaName + "\n")
		refFastq.writelines(max(conv, key=len) + "\n")
		refFastq.writelines("+\n")
		refFastq.writelines("".join(["d"] * (len(ref))) + "\n")
		### ref for concatinated fasta
		outRef.writelines(fastaName + "\t" + ref + "\t" + max(conv, key=len) + "\t" + str(end - start) + "\n")

	if count % 1000000 == 0:
		print(str(count) + " done.")

#restart matching for new group of reads 
	start = count
	name = [line.split("\t")[0].strip()]
	orig = [line.split("\t")[1].strip()]
	conv = [line.split("\t")[2].strip()]
	include = 1
	last_include = 0
infile.close()
refFastq.close()
open(os.path.join(working_dir,"makePreConsensus.done"), 'a')
