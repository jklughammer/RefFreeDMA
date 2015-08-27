#! /usr/bin/env python
__author__ = 'jklughammer'

import sys
import collections
import re
import os
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

infile_name = str(sys.argv[1])  # needs ordered input (by sequence length)
working_dir=str(sys.argv[2])

infile = open(infile_name, 'r')
outfile = open(infile_name + "-res", 'w')


# make 2 letter lists for RC and FW
print("making search dicts...")
RC = collections.OrderedDict()
RC_rev = collections.OrderedDict()
FW = collections.OrderedDict()
count = 0
for line in infile:
	seq = line.split("\t")[1]
	name = line.split("\t")[0]

	#RC
	rc_seq = str(Seq(seq, generic_dna).reverse_complement())
	rc_seq = rc_seq.replace("C", "T")
	rc_seq = rc_seq.replace("G", "A")
	RC[name] = [rc_seq, seq]
	if not rc_seq in RC_rev:
		RC_rev[rc_seq] = dict()
		RC_rev[rc_seq][name] = seq
	else:
		RC_rev[rc_seq][name] = seq
		#print("key " + str(rc_seq) + " already exists: " + str(RC_rev[rc_seq]) + ". Adding " + str(name))

	#FW
	fw_seq = seq.replace("C", "T")
	fw_seq = fw_seq.replace("G", "A")
	FW[name] = [fw_seq, seq]
	count = count + 1

print("now finding matches...")
origEntries = len(FW)
count = 0
hitcount = 0
already_popped = 0
keep_correct=dict()



for fw_key, fw_both in FW.iteritems():
	fw_val=fw_both[0]
	fw_orig=fw_both[1]
	count += 1
	if count % 10000 == 0:
		print(count)
	if fw_val in RC_rev and (not fw_key in RC_rev[fw_val] or len(RC_rev[fw_val])>1):
		hitcount += 1
		fw_orig_l=list(fw_orig)
		temp=RC_rev[fw_val]
		if len(RC_rev[fw_val])>1:
			similarity = dict()
			for RC_rev_key, RC_rev_val in RC_rev[fw_val].iteritems():
				if RC_rev_key == fw_key:
					continue
				similarity[RC_rev_key]= 0
				RC_rev_val_l = list(str(Seq(RC_rev_val, generic_dna).reverse_complement()))
				for i in range(0,len(fw_orig_l)-1):
					if (fw_orig_l[i] == "T" and rw_orig_l[i] in ["T","C"]) or (fw_orig_l[i] == "C" and rw_orig_l[i] in ["C"])\
							or (fw_orig_l[i] == "G" and rw_orig_l[i] in ["G","A"])or (fw_orig_l[i] == "A" and rw_orig_l[i] in ["A"]):
						similarity[RC_rev_key]+= 1
			rw_key = max(similarity)
		else:
			rw_key = RC_rev[fw_val].keys()[0]

		rw_orig = RC_rev[fw_val][rw_key]
		rw_orig_l=list(str(Seq(rw_orig, generic_dna).reverse_complement()))
		new_fw_l=list(fw_orig)
		for i in range(0,len(fw_orig_l)-1):
			if fw_orig_l[i] == "T" and rw_orig_l[i] == "C":
				new_fw_l[i] = "C"
		new_fw = "".join(new_fw_l)
		keep_correct[fw_key] = new_fw
		try: FW.pop(rw_key) #prevent also finding the inverse combination
		except: already_popped+=1

print(str(hitcount) + " RC-combies found")
print(str(already_popped) + " multihits found")
print(str(len(FW)) + " from " + str(origEntries) + " left")



print("now printing RC-reduced output")
printCount = 0
infile.seek(0)
for line in infile:
	spl_line=line.split("\t")
	name = spl_line[0]
	if name in FW:
		printCount += 1
		if name in keep_correct:
			spl_line[1] = keep_correct[name]
			join_line = "\t".join(spl_line)
			outfile.writelines(join_line)
			keep_correct.pop(name)
		else:
			outfile.writelines(line)
			FW.pop(name)
	if printCount % 10000 == 0:
		print(printCount)
print(str(printCount) + " lines printed")

open(os.path.join(working_dir,"mergeRevComp.done"), 'a')