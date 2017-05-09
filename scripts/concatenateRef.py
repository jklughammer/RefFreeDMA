#! /usr/bin/env python
__author__ = 'jklughammer'

import math
import time
import sys
import os

infile_name = str(sys.argv[1])
file_name = infile_name.split("_norc")[0]


infile = open(infile_name, 'r')
outfile = open(file_name+"_concat.fa", 'w')
temp = open(file_name+"_temp", 'w')
annot = open(file_name+"_concat.bed", 'w')
fasta = open(file_name+".fa", 'w')
fastq = open(file_name+".fq", 'w')

maxReadLen = int(sys.argv[2])
buffer = "".join(["N"]*(maxReadLen-1))+"C"
working_dir=str(sys.argv[3])
count = 0

for line in infile:
    ref = line.split("\t")[1].strip()[0:maxReadLen]# only use at max the first 48 bases
    if count == 0:
        start = maxReadLen

    else:
        start = maxReadLen + end

    end = start + len(ref)

    temp.writelines(buffer + ref)

    annot.writelines(file_name.split("/")[len(file_name.split("/"))-1]+"_concat"+"\t" + str(start) + "\t" + str(end) +
                     "\t" + line.split("\t")[0]+ "\t" + line.split("\t")[0].split("_")[2] + "\t" + line.split("\t")[2] + "\t" + line.split("\t")[3]
                     + "\t" + line.split("\t")[4] + "\t" + line.split("\t")[5].strip() + "\n")
    ##write fastq
    fastq.writelines("@"+line.split("\t")[0]+"\n")
    fastq.writelines(ref+"\n")
    fastq.writelines("+\n")
    fastq.writelines("".join(["I"]*(len(ref)))+"\n")
    ##write fasta
    fasta.writelines(">"+line.split("\t")[0]+"\n")
    fasta.writelines(ref+"\n")


    count = count + 1
    if count % 200000 == 0:
        print(count)

temp.close()

print("now writing fasta file")

temp = open(file_name+"_temp", 'r')


out = "".join(temp.readline())
outfile.writelines(">"+file_name.split("/")[len(file_name.split("/"))-1]+"_concat"+"\n")
for i in range(0,int(math.ceil(len(out)/50))+1):
    outfile.writelines(out[i*50:i*50+50] + "\n")
    if i % 200000 == 0:
        print(i)



open(os.path.join(working_dir,"concatenateRef.done"), 'a')
