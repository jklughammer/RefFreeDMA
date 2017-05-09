#! /usr/bin/env python

#Authors: Andreas Schoenegger, Johanna Klughammer, Nathan Sheffield, Christoph Bock

# Summary: iterates through alignment file sam or bam (see samtools parameter) and determines the number of concordant and discordant reads at every CpG
# writes bed file with the following columns: "chr","start", "strand", "MethylatedReadCount", "UnmethylatedReadCount", "StrangeReadCount",
#	                        "ConcordantMethylatedReadCount", "ConcordantUnmethylatedReadCount", "DiscordantReadCount",
#	                        "DiscordantReadCount", "MethylPlus", "MethylMinus", "UnmethylPlus", "UnmethylMinus"
#


import os
import math
import re
from Bio import SeqIO
import sys
import pysam


# definition of global variables
chromFiles = {}  # will be initialized on the first call to lookupSequence()

global genome
global debug
debug = True


# removes newline characters at the end of a string (similar to rstrip() but leaves tabs intact)
def removeLinefeed(s):
	while s[-1] in ["\n", "\r"]: s = s[0:-1]
	return s

# obtain the DNA sequence for a given genomic region
#JK: New function lookupSequence()
#Major changes:
#1. load chromosome sequences from an indexed genome fasta
#2. check CIGAR for insertions or deletions: if present use this information to retrieve the correct genomic sequence
#if not retrieve genomic sequence as was before.
def lookupSequence(genome, chrom, chromstart, chromend, cigar, lowercaseRepeats=False):
	global chromFiles

	# load the corresponding chromosome file
	try:
		chromFile = chromFiles[chrom]
	except KeyError:
		print "Reading chromosome file: " + chrom
		chromFiles[chrom] = genome[chrom].seq
		if debug: print "Reading chromosome file completed"

	if "I" in cigar or "D" in cigar:
		deletion = 0
		insertion = 0
		pos = 0
		cigar_str = re.findall('([MID])', cigar)  #JK: array containing the CIGAR IDs
		cigar_int = re.findall('(\d+)', cigar)  #JK: corresponding array containing the counts

		#JK: first calculate the new chromend by including insertions and deletions
		for i in range(len(cigar_str)):
			if cigar_str[i] == "D":
				deletion = deletion + int(cigar_int[i])
			if cigar_str[i] == "I":
				insertion = insertion + int(cigar_int[i])
		chromend = chromend + deletion - insertion
		#JK: get sequence using the updated chromend
		seq = list(str(chromFiles[chrom][chromstart:chromend]))
		#JK: perform insertions (as "N") and deletions in the retrieved genomic sequence according to the CIGAR info
		for i in range(len(cigar_str)):
			if cigar_str[i] == "D":
				for j in range(int(cigar_int[i])):
					seq.pop(pos)
			elif cigar_str[i] == "I":
				for j in range(int(cigar_int[i])):
					seq.insert(pos, "N")
					pos = pos + int(cigar_int[i])
			elif cigar_str[i] == "M":
				pos = pos + int(cigar_int[i])
		seq = "".join(seq)
	else:
		#retrieve genomic sequence with no modifications
		seq = str(chromFiles[chrom][chromstart:chromend])

	if not lowercaseRepeats: seq = seq.upper()
	if debug: print "Sequence of %s:%i-%i (%s) is: %s" % (chrom, chromstart, chromend, genome, seq)
	return seq


def passFilter(line):

	# parsing the current line of the SAM file
	cells = line.split("\t")
	if len(cells) < 11:
		raise Exception("WARNING: Skipping a line that is not consistent with BAM format: " + line)
	sam_flag = cells[1]
	chrom = cells[2]
	chromstart = int(cells[3]) - 1  # correct for the 0-based vs. 1-based difference between BAM and lookupSequence()
	cigar = cells[5]
	bisSeq = cells[9]
	chromend = chromstart + len(bisSeq)

	#JK:set passfilter to false for unmapped reads and reads that have strange CIGAR strings
	if int(sam_flag) & 4:
		print "JK: Exclude unmapped reads"
		return False

	if len(re.findall('[SHPX=]', cigar)) != 0:
		print "JK: skip if any other possible CIGAR characters except for MID (seldom and therefore not handled)"
		return False

	if "N" in cigar:
		return False

	return True  # keep everything that made it until here.


def checkConcordance(line, genome):
	# parsing the current line of the SAM file
	cells = line.split("\t")
	if len(cells) < 11:
		raise Exception("WARNING: Skipping a line that is not consistent with BAM format: " + line)

	sam_flag = cells[1]
	chrom = cells[2]
	chromstart = int(cells[3]) - 1  # fix the 0-based vs. 1-based difference between BAM and lookupSequence()
	cigar = cells[5]
	bisSeq = cells[9]

	#Introduced by JK 07.02.2017 to remove fill-in cs
	if int(sam_flag) & 16:	
		if bisSeq[:2] == "CA":
			bisSeq = bisSeq[2:]
			chromstart = chromstart + 2
	else:
		if bisSeq[-2:] == "TG":
			bisSeq = bisSeq[:-2]
		
				
	chromend = chromstart + len(bisSeq)


	if (int(sam_flag) & 16 and not (int(sam_flag) & 128)) or (not (int(sam_flag) & 16) and int(sam_flag) & 128):
		g2a = True
	else:
		g2a = False

	# obtaining the reference genome sequence
	try:
		refSeq = lookupSequence(genome, chrom, chromstart, chromend, cigar)

	except (Exception) as ex:
		raise Exception('Could not retrieve sequence for current line ("' + line[0:100] + '") due to the following exception: ' + str(type(ex)) + ": " + str(ex))

	# testing whether the current read shows sufficient evidence of correct bisulfite conversion
	#return [bisSeq, refSeq, chrom, chromstart, chromend, g2a]

	return concordanceAnalysis(chrom, chromstart, refSeq, bisSeq, sam_flag, g2a) # [true/false]


def concordanceAnalysis(chrom, chromstart, refSeq, bisSeq, sam_flag, g2a):
	#AS:
	concordant_read = False

	methylated_sites = 0
	unmethylated_sites = 0
	strange_sites = 0

	minimum_cgs_per_read = 4  # TODO: should be a parameter to the script (changed from 3 to 4 by JK 07.02.2017)
	actual_num_cgs = refSeq.count("CG")

	if actual_num_cgs >= minimum_cgs_per_read:

		cg_sites = list()

		# analyze the whole read for concordance/discordance
		for i in range(0, len(bisSeq)):	

			if refSeq[i:i + 2] == "CG":
				# deal with the strand of the read (related to g2a boolean). TODO: potential bug for WGBS PE
				strand = "+"
				if int(sam_flag) & 16:
					strand = "-"

				cpg_chrom_pos = chromstart + i

				# check what base the read has at the reference-genome CG position
				methylated = 0
				unmethylated = 0
				strange = 0

				methylated_plus = 0
				unmethylated_plus = 0

				methylated_minus = 0
				unmethylated_minus = 0

				# read mapped to forward strand
				if not g2a:
					if bisSeq[i] == "C":    # methylated
						methylated = 1
						methylated_sites += 1

						if strand == "+":
							methylated_plus = 1
						else:
							methylated_minus = 1

					elif bisSeq[i] == "T":    # unmethylated
						unmethylated = 1
						unmethylated_sites += 1

						if strand == "+":
							unmethylated_plus = 1
						else:
							unmethylated_minus = 1

					else:
						strange = 1
						strange_sites += 1

				# read mapped to reverse strand
				else:
					if i == len(bisSeq)-1:
						strange = 1
						strange_sites += 1
						break
					if bisSeq[i+1] == "G":    # methylated
						methylated = 1
						methylated_sites += 1

						if strand == "+":
							methylated_plus = 1
						else:
							methylated_minus = 1

					elif bisSeq[i+1] == "A":    # unmethylated
						unmethylated = 1
						unmethylated_sites += 1

						if strand == "+":
							unmethylated_plus = 1
						else:
							unmethylated_minus = 1

					else:
						strange = 1
						strange_sites += 1

				cg_sites.append([chrom, cpg_chrom_pos, methylated, unmethylated, strange, methylated_plus, methylated_minus, unmethylated_plus, unmethylated_minus])    # len = 9

		if (methylated_sites + unmethylated_sites + strange_sites) != actual_num_cgs:
			raise Exception("meth_cg + unmeth_cg + strange_cg do not add up. This should not happen! refSeq: " + refSeq + ", bisSeq: " + bisSeq)
			#TODO: sanity check especially because of snps

		else:

			# Assess whether the CpGs on this read are concordant or not:
			num_cgs = actual_num_cgs
			if strange_sites > 0:
				num_cgs -= strange_sites    # just pretend for now, that snps do not affect the concordant state. deal with snp problem later, maybe filter out the whole read?

			#concordant_read = False		# this is set to false on top of function
			if (methylated_sites or unmethylated_sites) == actual_num_cgs:	#TODO: this checks for the CG-sites in the reference. As soon as there is a SNP/mismatch the concordance fails
				concordant_read = True

			# Let's revisit the CpGs and write them into the global CpG dictionary:

			# add chrom to global dict if chromosome doesn't exist yet
			if chrom not in cpg_dict:
				cpg_dict[chrom] = dict()

			# cg_sites remembers all relevant information about the cpg sites discovered above:
			for one_cpg in cg_sites:
				cpg_chrom_pos = one_cpg[1]
				methylated = one_cpg[2]
				unmethylated = one_cpg[3]
				strange = one_cpg[4]
				methylated_plus = one_cpg[5]
				methylated_minus = one_cpg[6]

				unmethylated_plus = one_cpg[7]
				unmethylated_minus = one_cpg[8]

				# variable assignment very inefficient but code more readable

				# add chrom_pos to this chromosome if it doesn't exist
				if cpg_chrom_pos not in cpg_dict[chrom]:
					cpg_dict[chrom][cpg_chrom_pos] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # len = 11
					# 0: methylated readcount
					# 1: unmethylated readcount
					# 2: strange readcount (snp?)
					# 3: methylated concordant
					# 4: unmethylated concordant
					# 5: concordant
					# 6: discordant
					# 7: methylated_plus
					# 8: methylated_minus
					# 9: unmethylated_plus
					# 10: unmethylated_minus

				# increase meth/unmeth/strange counters (based on values in first for loop
				cpg_dict[chrom][cpg_chrom_pos][0] += methylated
				cpg_dict[chrom][cpg_chrom_pos][1] += unmethylated
				cpg_dict[chrom][cpg_chrom_pos][2] += strange

				if concordant_read:
						cpg_dict[chrom][cpg_chrom_pos][5] += 1    # +1 if seen in concordant read regardless whether methylated or not

						if methylated == 1: cpg_dict[chrom][cpg_chrom_pos][3] += 1
						if unmethylated == 1: cpg_dict[chrom][cpg_chrom_pos][4] += 1

				else:
						cpg_dict[chrom][cpg_chrom_pos][6] += 1    # +1 if part of a discordant read

				# add counts for plus and minus strand:
				cpg_dict[chrom][cpg_chrom_pos][7] += methylated_plus
				cpg_dict[chrom][cpg_chrom_pos][8] += methylated_minus
				cpg_dict[chrom][cpg_chrom_pos][9] += unmethylated_plus
				cpg_dict[chrom][cpg_chrom_pos][10] += unmethylated_minus

	return concordant_read



# main analysis procedure
def performAnalysis(options):
	# prepare output file
	outfile = open(options.outfile, 'w')
	if options.produce_sam:
		sam_concordance_outfile = open(options.concordantOutfile, "w")
		sam_discordance_outfile = open(options.discordantOutfile, "w")

	if options.skipped != "": skipped = open(options.skipped, 'w')

	# download genome sequence (if necessary)
	if not os.path.exists(options.genomeDir + os.sep + options.genome + os.sep + options.genome + ".fa"):
		print("Genome assembly not found for '" + options.genome )
		raise SystemExit

	# JK Initialize genome database (using SeqIO.index_db makes it possible to use a single fasta file instead of having all chroms in separate fasta files)
	idx_file=os.path.join(options.genomeDir, options.genome, options.genome + ".idx")
	tmp_idx_file=os.path.join(options.genomeDir, options.genome, options.genome + ".idx_tmp")	
	fasta=os.path.join(options.genomeDir, options.genome, options.genome + ".fa")

	if not os.path.exists(idx_file):	
		genome_temp = SeqIO.index_db(tmp_idx_file,fasta,"fasta")
		os.system("mv " + tmp_idx_file + " " + idx_file)
	
	genome = SeqIO.index_db(idx_file,fasta,	"fasta")
	print("Genome " + options.genome + " indexed")

	

	#JK convert bam to samfile if samtools is set, otherwise assume input is already sam

	if not options.samtools == "":
		in_sam=	re.sub(".bam$",".sam",options.infile)
		cmd = options.samtools + "/samtools view " + options.infile +">"+ in_sam
		print(cmd)	
		os.system( cmd )	
		infile = open(in_sam, 'r')	
	else:
		infile = open(options.infile, 'r')

	# iterate through all lines of the SAM file

	count = 0
	warnings = 0
	acceptedRecords = 0
	unmappedRecords = 0
	skippedRecords = 0
	splicedRecords = 0
	cigarskippedRecords = 0
	unmatchedReads_accepted = 0
	unmatchedReads_skipped = 0
	isMate1 = True
	prevLine = ""
	print("Starting the read filtering")
	if not options.pairedEnd:
		print "single end mode"
	else:
		raise NotImplementedError("Paired-end mode not yet implemented for RRBS")
	
	for line in infile:
	

		if count % 1000000 == 0: sys.stdout.write(".")
		if count % 100000000 == 0: sys.stdout.write("\n.")


		if warnings > 25:
			print("Maximum number of warning messages reached. Aborting.")
			break

		count +=1
		# optionally process header lines
		if count <= options.headerLines:
			if debug: print("Header line: " + line)
			outfile.write(line + "\n")
			if options.skipped != "": skipped.write(line + "\n")  #JK: write header lines in out AND skipped file
			continue
		
		line = removeLinefeed(line)
		if debug: print(line)

		line_spl = line.split("\t")
		chrom = line_spl[2]

		if chrom in [options.genome]:

			# testing whether the current read or read pair shows sufficient evidence of correct bisulfite conversion
			try:
				# AS: pairedEnd completely ignored for RRBS (therefore the code below is out of date), continue reading at else!
				if options.pairedEnd:
						if isMate1:
							isMate1 = False
							prevLine = line
							continue
						else:
							if debug: print("Processing paired-end read")
							if prevLine.split("\t")[0] != line_spl[0]:
								# process unmatched read as a singleton and proceed
								if passFilter(prevLine, genome):
									outfile.write(prevLine + "\n")
									acceptedRecords += 1
									unmatchedReads_accepted += 1
								else:
									if options.skipped != "": skipped.write(prevLine + "\n")
									skippedRecords += 1
									unmatchedReads_skipped += 1
								isMate1 = False
								prevLine = line
							else:
								if passFilter(prevLine, genome) and passFilter(line, genome):
									outfile.write(prevLine + "\n" + line + "\n")
									acceptedRecords += 2
								else:
									if options.skipped != "": skipped.write(prevLine + "\n" + line + "\n")
									skippedRecords += 2
								isMate1 = True
								prevLine = ""

				# single-end mode:
				else:
					if passFilter(line):
						acceptedRecords += 1

						if checkConcordance(line, genome):
							if options.produce_sam:
								sam_concordance_outfile.write(line + "\n")
						else:
							if options.produce_sam:
								sam_discordance_outfile.write(line + "\n")

						if acceptedRecords % 1000000 == 0:
							print acceptedRecords


					else:
						if options.skipped != "": skipped.write(line + "\n")
						skippedRecords += 1

			except Exception as ex:
				print "Exception!"
				print ex
				warnings += 1
				raise

	print "Writing bed file"
	write_cpg_info(outfile)

	# print summary statistics
	print("Total number of processed lines:   " + str(count))
	print("Total number of accepted records: " + str(acceptedRecords))
	if options.pairedEnd:
		print("Total number of accepted unmatched reads: " + str(unmatchedReads_accepted))

	# cleaning up
	if options.produce_sam:
		sam_concordance_outfile.close()
		sam_discordance_outfile.close()

	infile.close()
	outfile.close()

	if options.skipped != "": skipped.close()
	os.remove(in_sam)


def write_cpg_info(outfile):

	# print header line for output file:
	header = "\t".join(["chr","start", "strand", "MethylatedReadCount", "UnmethylatedReadCount", "StrangeReadCount",
							"ConcordantMethylatedReadCount", "ConcordantUnmethylatedReadCount",
							"DiscordantReadCount", "MethylPlus", "MethylMinus", "UnmethylPlus", "UnmethylMinus"])
	outfile.write(header + "\n")
	outfile.flush()

	#reminder for list content:
	#
	#cpg_dict[chrom][cpg_chrom_pos] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
	# 0: methylated readcount
	# 1: unmethylated readcount
	# 2: strange readcount (snp?)
	# 3: methylated concordant
	# 4: unmethylated concordant
	# 5: concordant
	# 6: discordant
	# 7: methylated_plus
	# 8: methylated_minus
	# 9: unmethylated_plus
	# 10: unmethylated_minus

	# TODO: proper sorting of chromosome names
	for chrom in sorted(cpg_dict.keys()):
		for chromstart in sorted(cpg_dict[chrom].keys()):
			one_cpg = cpg_dict[chrom][chromstart]

			MethylatedReadCount = str(one_cpg[0])
			UnmethylatedReadCount = str(one_cpg[1])
			StrangeReadCount = str(one_cpg[2])
			ConcordantMethylatedReadCount = str(one_cpg[3])
			ConcordantUnmethylatedReadCount = str(one_cpg[4])
			DiscordantReadCount = str(one_cpg[6])

			MethylPlus = str(one_cpg[7])
			MethylMinus = str(one_cpg[8])
			UnmethylPlus = str(one_cpg[9])
			UnmethylMinus = str(one_cpg[10])

			strand = "-"
			# Convention: if there is at least one read mapped on the forward strand -> plus strand
			if one_cpg[7] > 0 or one_cpg[9]:
				strand = "+"

			line = "\t".join([chrom, str(chromstart), strand, MethylatedReadCount, UnmethylatedReadCount, StrangeReadCount, ConcordantMethylatedReadCount, ConcordantUnmethylatedReadCount, DiscordantReadCount, MethylPlus, MethylMinus, UnmethylPlus, UnmethylMinus])
			#print line
			outfile.write(line + "\n")


# AS:
cpg_dict = dict()

if __name__ == '__main__':
	print "Starting program..."
	# constructing command line parser
	import optparse

	parser = optparse.OptionParser()
	parser.add_option('--samtools', action='store', type='string', dest='samtools', help='Specify the path to samtools. If set expects bam as input; if not sam', default="")
	parser.add_option('--infile', action='store', type='string', dest='infile', help='Specify the name of the input file', default="")
	parser.add_option('--outfile', action='store', type='string', dest='outfile', help='Specify the name of the output file', default="")

	parser.add_option('--concordantOutfile', action='store', type='string', dest='concordantOutfile', help='Specify the name of the concordant output file', default="")
	parser.add_option('--discordantOutfile', action='store', type='string', dest='discordantOutfile', help='Specify the name of the discordant output file', default="")

	parser.add_option('--skipped', action='store', type='string', dest='skipped', help='Specify the name of the file to store skipped reads (optional)', default="")
	parser.add_option('--skipHeaderLines', action='store', type='int', dest='headerLines', help='Specify the number of header lines that should be skipped (default: 0)', default=0)
	parser.add_option('--pairedEnd', action='store_true', dest='pairedEnd', help='Specify paired-end mode (default: False)', default=False)
	parser.add_option('--produce_sam', action='store_true', dest='produce_sam', help='Specify whether samfiles for concordant and disconcordant reads should be produced (default: False)', default=False)
	parser.add_option('--genome', action='store', type='string', dest='genome', help='Specify the genome for which the DNA sequence should be retrieved')
	parser.add_option('--genomeDir', action='store', type='string', dest='genomeDir', help='Specify the name of the directory containing the genome sequence data (will be downloaded from USCS Genome Browser if not available locally)',
					  default=".")
	parser.add_option('-v', '--verbose', action='store_true', dest='verbose', help='Print debugging information (default: False)', default=False)

	(options, args) = parser.parse_args()
	debug = options.verbose
	if options.outfile == "": options.outfile = options.infile + ".output.txt"
	if not options.infile or not options.genome:
		print("Mandatory parameters missing. Program will terminate now.")
		print("\nYour parameter settings:")
		print(options)
		raise SystemExit

	performAnalysis(options)

	print "Program successfully terminating...."
