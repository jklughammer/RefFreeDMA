#! /bin/bash

#Author: Johanna Klughammer
#Date: 26.07.2015


# -----------------------------------------
# Preparations
# -----------------------------------------
analysis_dir=$1
input_fastq=$2
ref_genome_fasta=$3
ref_genome_name=$4
sample=$5
mismatchRate=$6
nProcesses=$9

#------------------------------
#tools
#------------------------------
samtools_path=$7
bsmap_path=$8

echo $bsmap_path
echo $samtools_path

# -----------------------------------------
# bsmap
# -----------------------------------------
function run_bsmap {

	echo "BSMAP Input: $input_fastq"
	echo "BSMAP ref genome: $ref_genome_fasta"
	
	#BSMAP default params:
		
	# Important options:
	r=1     # -r  [0,1]   how to report repeat hits, 0=none(unique hit/pair only); 1=random one, default:1.
	A=""    # adapter trimming
	D=" "   # "-D C-CGG"    # " -D C-CGG "  # set to " " to disable!
	p=$nProcesses     # -p  <int>   number of processors to use, default=8

	# Other options:
	s=16    # -s  <int>   seed size, default=16(WGBS mode), 12(RRBS mode). min=8, max=16.
	v=$mismatchRate  # -v  <float> if this value is between 0 and 1, it's interpreted as the mismatch rate w.r.t to the read length.
	#                     otherwise it's interpreted as the maximum number of mismatches allowed on a read, <=15.
	#                     example: -v 5 (max #mismatches = 5), -v 0.1 (max #mismatches = read_length * 10%)
	#                     default=0.08.
	g=0     # -g  <int>   gap size, BSMAP only allow 1 continuous gap (insert or deletion) with up to 3 nucleotides
	#                     default=0
	w=100   # -w  <int>   maximum number of equal best hits to count, <=1000
		# -3          using 3-nucleotide mapping approach
	q=0     # -q  <int>   quality threshold in trimming, 0-40, default=0 (no trim)

	S=1    #1     # -S  <int>   seed for random number generation used in selecting multiple hits
	f=5     # -f  <int>   filter low-quality reads containing >n Ns, default=5

	n=0     # -n  [0,1]   set mapping strand information. default: 0
		#             -n 0: only map to 2 forward strands, i.e. BSW(++) and BSC(-+),
		#             for PE sequencing, map read#1 to ++ and -+, read#2 to +- and --.
		#             -n 1: map SE or PE reads to all 4 strands, i.e. ++, +-, -+, --

	output_bam=$align_dir/$sample.bsmap.$mismatchRate.mism.$r.r.$n.n.bam
	echo "BSMAP outout: $output_bam"

	command="$bsmap_path/bsmap -a $input_fastq -d $ref_genome_fasta -o $output_bam $D -w $w -v $v -r $r -p $p -n $n -S $S -f $f $A -u"
	echo "COMMAND TO BE RUN: "
	echo $command
	`$command`
	$samtools_path/samtools flagstat $output_bam > $output_bam.flagstat
}



# -----------------------------------------
# Main
# -----------------------------------------

# --- Alignment with BSMAP: ---
align_dir=$analysis_dir/crossMapping/$ref_genome_name
mkdir -p $align_dir

run_bsmap


