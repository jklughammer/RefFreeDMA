#! /bin/bash

#Author: Johanna Klughammer
#Date: 26.07.2015


# -----------------------------------------
# Preparations
# -----------------------------------------

analysis_dir=$1
unmapped_fastq=$2
ref_genome_fasta=$3
ref_genome_name=$4
sample=$5
nProcesses=$9
nonCpG=${10}
#------------------------------
#tools
#------------------------------
samtools_path=$6
bsmap_path=$7
biseq_path=$8
bedtools_path=${12}
scripts=${11}
# -----------------------------------------
# bsmap
# -----------------------------------------
function run_bsmap {
	input_fastq=$1
	output_bam=$2
	ref_genome_fasta=$3
	mismatches=$4

	echo "BSMAP Input: $input_fastq"
	echo "BSMAP outout: $output_bam"
	echo "BSMAP ref genome: $ref_genome_fasta"
	
	#BSMAP default params:
		
	# Important options:
	r=1     # -r  [0,1]   how to report repeat hits, 0=none(unique hit/pair only); 1=random one, default:1.
	A=""    # adapter trimming
	D=" "   # "-D C-CGG"    # " -D C-CGG "  # set to " " to disable!
	p=$nProcesses     # -p  <int>   number of processors to use, default=8

	# Other options:
	s=12    # -s  <int>   seed size, default=16(WGBS mode), 12(RRBS mode). min=8, max=16.
	v=$mismatches  # -v  <float> if this value is between 0 and 1, it's interpreted as the mismatch rate w.r.t to the read length.
	#                     otherwise it's interpreted as the maximum number of mismatches allowed on a read, <=15.
	#                     example: -v 5 (max #mismatches = 5), -v 0.1 (max #mismatches = read_length * 10%)
	#                     default=0.08.
	g=0     # -g  <int>   gap size, BSMAP only allow 1 continuous gap (insert or deletion) with up to 3 nucleotides
	#                     default=0
	w=100   # -w  <int>   maximum number of equal best hits to count, <=1000
		# -3          using 3-nucleotide mapping approach
	q=0     # -q  <int>   quality threshold in trimming, 0-40, default=0 (no trim)

	S=1    #1     # -S  <int>   seed for random number generation used in selecting multiple hits normally set to 1, for seed_control set to 5
	f=5     # -f  <int>   filter low-quality reads containing >n Ns, default=5

	n=0     # -n  [0,1]   set mapping strand information. default: 0
		#             -n 0: only map to 2 forward strands, i.e. BSW(++) and BSC(-+),
		#             for PE sequencing, map read#1 to ++ and -+, read#2 to +- and --.
		#             -n 1: map SE or PE reads to all 4 strands, i.e. ++, +-, -+, --

	command="$bsmap_path/bsmap -a $input_fastq -d $ref_genome_fasta -o $output_bam $D -w $w -v $v -r $r -p $p -n $n -S $S -f $f $A -u"	
	echo "COMMAND TO BE RUN: "
	echo $command
	`$command`
	$samtools_path/samtools flagstat $output_bam > $output_bam.flagstat
	$samtools_path/samtools sort -f $output_bam $output_bam
	$samtools_path/samtools index $output_bam
}


function biseq {
	echo "Starting biseqMethilationCalling..."
	
	# Software Requirements:
	# - biopython: wget http://biopython.org/DIST/biopython-1.63.zip
	# - bitarray: https://pypi.python.org/packages/source/b/bitarray/bitarray-0.8.1.tar.gz
	# - guppy: https://pypi.python.org/packages/source/g/guppy/guppy-0.1.10.tar.gz
	# - pysam: https://code.google.com/p/pysam/downloads/detail?name=pysam-0.7.5.tar.gz
	# install locally (to ~/.local/) with: python setup.py install --user 
	
	sample_name=$1
	file=$2
	genome_id=$3
	align_dir=$4
	genome_dir=$(dirname $(dirname $5))
	echo $genome_dir
	mkdir -p $align_dir/biseqMethcalling
	cd $align_dir/biseqMethcalling

command="python $biseq_path/biseqMethCalling.py \
	--sampleName=$sample_name \
	--alignmentFile=$file \
	--methodPrefix=RRBS \
	--includedChromosomes=$genome_id \
	--minFragmentLength=20 \
	--maxFragmentLength=1000 \
	--pfStatus=All \
	--maxMismatches=0.1 \
	--baseQualityScoreC=20 \
	--baseQualityScoreNextToC=10 \
	--laneSpecificStatistics \
	--appendStatisticsOutput=$analysis_dir/$ref_genome_name/RRBS_biseq_statistics.txt \
	--deleteTemp \
	--webOutputDir=RRBS_web \
	--outputDir=. \
	--tempDir=RRBS_temp \
	--timeDelay=0 \
	--genomeFraction=50 \
	--smartWindows=250000 \
	--maxProcesses=$nProcesses \
	--genomeDir=$genome_dir \
	--inGenome=$genome_id \
	--outGenome=$genome_id"

if [ $nonCpG = "TRUE" ]; then
echo "nonCpG mode"
command="$command \
	--processCpHpG \
	--processCpHpH"
fi
echo $command
$command

}


# -----------------------------------------
# Main
# -----------------------------------------

# --- Alignment with BSMAP: ---
align_dir=$analysis_dir/$ref_genome_name/$sample

mkdir -p $align_dir
mismatches=0.08
allow_multimapper=1
bsmap_aligned_bam=$align_dir/$sample.all.aligned.bsmap.$mismatches.mism.$allow_multimapper.r.bam
if [ ! -f $bsmap_aligned_bam ]; then
	run_bsmap $unmapped_fastq $bsmap_aligned_bam $ref_genome_fasta $mismatches
fi

# --- Methylation Analysis: biseq ---
if [ ! -f $align_dir/biseqMethcalling/RRBS_cpgMethylation_$sample.bed ]; then
	biseq  $sample $bsmap_aligned_bam $ref_genome_name $align_dir $ref_genome_fasta
fi

echo "" > $analysis_dir/${sample}.done

#---check read ristribution---
Rscript $scripts/checkReadDistribution.R $align_dir/biseqMethcalling $sample

#---check % fragment coverage
in_bed=$(dirname $ref_genome_fasta).bed
out=${bsmap_aligned_bam/".bam"/".cov"}
#$bedtools_path/bedtools coverage -a $in_bed -b $bsmap_aligned_bam > $out
bedtools coverage -a $in_bed -b $bsmap_aligned_bam -sorted > $out


no_cov=`cut -f11,11 $out|grep "0"|wc -l`
all=`cut -f11,11 $out|wc -l`
ratio=`printf "%.2f" $(echo "100*$no_cov/$all"|bc -l)`

echo -e "$all\t$no_cov\t$ratio" > ${bsmap_aligned_bam/".bam"/".covstats"}







