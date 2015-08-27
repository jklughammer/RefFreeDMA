#!/bin/bash

#Author: Johanna Klughammer
#Date: 26.07.2015


#-----------------------PATHS_START----------------------------
reduced_dir=$1
ref_fasta=$2
in_fastq=$3
working_dir=$4
bowtie2_path=$5
nProcesses=$6
#-----------------------PATHS_END----------------------------


#-----------------------FUNCTIONS_START------------------------

#________bowtie2_______
function make_bt2idx {
	input_fasta=$1
	prefix=$2

	$bowtie2_path/bowtie2-build $input_fasta $prefix
}

function run_bt2 {
	reference=$1
	fastq=$2
	out_sam=$3

	N=1
	L=22
	nCeil="L,0,0.2"
	mp=3
	np=0
	scoreMin="L,-0.6,-0.6"
	k=300
	D=3
	rdg="20,20"
	rfg="20,20"
	p=$nProcesses

	command="bowtie2 -t -q --phred33 --end-to-end -N $N -L $L --norc --n-ceil $nCeil --mp $mp  --np $np --score-min $scoreMin -k $k -D $D --rdg $rdg --rfg $rfg  -p $p  -x $reference -U $fastq -S $out_sam"
#	command="bowtie2 -t -q --phred33 --end-to-end -N $N -L $L --n-ceil $nCeil --mp $mp  --np $np --score-min $scoreMin -k $k -D $D --rdg $rdg --rfg $rfg  -p $p  -x $reference -U $fastq -S $out_sam"

`$command`
}


#-----------------------FUNCTIONS_END------------------------



#----------------EXECUTE_BT2--------------------------------
mkdir -p $reduced_dir/forBowtie2
cp $ref_fasta $reduced_dir/forBowtie2/
prefix=$reduced_dir/forBowtie2/$(basename $ref_fasta .fa)

make_bt2idx $prefix.fa $prefix

run_bt2 $prefix $in_fastq $reduced_dir/$(basename $in_fastq .fq)_toSelf.sam || exit 1

echo "" > $working_dir/redToSelf.done







