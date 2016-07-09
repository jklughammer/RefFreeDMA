#!/bin/bash

#Author: Johanna Klughammer
#Date: 26.07.2015

#-----------------------PATHS_START----------------------------
working_dir=$1
bam_dir=$working_dir/unmapped_bam
in_file=$2
maxReadLen=$3
picard_path=$4
trim_galore_path=$5
cutadapt_path=$6
nameSeparator=$7
resMotifs=$8
samtools_path=$9
selected=${10}
#-----------------------PATHS_END----------------------------



#-----------------------PREPS_START----------------------------
#create and set tempdir, sothat processes don't clash
tempdir=$working_dir/TEMP/$(basename $in_file .bam)
mkdir -p $tempdir
export TMPDIR=$tempdir

#-----------------------PREPS_START----------------------------




#-----------------------FUNCTIONS_START------------------------

#picard_sam_to_fastq

function picard_sam_to_fastq {
        input_bam=$1
        output_fastq=$2
        java -Xmx2g -jar $picard_path/SamToFastq.jar I=$input_bam F=$output_fastq INCLUDE_NON_PF_READS=false
}

#trim_galore

function trim_galore {
        input_fastq=$1
        output_dir=$2
        #trimGalore parameters
        #Adapter
        a="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
        #quality trimming
        q=20
        #stringency: Overlap with adapter sequence required to trim a sequence
        s=1
        #Maximum allowed error rate
        e=0.1
        #Minimum Read length
        l=16
        echo $q $s $e $l
        $trim_galore_path/trim_galore -q $q --phred33 -a $a --stringency $s -e $e --length $l --output_dir $output_dir $input_fastq
		
		#dont use --rrbs mode because biseqmethcalling takes care of fill-in Cs
}

#-----------------------FUNCTIONS_END------------------------


#----------------DETERMINE_READ_LENGTH--------------------------------
read=`$samtools_path/samtools view $in_file|head -n1|cut -f10`
readlen="${#read}"
echo "Determined read lenth is $readlen" 


#----------------PREPARE_READS--------------------------------
mkdir -p $working_dir/fastq
# replace # by __ in sample name if there is one (# can cause problems)
old_name=$(basename $in_file .bam)
new_name=${old_name/"$nameSeparator"/"__"}
echo $new_name
out_fastq=$working_dir/fastq/$new_name.fastq
trimmed_fastq=$working_dir/fastq/$(basename $out_fastq .fastq)_trimmed.fq

#bam to fastq conversion
if [ ! -f $out_fastq ]; then
	picard_sam_to_fastq $in_file $out_fastq
fi

#read trimming
if [ ! -f $trimmed_fastq ]; then
	trim_galore $out_fastq $working_dir/fastq/
	
	# additionnally top reads at desired maxReadLen
	if [ $maxReadLen -lt $readlen ]; then
		awk -v top=$maxReadLen '{if(NR%4==2||NR%4==0){print substr($0,1,top)}else{print}}' $trimmed_fastq > $working_dir/fastq/${new_name}_temp
		mv $working_dir/fastq/${new_name}_temp $trimmed_fastq
	fi
fi

truncate -s0 $out_fastq

#check which samples should be included for reference building

sample_name=`echo $new_name|awk 'BEGIN{FS="__"}{print $2}'`
echo $sample_name
printf $selected


#if [ `echo $selected|grep -c $sample_name` == 1 ]
if [[ $selected =~ .*"|$sample_name|".* ]];then
	echo "$sample_name is selected for reference generation!"
	mkdir -p $working_dir/reduced

	#make unique based on unconverted sequences and only include reads that start with TGG in the converted form
	convMotifs=`echo $resMotifs| awk '{gsub(/C/,"T");gsub("[|]","|^");print "^"$1}'`
	echo $convMotifs
	paste <(awk 'NR%4==1' $trimmed_fastq) <(awk 'NR%4==2' $trimmed_fastq) <(awk 'NR%4==2 {gsub(/C/,"T",$1);print $1}' $trimmed_fastq)|awk -v m=$convMotifs '$3 ~ m && $2 !~ /N/'|sort -k 2,2 | uniq -f 1 > $working_dir/reduced/$(basename $out_fastq .fastq)_uniq.ref || exit 1
	if [ ! -s $working_dir/reduced/$(basename $out_fastq .fastq)_uniq.ref ];then
		rm $working_dir/reduced/$(basename $out_fastq .fastq)_uniq.ref		
	fi

else
	echo "$sample_name not selected for reference generation!"
fi


echo "" > $working_dir/$new_name.done

#now check restriction digest efficiency (this is independent of the rest of the pipeline)

total_reads=`$samtools_path/samtools view $in_file|wc -l`
echo "total_reads\t$total_reads" > $working_dir/fastq/$new_name.stats

IFS='|' read -a resMotif_array <<<"$resMotifs"

all_num="counts"
all_motifs="motifs"
all_sum=0
for motif in "${resMotif_array[@]}";do 
	motifComplete="${motif//C/[CT]}"
	motifT="${motif//C/T}"
num=`$samtools_path/samtools view $in_file|grep -P "\t$motifComplete"|wc -l`

num_c=`$samtools_path/samtools view $in_file|grep -P "\t$motif"|wc -l`
num_t=`$samtools_path/samtools view $in_file|grep -P "\t$motifT"|wc -l`

perc_c=`printf "%.2f" $(echo "100*$num_c/$total_reads"|bc -l)`
perc_t=`printf "%.2f" $(echo "100*$num_t/$total_reads"|bc -l)`

all_sum=`echo "$all_sum+$num"|bc`

all_num="$all_num\t$num\t$perc_c\t$perc_t"
all_motifs="$all_motifs\t$motifComplete\t$motif\t$motifT"
done

all_sum_perc=`printf "%.2f" $(echo "100*($total_reads-$all_sum)/$total_reads"|bc -l)`

all_num="$all_num\t$all_sum_perc"
all_motifs="$all_motifs\tothers"

echo $all_motifs>> $working_dir/fastq/$new_name.stats
echo $all_num>> $working_dir/fastq/$new_name.stats

rm -r $tempdir


