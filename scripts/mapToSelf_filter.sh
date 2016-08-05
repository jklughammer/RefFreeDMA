#! /bin/bash

#Author: Johanna Klughammer
#Date: 26.07.2015


in_fastq=$1
reduced_dir=$2
filter=$3
working_dir=$4

sort -k1,1 -t$'\t' ${in_fastq/.fq/.ref}> ${in_fastq/.fq/_sorted.ref}

samtools view -S -F 4 $reduced_dir/$(basename $in_fastq .fq)_toSelf.sam|awk -v filter=$filter -F'\tXM:i:' '{split($1,a,"\t") split($2,b,"\t"); if(a[2]&&0x100==256) type="sec";else type="prim"; if(b[1]/length(a[10])<=filter) print a[1]"\t"a[10]"\t"a[3]"\t"a[4]"\t"a[6]"\t"b[1]"\t"length(a[10])"\t"type"\t"a[2]}'|sort -T $working_dir/TEMP -k1,1 -t$'\t'|join -j 1 -t$'\t' - ${in_fastq/.fq/_sorted.ref}| awk 'BEGIN {OFS="\t"} { $2 = $10 ; print }'|cut -f1,2,3,4,5,6,7,8,9>$reduced_dir/toSelf_filtered_${filter}mm


touch $working_dir/mapToSelf_filter.done	
