#! /bin/bash

#Author: Johanna Klughammer
#Date: 26.07.2015


cons_dir=$1
sample=$2
working_dir=$3
resMotifs_orig=$4
resMotifs=${resMotifs_orig//"|"/"\|"}

#grep -v -P "\t[CT]GG.*[CT][CT]G\t" $cons_dir/${sample}_final>$cons_dir/${sample}_final_norc
#grep -P "\t[CT]GG.*[CT][CT]G\t" $cons_dir/${sample}_final| sort -k 2,2 -t$'\t'|awk '{ print length($2), $0 | "sort -n" }'|cut -d" " -f2 >$cons_dir/${sample}_final_rc

pattern=$(echo $resMotifs| awk '{"echo " $0|getline v0; "echo " $0 "| rev | tr ATGC TACG"|getline v1; v2="("v0").*("v1")"; v3=gensub("C","[CT]","g",v2);print "\\t"v3"\\t" }')

echo "pattern: " $pattern
grep -v -P $pattern $cons_dir/${sample}_final>$cons_dir/${sample}_final_norc
grep -P $pattern $cons_dir/${sample}_final| sort -k 2,2 -t$'\t'|awk '{ print length($2), $0 | "sort -n" }'|cut -d" " -f2 >$cons_dir/${sample}_final_rc



touch $working_dir/findRevComp.done	
