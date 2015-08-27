#! /bin/bash

#Author: Johanna Klughammer
#Date: 26.07.2015


cons_dir=$1
sample=$2
working_dir=$3

grep -v -P "\t[CT]GG.*[CT][CT]G\t" $cons_dir/${sample}_final>$cons_dir/${sample}_final_norc
grep -P "\t[CT]GG.*[CT][CT]G\t" $cons_dir/${sample}_final| sort -k 2,2 -t$'\t'|awk '{ print length($2), $0 | "sort -n" }'|cut -d" " -f2 >$cons_dir/${sample}_final_rc


touch $working_dir/findRevComp.done	