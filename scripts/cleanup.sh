#!/bin/bash
#command to run this on directory containing may sub projects (RefFreeDMA run folders):
#find . -name '*.cfg' -exec /data/groups/lab_bock/jklughammer/gitRepos/RefFreeDMA_public/scripts/cleanup.sh {} \;


if [ $# -eq 0 ]
then
	echo "Please pass the configuration file as parameter to cleanup.sh!"
	exit 1
fi
source $1

reduced_dir="$working_dir/reduced"
echo $species
rm $working_dir/${species}_incomplete.flag

if [ ! -s $working_dir/toSelf_filtered*/diffMeth_cpg/*overlap_stats.pdf ]
then
	touch $working_dir/${species}_incomplete.flag
	exit 0
fi



rm $reduced_dir/* 2> /dev/null
rm -r $reduced_dir/forBowtie2
rm $working_dir/conversionCtr/*/*bismark.bam 2> /dev/null
rm -r $working_dir/conversionCtr/*/bismarkTemp 2> /dev/null

find  $reduced_dir/consensus -maxdepth 1 -type f -not -name 'toSelf_filtered_*mm_final.fa' -not -name 'toSelf_filtered_*mm_final.fq' -not -name 'toSelf_filtered_*mm_final_concat.bed' -delete


while [  `sacct --parsable --long|grep '|PENDING|\||RUNNING|'|wc -l` -gt 10 ]; do echo "Too many jobs! --> waiting"; sleep 10m; done
sbatch --export=NONE --get-user-env --job-name="cleanup_fq_$species" --ntasks=1 --cpus-per-task=1 --mem-per-cpu=4000 --partition=longq --time=4-00:00:00 -o /dev/null --wrap="find $working_dir/fastq -name '*trimmed.fq' -exec gzip --best {} \;"
sbatch --export=NONE --get-user-env --job-name="cleanup_bed_$species" --ntasks=1 --cpus-per-task=1 --mem-per-cpu=4000 --partition=longq --time=4-00:00:00 -o /dev/null --wrap="find $working_dir/toSelf_filtered*/ -name '*.bed' -exec gzip -f --best {} \;"

sleep 0.1m
