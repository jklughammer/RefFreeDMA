#!/bin/bash

if [ $# -eq 0 ]
then
	echo "Please pass the configuration file as parameter to parse_stats.sh!"
	exit 1
fi
source $1

analysis_dir="toSelf_filtered_${mapToSelf_filter}mm_final_concat"
convCtr_dir="conversionCtr"

summary_dir=$working_dir


dir=$working_dir/$analysis_dir
for sample in `ls -d $dir/BSF*`; do
sample=$(basename $sample)
sample_name=`echo $sample | awk -F'__' '{print $2}'`

#mapped flagstat
total_reads=`awk 'NR==1' $dir/$sample/*flagstat|cut -f1 -d " "`
mapped_reads=`awk 'NR==5' $dir/$sample/*flagstat|cut -f1 -d " "`
mapping_eff=`awk 'NR==5' $dir/$sample/*flagstat|cut -f5 -d " "|cut -f2 -d "("|cut -f1 -d ":"`
mapping_eff=${mapping_eff/"%"/}

#biseq stats
informative_reads=`awk 'NR==2' $dir/$sample/biseqMethcalling/RRBS_statistics*|cut -f6`
CpG_meth=`awk 'NR==2' $dir/$sample/biseqMethcalling/RRBS_statistics*|cut -f13` #don't use firstMotifMethylationMean (12) but uniqueMotifMethylationMean 
avg_meth=`awk 'NR==2' $dir/$sample/biseqMethcalling/RRBS_statistics*|cut -f10`
coveredCpGs=`awk 'NR==2' $dir/$sample/biseqMethcalling/RRBS_statistics*|cut -f8`
conversionRate=`awk 'NR==2' $dir/$sample/biseqMethcalling/RRBS_statistics*|cut -f9`
CpGMeasurements=`awk 'NR==2' $dir/$sample/biseqMethcalling/RRBS_statistics*|cut -f20`

#conversion controls
if [ -d $working_dir/$convCtr_dir/$sample/ ]; then
k1_unmeth=`cut -f 3 -d , $working_dir/$convCtr_dir/$sample/*K1_unmethylated.txt`
k3_meth=`cut -f 3 -d , $working_dir/$convCtr_dir/$sample/*K3_methylated.txt`
totalMeasurements_k1=`cut -f 20 -d , $working_dir/$convCtr_dir/$sample/*K1_unmethylated.txt`
totalMeasurements_k3=`cut -f 20 -d , $working_dir/$convCtr_dir/$sample/*K3_methylated.txt`
else
echo "No external conversion contol found for $sample."
k1_unmeth="NA"
k3_meth="NA"
totalMeasurements_k1="NA"
totalMeasurements_k3="NA"
fi

#restriction digest
total_reads_untrimmed=`grep "total" $working_dir/fastq/$sample.stats`
motifs=`grep "motif" $working_dir/fastq/$sample.stats`
counts=`grep "counts" $working_dir/fastq/$sample.stats`

total_reads_untrimmed=${total_reads_untrimmed/"total_reads\t"/}
motifs=${motifs/"motifs\t"/}
counts=${counts/"counts\t"/}

if [ ! -f $summary_dir/summary.txt ]; then
	echo -e "sample\tspecies\ttotal_reads\tmapped_reads\tmapping_efficiency\tinformative_reads\tCpG_meth\tavg_meth\tCpG_measurements\tcoveredCpGs\tconversionRate\tk1_unmeth\tk3_meth\ttotalMeasurements_k1\ttotalMeasurements_k3\ttotal_reads_untrimmed\t$motifs" >$summary_dir/summary.txt
fi

echo -e "$sample_name\t$species\t$total_reads\t$mapped_reads\t$mapping_eff\t$informative_reads\t$CpG_meth\t$avg_meth\t$CpGMeasurements\t$coveredCpGs\t$conversionRate\t$k1_unmeth\t$k3_meth\t$totalMeasurements_k1\t$totalMeasurements_k3\t$total_reads_untrimmed\t$counts" >>$summary_dir/summary.txt

done
