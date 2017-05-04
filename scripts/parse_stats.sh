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

#samples=`awk 'NR==1{for(i=1;i<=NF;i++){if($i=="Sample_Name"){n=i}}} NR>1{print $n}' $sample_annotation` 
#echo "Samples: "$samples

#remove preexisting file
rm $summary_dir/summary.txt

#get number of fragments in the deduced genome
frag_file=$working_dir/reduced/consensus/toSelf_filtered_${mapToSelf_filter}mm_final.fa
total_frags=`grep -v ">" $frag_file| wc -l`

IFS='|' read -a resMotif_array <<<"$restrictionSites"

all_num="counts"
all_motifs="motifs"
all_sum=0
for motif in "${resMotif_array[@]}";do 
	motifComplete="${motif//C/[CT]}"
	motifT="${motif//C/T}"
num=`grep -P "^$motifComplete" $frag_file|wc -l`

num_c=`grep -P "^$motif" $frag_file|wc -l`
num_t=`grep -P "^$motifT" $frag_file|wc -l`

perc_c=`printf "%.2f" $(echo "100*$num_c/$total_frags"|bc -l)`
perc_t=`printf "%.2f" $(echo "100*$num_t/$total_frags"|bc -l)`

all_sum=`echo "$all_sum+$num"|bc`

all_num="$all_num\t$num\t$perc_c\t$perc_t"
all_motifs="$all_motifs\t$motifComplete\t$motif\t$motifT"
done

all_sum_perc=`printf "%.2f" $(echo "100*($total_frags-$all_sum)/$total_frags"|bc -l)`

#samples used for reference generation
shopt -s extglob
num_samples=`ls $working_dir/reduced/!(merged)_uniq.ref|wc -l`
shopt -u extglob

all_num="$all_num\t$total_frags\t$all_sum_perc\t$num_samples\t$species"
all_motifs="$all_motifs\ttotal\tothers\tsamples\tspecies"

echo -e $all_motifs> $summary_dir/ref_summary.txt
echo -e $all_num>> $summary_dir/ref_summary.txt

#now analyse each sample
dir=$working_dir/$analysis_dir
shopt -s extglob
subdirs=`ls -d $dir/!(diffMeth_*)/`
shopt -u extglob

for sample in $subdirs; do
sample=$(basename $sample)
sample_name=${sample/*__/""}
echo $sample
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

#fragment coverage
fragment_cov=`cat $working_dir/$analysis_dir/$sample/*.covstats`

#contamination stats
if [ -f $working_dir/fastq/${sample}_trimmed_bwaMeth_decon_summary.txt ]; then
decon_res=$(Rscript -e "library(data.table);library(splitstackshape);decon_stats=cSplit(fread('$working_dir/fastq/${sample}_trimmed_bwaMeth_decon_summary.txt'),'V1',' ');decon_stats[,V4:=gsub('gi.*ref','',V1_2),];decon_stats_col=decon_stats[-1,][,sum(V1_1),by=V4];max_cont=as.character(decon_stats_col[which.max(V1),'V1',with=FALSE]);max_cont_sp=as.character(decon_stats_col[which.max(V1),'V4',with=FALSE]);cont=as.numeric(decon_stats[-1,sum(V1_1)]);cont_rat=cont/(cont+as.numeric(decon_stats[1,'V1_1',with=FALSE]));cat(paste0(c(max_cont_sp,max_cont,cont,cont_rat),collapse='\t'))")
else
decon_res="NA\tNA\tNA\tNA"
fi


if [ ! -f $summary_dir/summary.txt ]; then
	echo -e "sample\tspecies\ttotal_reads\tmapped_reads\tmapping_efficiency\tinformative_reads\tCpG_meth\tavg_meth\tCpG_measurements\tcoveredCpGs\tconversionRate\tk1_unmeth\tk3_meth\ttotalMeasurements_k1\ttotalMeasurements_k3\ttotal_reads_untrimmed\t$motifs\tfragments_ref\tfragments_uncovered\tfragments_uncovered_perc\tmax_cont_sp\tmax_cont\tcont\tcont_rat" >$summary_dir/summary.txt
fi

echo -e "$sample_name\t$species\t$total_reads\t$mapped_reads\t$mapping_eff\t$informative_reads\t$CpG_meth\t$avg_meth\t$CpGMeasurements\t$coveredCpGs\t$conversionRate\t$k1_unmeth\t$k3_meth\t$totalMeasurements_k1\t$totalMeasurements_k3\t$total_reads_untrimmed\t$counts\t$fragment_cov\t$decon_res" >>$summary_dir/summary.txt

done
