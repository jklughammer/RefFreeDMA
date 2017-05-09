#! /bin/bash

#Author: Johanna Klughammer
#Date: 26.07.2015


#-----------------------LOAD_CONFIGURATION_START---------------
if [ $# -eq 0 ]
then
	echo "Please pass the configuration file as parameter to RefFreeDMA.sh!"
	exit 1
fi
#turn decontamination off by default (overwrite through configuration file)
decon=FALSE
source $1

if [ ! -n "$unconv_tag" ];then 
	unconv_tag="no unconv_tag"
fi

#-----------------------LOAD_CONFIGURATION_END---------------

#-----------------------PATHS_START----------------------------
bam_dir=$working_dir/unmapped_bam
scripts=$(cd "$(dirname $0)"; pwd)/scripts
logdir=$working_dir/log
mkdir -p $logdir
#-----------------------PATHS_END------------------------------

#-----------------------TOOLS_START----------------------------
biseq_path=$scripts/
export PATH=$cutadapt_path/bin:$picard_path:$trim_galore_path:$bowtie2_path:$bsmap_path:$samtools_path:$bedtools_path:$bwa_path:$PATH
export PYTHONPATH=$tool_path/python2.7:$cutadapt_path/lib/python2.7/site-packages:~/.local/lib/python2.7/site-packages/:$PYTHONPATH

#check if python libraries are there
printf "\nChecking for required python libraries...\n"
fail=0
python -c 'import Bio' 2>/dev/null && echo "python Bio ... OK" ||{ echo "python Bio ...  FAIL";fail=1; }
python -c 'import pysam' 2>/dev/null && echo "python pysam ... OK" ||{ echo "python pysam ...  FAIL";fail=1; }
python -c 'import bitarray' 2>/dev/null && echo "python bitarray ... OK" ||{ echo "python bitarray ...  FAIL";fail=1; }
python -c 'import guppy' 2>/dev/null && echo "python guppy ... OK" ||{ echo "python guppy ...  FAIL";fail=1; }
if [ $decon = "TRUE" ];then
python -c 'import toolshed' 2>/dev/null && echo "python toolshed ... OK" ||{ echo "python toolshed ...  FAIL";fail=1; }
fi
#check if tools are there
printf "\nChecking for required tools...\n"
which samtools &>/dev/null && echo "samtools ... OK" ||{ echo "samtools ...  FAIL";fail=1; }
which bedtools &>/dev/null && echo "bedtools ... OK" ||{ echo "bedtools ...  FAIL";fail=1; }
which trim_galore &>/dev/null && echo "trim_galore ... OK" ||{ echo "trim_galore ...  FAIL";fail=1; }
which cutadapt &>/dev/null && echo "cutadapt ... OK" ||{ echo "cutadapt ...  FAIL";fail=1; }
which bowtie2 &>/dev/null && echo "bowtie2 ... OK" ||{ echo "bowtie2 ...  FAIL";fail=1; }
which bsmap &>/dev/null && echo "bsmap ... OK" ||{ echo "bsmap ...  FAIL";fail=1; }
[ -e $picard_path/SamToFastq.jar ]  && echo "samToFastq.jar ... OK" ||{ echo "samToFastq.jar ...  FAIL";fail=1; } 
if [ $decon = "TRUE" ];then
which bwa &>/dev/null && echo "bwa ... OK" ||{ echo "bwa ...  FAIL";fail=1; }
python $bwameth_path/bwameth.py -h &>/dev/null && echo "bwameth ... OK" ||{ echo "bwameth ...  FAIL";fail=1; }
fi


if [ $fail -eq 1 ];then
	printf "\nPlease provide missing python libraries and tools.\nIf you are using the external software bundle, did you set the path in the configuration file?\n"
	exit 0
fi
#-----------------------TOOLS_END----------------------------


#-----------------------FUNCTIONS_START------------------------
function wait_for_slurm {
	wait_time=$1
	submitted=$2
	working_dir=$3
	
	echo "wait_time $wait_time"
	echo "submitted $submitted"
#	echo "wd $working_dir"
	
	
	while [ `ls $working_dir/*.done 2>/dev/null|wc -w` -lt 1 ]
	do
		echo "waiting for first job to finish..."
		sleep ${wait_time}m
	done

	no=`ls $working_dir/*.done|wc -w`
	while [ $no -lt $submitted ]
	do
		no=`ls $working_dir/*.done|wc -w`
		echo "waiting till all jobs are finished... $no done"
		sleep ${wait_time}m

	done
	
	rm $working_dir/*.done
}

function get_proc_stats {
echo $2 >> $logdir/time
eval "/usr/bin/time -o $logdir/time --append -f 'command: %C\ntime(sec): %e(real_time) %S(system) %U(user)\nmemory(Kbytes): %M(peak) %K(avg_total)\n' $1"

}

#-----------------------FUNCTIONS_END------------------------

#-----------------------BASIC_CHECKS_START-------------------------

if [[ $working_dir =~ "__" ]]; then 
echo "The string '__' was found in the working_dir path. This is not allowed!"
exit 1
fi

#-----------------------BASIC_CHECKS_END-------------------------


#------------------READ_PREPARATION_START-----------------------


#convert .bam files to .fastq files and trimm reads (quality and adapter)
step="\n-------Read preparation-------\n"
printf "$step"

selected=`awk 'NR==1{for(i=1;i<=NF;i++){if($i=="Select"){s=i};if($i=="Sample_Name"){n=i}}} NR>1{if($s==1) {printf "|"$n"|"}}' $sample_annotation`
sel_count=`echo $selected|grep -o "|"|wc -l`;number_selected=$((sel_count / 2))

printf "$number_selected samples selected for reference generation:\n$selected\n"
printf "Unconverted tag: $unconv_tag\n"


if [ $decon = "TRUE" ];then
	echo "Running in decontamination mode. This will take a lot of memory and most likely fail if run locally."
	mem=100000
	queue="longq"
	time="2-00:00:00"
	if [ ! -s $decon_reference.bwameth.c2t.bwt ]; then
	echo "Cound not find decon_reference. Exiting!"
	exit 1
	fi
else
	mem=6000
	queue="shortq"
	time="08:00:00"
	decon_reference="NA"
	bwameth_path="NA"
fi

count=0
submitted=0
rm $working_dir/*.done 2>/dev/null
for file in `ls $bam_dir/*.bam`; do
	sample=$(basename $file .bam)
	sample=${sample/"$nameSeparator"/"__"}
	echo $sample
	if [ ! -s $working_dir/reduced/*${sample}_uniq.ref* ] ; then
		#create tempdir, sothat processes don't clash
		tempdir=$working_dir/TEMP/$sample
		mkdir -p $tempdir
		echo "submitted"
		if [ $parallel = "TRUE" ]; then
			sbatch --export=ALL --get-user-env --job-name=prepareReads_$sample --workdir $tempdir --ntasks=1 --cpus-per-task=1 --mem-per-cpu=$mem --partition=$queue --time=$time -e $logdir/prepareReads_${sample}_%j.err -o $logdir/prepareReads_${sample}_%j.log $scripts/prepareReads.sh $working_dir $file $maxReadLen $picard_path $trim_galore_path $cutadapt_path "$nameSeparator" $restrictionSites $samtools_path "$selected" $decon $bwameth_path $decon_reference $tempdir $unconv_tag
			sleep 0.01m
			((submitted++))
		else
			get_proc_stats "$scripts/prepareReads.sh $working_dir $file $maxReadLen $picard_path $trim_galore_path $cutadapt_path '$nameSeparator' '$restrictionSites' '$samtools_path' '$selected' $decon $bwameth_path $decon_reference $tempdir $unconv_tag &> $logdir/prepareReads_${sample}.log" "$step"
		fi
	else
		echo "${sample}_uniq.ref* exists. Not submitted!"
	fi
	((count++))
done
if [ ! $submitted = 0 ]; then
	wait_for_slurm $wait_time $submitted $working_dir
fi
shopt -s extglob
if [ ! `ls $working_dir/reduced/!(merged)_uniq.ref*|wc -l` = $number_selected ] || [ ! `ls $working_dir/fastq/*_trimmed.fq|wc -l` = $count ]; then
	echo "Didn't find the expected number of uniq.ref files ($number_selected) or trimmed.fq files ($count). Exiting!"
	exit 1
fi
shopt -u extglob

#------------------READ_PREPARATION_END-----------------------


#------------------REFERENCE_GENERATION_START-----------------------

step="\n-------Noise reduction-------\n"
count=0
printf "$step"
if [ ! -f $working_dir/reduced/merged_dupl.ref ]; then
	if [ $parallel = "TRUE" ]; then
		sbatch --export=ALL --get-user-env --job-name=noiseReduction --ntasks=1 --cpus-per-task=1 --mem-per-cpu=4000 --partition=shortq --time=08:00:00 -e $logdir/noiseReduction_%j.err -o $logdir/noiseReduction_%j.log $scripts/noiseReduction.sh $working_dir $filtLim $maxSamples $unconv_tag
		((count++))
	else
		get_proc_stats "$scripts/noiseReduction.sh $working_dir $filtLim $maxSamples $unconv_tag" "$step"
	fi
else
	echo "merged_dupl.ref exists. Skipping!"
fi
if [ ! $count = 0 ]; then
	wait_for_slurm $wait_time $count $working_dir
fi

if [ ! -s $working_dir/reduced/merged_dupl.ref ]; then
	echo "Noise reduction failed. Exiting!"
	exit 1
fi


#reduce redundency further by including smaller in larger reads and by merging reads with the same converted sequence
step="\n-------Basic consensus finding-------\n"
printf "$step"

count=0
if [ `ls $working_dir/reduced/*collapse* 2>/dev/null|wc -l` -lt 1 ]; then
	
	if [ $parallel = "TRUE" ]; then
		sbatch --export=ALL --get-user-env --job-name=makePreConsensus --ntasks=1 --cpus-per-task=1 --mem-per-cpu=4000 --partition=shortq --time=08:00:00 -e $logdir/makePreConsensus_%j.err -o $logdir/makePreConsensus_%j.log $scripts/makePreConsensus.py $working_dir/reduced/merged_uniq.ref $maxReadLen $maxReadLen $working_dir $cLimit $unconv_tag
		((count++))
	else
		get_proc_stats "python $scripts/makePreConsensus.py $working_dir/reduced/merged_uniq.ref $maxReadLen $maxReadLen $working_dir $cLimit $unconv_tag" "$step"
	fi
else
	echo "At least one collapse file exists. Skipping!"
fi
if [ ! $count = 0 ]; then
	wait_for_slurm $wait_time $count $working_dir
fi

if [ ! -s $working_dir/reduced/*collapse.fa ] || [ ! -s $working_dir/reduced/*collapse.fq ]; then
	echo "Basic consensus finding failed. Exiting!"
	exit 1
fi

#reduce redundency of deduced_references further by now also including SNPs
#prepare bowtie2 output to merge similar refs in order to reduce the redundancy --> makeFinalConsensus.R
#prduce table of form
#1.     Mapped_ref_name
#2.     mapped_ref_seq
#3.     target_ref_name
#4.     target_ref_pos
#5.     CIGAR
#6.     Mismatches
#7.     mapped_ref_length
#8.     alignment(primary|secondary)
#END
step="\n-------Map to self with Bowtie2-------\n"
printf "$step"
sample="merged_uniq_collapse"
reduced_dir=$working_dir/reduced/
ref_fasta=`ls $reduced_dir/*collapse.fa` || exit 1
in_fastq=`ls $reduced_dir/*collapse.fq` || exit 1


count=0
if [ ! -f $reduced_dir/*_toSelf.sam ]; then
	rm $working_dir/*.done 2>/dev/null
	if [ $parallel = "TRUE" ]; then
		sbatch --export=ALL --get-user-env --job-name=redToSelf --ntasks=1 --cpus-per-task=$nProcesses --mem-per-cpu=4000 --partition=shortq --time=08:00:00 -e $logdir/redToSelf_%j.err -o $logdir/redToSelf_%j.log $scripts/mapToSelf_bt2.sh $reduced_dir $ref_fasta $in_fastq $working_dir $bowtie2_path $nProcesses
		count=1
	else
		get_proc_stats "$scripts/mapToSelf_bt2.sh $reduced_dir $ref_fasta $in_fastq $working_dir $bowtie2_path $nProcesses &> $logdir/redToSelf_${sample}.log" "$step"
	fi
else
	echo "toSelf.sam exists. Skipping!"
fi

if [ ! $count = 0 ]; then
	wait_for_slurm $wait_time $count $working_dir
fi

if [ ! -s $reduced_dir/*_toSelf.sam ]; then
	echo "Map to self failed. Exiting!"
	exit 1
fi

# filter bt2 alignment (max no of mismatches allowed) and produce .ref file Takes about 5 h for 800M lines (~12M refs)
step="\n-------MapToSelf filter-------\n"
count=0
printf "$step"
if [ ! -f $reduced_dir/toSelf_filtered_${mapToSelf_filter}mm ]; then
	if [ $parallel = "TRUE" ]; then
		sbatch --export=ALL --get-user-env --job-name=mapToSelf_filter --ntasks=1 --cpus-per-task=1 --mem-per-cpu=6000 --partition=mediumq --time=1-00:00:00 -e $logdir/mapToSelf_filter_%j.err -o $logdir/mapToSelf_filter_%j.log $scripts/mapToSelf_filter.sh $in_fastq $reduced_dir $mapToSelf_filter $working_dir
		((count++))
	else
		get_proc_stats "$scripts/mapToSelf_filter.sh $in_fastq $reduced_dir $mapToSelf_filter $working_dir" "$step"
	fi
else
	echo "toSelf_filtered_${mapToSelf_filter}mm exists. Skipping!" 
fi

if [ ! $count = 0 ]; then
	wait_for_slurm $wait_time $count $working_dir
fi

if [ ! -s $reduced_dir/toSelf_filtered_${mapToSelf_filter}mm ]; then
	echo "Map to self filtering failed. Exiting!"
	exit 1
fi

# Account for SNPs for ~300GB give 90000 mem per cpu
step="\n-------Accurate consensus finding-------\n"
printf "$step"
sample=toSelf_filtered_${mapToSelf_filter}mm
cons_dir=$reduced_dir/consensus/
mkdir -p $cons_dir
count=0
if [ ! -f $cons_dir/toSelf_*_final ]; then
	rm $working_dir/*.done 2>/dev/null
	if [ $parallel = "TRUE" ]; then
		sbatch --export=ALL --get-user-env --job-name=makeConsensus --ntasks=1 --cpus-per-task=1 --mem-per-cpu=90000 --partition=mediumq --time=42:00:00 -e $logdir/makeConsensus_%j.err -o $logdir/makeConsensus_%j.log $scripts/makeFinalConsensus.R $sample $cons_dir $reduced_dir $working_dir $consensus_dist $cLimit
		count=1
	else
		get_proc_stats "$scripts/makeFinalConsensus.R $sample $cons_dir $reduced_dir $working_dir $consensus_dist $cLimit &> $logdir/makeConsensus_${sample}.log" "$step"
	fi
else
	echo "final file exists. Skipping!"
fi

if [ ! $count = 0 ]; then
	wait_for_slurm $wait_time $count $working_dir
fi

if [ ! -s $cons_dir/toSelf_*_final ]; then
	echo "Accurate consensus finding failed. Exiting!"
	exit 1
fi

#deal with reverse complements
step="\n-------Identify potential reverse complements-------\n"
count=0
printf "$step"
if [ ! -f $cons_dir/${sample}_final_rc ]; then
	if [ $parallel = "TRUE" ]; then
		sbatch --export=ALL --get-user-env --job-name=findRevComp_filter --ntasks=1 --cpus-per-task=1 --mem-per-cpu=4000 --partition=shortq --time=08:00:00 -e $logdir/findRevComp_%j.err -o $logdir/findRevComp_%j.log $scripts/findRevComp.sh $cons_dir $sample $working_dir $restrictionSites
		((count++))
	else
		get_proc_stats "$scripts/findRevComp.sh $cons_dir $sample $working_dir '$restrictionSites'" "$step"
	fi
else
	echo "${sample}_final_rc exists. Skipping!"
fi
if [ ! $count = 0 ]; then
	wait_for_slurm $wait_time $count $working_dir
fi

if [ ! -s $cons_dir/${sample}_final_rc ]; then
	echo "Identifying potential reverese complements failed. Exiting!"
	exit 1
fi

step="\n-------Merge reverse complements-------\n"
count=0
printf "$step"
if [ ! -f $cons_dir/${sample}_final_rc-res ]; then
	if [ $parallel = "TRUE" ]; then
		sbatch --export=ALL --get-user-env --job-name=mergeRevComp --ntasks=1 --cpus-per-task=1 --mem-per-cpu=4000 --partition=shortq --time=08:00:00 -e $logdir/mergeRevComp_%j.err -o $logdir/mergeRevComp_%j.log $scripts/mergeRevComp.py $cons_dir/${sample}_final_rc $working_dir
		((count++))
	else
		get_proc_stats "python $scripts/mergeRevComp.py $cons_dir/${sample}_final_rc $working_dir" "$step"
	fi
	if [ ! $count = 0 ]; then
		wait_for_slurm $wait_time $count $working_dir
	fi

	cat $cons_dir/${sample}_final_rc-res >> $cons_dir/${sample}_final_norc
else
	echo "${sample}_final_rc-res exists. Skipping!"
fi	

if [ ! -s $cons_dir/${sample}_final_rc-res ]; then
	echo "Merging potential reverese complements failed. Exiting!"
	exit 1
fi

#concatenate give 10000 for genomes with ~5M fragments
step="\n-------Concatenating deduced genome fragments-------\n"
count=0
printf "$step"
if [ ! -f $cons_dir/${sample}_final_concat/${sample}_final_concat.fa ]; then
	if [ $parallel = "TRUE" ]; then
		sbatch --export=ALL --get-user-env --job-name=concatenateRef --ntasks=1 --cpus-per-task=1 --mem-per-cpu=10000 --partition=shortq --time=08:00:00 -e $logdir/concatenateRef_%j.err -o $logdir/concatenateRef_%j.log $scripts/concatenateRef.py $cons_dir/${sample}_final_norc $maxReadLen $working_dir
		((count++))
	else
		get_proc_stats "python $scripts/concatenateRef.py $cons_dir/${sample}_final_norc $maxReadLen $working_dir" "$step"
	fi
	
	if [ ! $count = 0 ]; then
	wait_for_slurm $wait_time $count $working_dir
	fi
	
	mkdir -p $cons_dir/${sample}_final_concat
	if [ -f $cons_dir/${sample}_final_concat.fa ]; then
		mv $cons_dir/${sample}_final_concat.fa $cons_dir/${sample}_final_concat/
		
	fi 
else
	echo "${sample}_final_concat.fa exists. Skipping!"
fi

if [ ! -s $cons_dir/${sample}_final_concat/${sample}_final_concat.fa ]; then
	echo "Concatenating failed. Exiting!"
	exit 1
fi

#------------------REFERENCE_GENERATION_END-----------------------

#---------------------CrossMapping_START--------------------------
if [ ! $cross_genome_fa = "-" ]; then
	ref_genome_fasta=$cross_genome_fa
	cross_genome_id=$(basename $ref_genome_fasta .fa)
	step="\n-------Cross mapping to $cross_genome_id-------\n"
	printf "$step"
	if [ ! -s $working_dir/crossMapping/$cross_genome_id/$sample*.bam ]; then
		unmapped_fastq=$working_dir/reduced/consensus/${sample}_final.fq
		if [ $parallel = "TRUE" ]; then
			sbatch --export=ALL --get-user-env --job-name=crossMapping_$cross_genome_id --ntasks=1 --cpus-per-task=4 --mem-per-cpu=4000 --partition=shortq --time=08:00:00 -e "$logdir/crossMapping_${cross_genome_id}_%j.err" -o "$logdir/crossMapping_${cross_genome_id}_%j.log" $scripts/crossMapping.sh $working_dir $unmapped_fastq $ref_genome_fasta $cross_genome_id $sample $crossMap_mismatchRate $samtools_path $bsmap_path $nProcesses
		else
			get_proc_stats "$scripts/crossMapping.sh $working_dir $unmapped_fastq $ref_genome_fasta $cross_genome_id $sample $crossMap_mismatchRate $samtools_path $bsmap_path $nProcesses &> $logdir/crossMapping_${cross_genome_id}.log" "$step"
		fi
	else
		echo "Cross-mapping flagstat file already exists. Skipping!"
	fi
fi
#------------------CrossMapping_END-----------------------

#---------------------Sample CrossMapping_START--------------------------
if [ ! $cross_genome_fa = "-" ]; then
	ref_genome_fasta=$cross_genome_fa
	cross_genome_id=$(basename $ref_genome_fasta .fa)
	step="\n-------Sample Cross mapping to $cross_genome_id-------\n"
	printf "$step"
	
	for unmapped_fastq in `ls $working_dir/fastq/*trimmed.fq`; do
	run_sample=$(basename $unmapped_fastq .fq)
	run_sample=${run_sample//_trimmed/}
	echo $run_sample
	if [ ! -s $working_dir/crossMapping/$cross_genome_id/$run_sample*.bam ]; then
		if [ $parallel = "TRUE" ]; then
			sbatch --export=ALL --get-user-env --job-name=crossMapping_$run_sample --ntasks=1 --cpus-per-task=4 --mem-per-cpu=4000 --partition=shortq --time=08:00:00 -e "$logdir/crossMapping_${run_sample}_%j.err" -o "$logdir/crossMapping_${run_sample}_%j.log" $scripts/crossMapping.sh $working_dir $unmapped_fastq $ref_genome_fasta $cross_genome_id $run_sample $crossMap_mismatchRate $samtools_path $bsmap_path $nProcesses
		else
			get_proc_stats "$scripts/crossMapping.sh $working_dir $unmapped_fastq $ref_genome_fasta $cross_genome_id $run_sample $crossMap_mismatchRate $samtools_path $bsmap_path $nProcesses &> $logdir/crossMapping_${run_sample}.log" "$step"
		fi
	else
		echo "Cross-mapping bam file already exists. Skipping!"
	fi
	done
fi
#------------------Sample CrossMapping_END-----------------------


#---------------------READ_Mapping_and_Methcalling_START--------------------------
step="\n-------Mapping to the deduced genome and methylation calling-------\n"
printf "$step"
mode="ded_ref"
genome_id=${sample}_final_concat
ref_genome_fasta=$working_dir/reduced/consensus/$genome_id/$genome_id.fa 
 
count=0
submitted=0
rm $working_dir/*.done 2>/dev/null
for unmapped_fastq in `ls $working_dir/fastq/*trimmed.fq`; do
	sample=$(basename $unmapped_fastq .fq)
	sample=${sample//_trimmed/}
	echo $sample
	if [ ! -f $working_dir/$genome_id/$sample/biseqMethcalling/*cpgMethylation*.bed* ]; then
		echo submitted
		if [ $parallel = "TRUE" ]; then
			sbatch --export=ALL --get-user-env --job-name=meth_calling_$sample --ntasks=1 --cpus-per-task=$nProcesses --mem-per-cpu=8000 --partition=mediumq --time=2-00:00:00 -e "$logdir/meth_calling_${sample}_%j.err" -o "$logdir/meth_calling_${sample}_%j.log" $scripts/getMeth_deduced.sh $working_dir $unmapped_fastq $ref_genome_fasta $genome_id $sample $samtools_path $bsmap_path $biseq_path $nProcesses $nonCpG $scripts $bedtools_path
			((submitted++))
		else
			get_proc_stats "$scripts/getMeth_deduced.sh $working_dir $unmapped_fastq $ref_genome_fasta $genome_id $sample $samtools_path $bsmap_path $biseq_path $nProcesses $nonCpG $scripts $bedtools_path &> $logdir/meth_calling_${sample}.log" "$step"
		fi
	else
		echo "$sample already processed. Not submitted!"
	fi
	((count++))
done

if [ ! $submitted = 0 ]; then
	wait_for_slurm $wait_time $submitted $working_dir
fi

if [ ! `ls $working_dir/$genome_id/*/biseqMethcalling/*cpgMethylation*.bed*|wc -l` = $count ]; then
	echo "Didn't find the expected number of cpgMethylation.bed files ($count). Exiting!"
	exit 1
fi

fail=0
if [ $nonCpG = "TRUE" ]; then
	if [ ! `ls $working_dir/$genome_id/*/biseqMethcalling/*cphpgMethylation*.bed*|wc -l` = $count ]; then
		echo "Didn't find the expected number of cphpgMethylation.bed files ($count). Exiting!"
		fail=1
	fi
	if [ ! `ls $working_dir/$genome_id/*/biseqMethcalling/*cphphMethylation*.bed*|wc -l` = $count ]; then
		echo "Didn't find the expected number of cphphMethylation.bed files ($count). Exiting!"
		fail=1
	fi
	if [ $fail = 1 ]; then
		exit 1
	fi
fi

#---------------------READ_Mapping_and_Methcalling_END--------------------------

#------------------------------CALCULATE_PDR------------------------------------

genome_dir=$working_dir/reduced/consensus/ 

step="\n-------PDR calculation-------\n"
printf "$step"
	
for aligned_bam in `ls $working_dir/$genome_id/*/*.bam`; do
	run_sample=$(basename $aligned_bam .bam)
	sample=${run_sample//.all.aligned*/}	
	echo $run_sample
	PDR_bed=$working_dir/$genome_id/$sample/${run_sample}_PDR.bed
	echo $PDR_bed
	if [ ! -s $PDR_bed ]; then
		if [ $parallel = "TRUE" ]; then
			sbatch --export=ALL --get-user-env --job-name=PDR_$run_sample --ntasks=1 --cpus-per-task=1 --mem-per-cpu=4000 --partition=shortq --time=08:00:00 -e "$logdir/PDR_${run_sample}_%j.err" -o "$logdir/PDR_${run_sample}_%j.log" $scripts/bisulfiteReadConcordanceAnalysis.py --infile=$aligned_bam  --outfile=$PDR_bed --skipHeaderLines=0 --genome=$genome_id --genomeDir=$genome_dir --samtools=$samtools_path
		else
			get_proc_stats "$scripts/bisulfiteReadConcordanceAnalysis.py --infile=$aligned_bam  --outfile=$PDR_bed --skipHeaderLines=0 --genome=$genome_id --genomeDir=$genome_dir --samtools=$samtools_path &> $logdir/PDR_${run_sample}.log" "$step"
		fi
	else
		echo "PDR bed file already exists. Skipping!"
	fi
	while [ ! -s $genome_dir/$genome_id/$genome_id.idx ]; do
		sleep 0.1m
	done
done



#------------------------------CALCULATE_PDR_END--------------------------------


#---------------------DFFERENTIAL_METH_ANALYSIS_START--------------------------
count=0
submitted=0
rm $working_dir/*.done 2>/dev/null
step="\n-------Differential methylation analysis for gpg motifs-------\n"
printf "$step"
((count++))
motif="cpg"
if [ ! -f $working_dir/$genome_id/diffMeth_$motif/*_diff_meth.tsv ]; then
	if [ $parallel = "TRUE" ]; then
		sbatch --export=ALL --get-user-env --job-name=diffMeth --ntasks=1 --cpus-per-task=1 --mem-per-cpu=50000 --partition=shortq --time=12:00:00 -e "$logdir/diffMeth_${motif}_%j.err" -o "$logdir/diffMeth_${motif}_%j.log" $scripts/diffMeth.R $working_dir $genome_id $species $genome_id $sample_annotation $compCol $groupsCol $nTopDiffMeth $scripts $motif $unconv_tag
		((submitted++))
	else
		get_proc_stats "$scripts/diffMeth.R $working_dir $genome_id $species $genome_id $sample_annotation $compCol $groupsCol $nTopDiffMeth $scripts $motif $unconv_tag &> $logdir/diffMeth_$motif.log" "$step"
	fi
else
	echo "$motif diff_meth.tsv already exists. Skipping..."
fi

if [ $nonCpG = "TRUE" ]; then

	step="\n-------Differential methylation analysis for gphpg motifs-------\n"
	printf "$step"
	((count++))
	motif="cphpg"
	if [ ! -f $working_dir/$genome_id/diffMeth_$motif/*_diff_meth.tsv ]; then
		if [ $parallel = "TRUE" ]; then
			sbatch --export=ALL --get-user-env --job-name=diffMeth --ntasks=1 --cpus-per-task=1 --mem-per-cpu=50000 --partition=shortq --time=12:00:00 -e "$logdir/diffMeth_${motif}_%j.err" -o "$logdir/diffMeth_${motif}_%j.log" $scripts/diffMeth.R $working_dir $genome_id $species $genome_id $sample_annotation $compCol $groupsCol $nTopDiffMeth $scripts $motif $unconv_tag
			((submitted++))		
		else
			get_proc_stats "$scripts/diffMeth.R $working_dir $genome_id $species $genome_id $sample_annotation $compCol $groupsCol $nTopDiffMeth $scripts $motif $unconv_tag &> $logdir/diffMeth_$motif.log" "$step"
		fi
	else
		echo "$motif diff_meth.tsv already exists. Skipping..."
	fi

	step="\n-------Differential methylation analysis for gphph motifs-------\n"
	printf "$step"
	((count++))
	motif="cphph"
	if [ ! -f $working_dir/$genome_id/diffMeth_$motif/*_diff_meth.tsv ]; then
		if [ $parallel = "TRUE" ]; then
			sbatch --export=ALL --get-user-env --job-name=diffMeth --ntasks=1 --cpus-per-task=1 --mem-per-cpu=50000 --partition=shortq --time=12:00:00 -e "$logdir/diffMeth_${motif}_%j.err" -o "$logdir/diffMeth_${motif}_%j.log" $scripts/diffMeth.R $working_dir $genome_id $species $genome_id $sample_annotation $compCol $groupsCol $nTopDiffMeth $scripts $motif $unconv_tag
			((submitted++))	
		else
			get_proc_stats "$scripts/diffMeth.R $working_dir $genome_id $species $genome_id $sample_annotation $compCol $groupsCol $nTopDiffMeth $scripts $motif $unconv_tag &> $logdir/diffMeth_$motif.log" "$step"
		fi
	else
		echo "$motif diff_meth.tsv already exists. Skipping..."
	fi
fi


if [ ! $submitted = 0 ]; then
	wait_for_slurm $wait_time $submitted $working_dir
	fi
rm $working_dir/*.done 2>/dev/null

if [ ! `ls $working_dir/$genome_id/diffMeth_*/*_diff_meth.tsv|wc -l` = $count ]; then
	echo "Differential methylation analysis failed. Didn't find the expected number of diff_meth.tsv files ($count). Exiting!"
	exit 1
fi



#---------------------DFFERENTIAL_METH_ANALYSIS_END--------------------------
echo "Done."
