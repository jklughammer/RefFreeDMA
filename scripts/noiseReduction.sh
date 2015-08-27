#! /bin/bash

#Author: Johanna Klughammer
#Date: 26.07.2015


working_dir=$1
filtLim=$2
maxSamples=$3


cd $working_dir/reduced/

filesCount=0
totalCount=0
for file in `ls *_uniq.ref`; do
	((filesCount++))
	((totalCount++))
	echo $file
	cat $file >> merged.ref
	if [ $filesCount = $filtLim ]; then
		echo "reducing..."
		sort -k 3,3 merged.ref | uniq -f 2 -D >> merged_dupl.ref
		rm merged.ref
		filesCount=0
	fi
	if [ $totalCount = $maxSamples ]; then
		break
	fi
done 

if [ -f merged.ref ]; then
	allFiles=(`ls *_uniq.ref`)
	missing=`echo "$filtLim-$filesCount"|bc`
	for i in `seq 1 $missing`; do
		echo ${allFiles[i]}
		cat ${allFiles[i]} >> merged.ref
	done
	sort -k 3,3 merged.ref | uniq -f 2 -D >> merged_dupl.ref
	rm merged.ref
	filesCount=0
fi
sort -k 2,2 merged_dupl.ref | uniq -f 1| sort -k 3,3 > merged_uniq.ref

cd -

touch $working_dir/noiseReduction.done