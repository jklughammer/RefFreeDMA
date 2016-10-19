#! /bin/bash

#Author: Johanna Klughammer
#Date: 26.07.2015


working_dir=$1
filtLim=$2
maxSamples=$3
unconv_tag=$4


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
	echo "filling up to filtLim"
	allFiles=(`ls *_uniq.ref`)
	missing=`echo "$filtLim-$filesCount-1"|bc`
	for i in `seq 0 $missing`; do
		echo ${allFiles[i]}
		cat ${allFiles[i]} >> merged.ref
	done
	sort -k 3,3 merged.ref | uniq -f 2 -D >> merged_dupl.ref
	rm merged.ref
	filesCount=0
fi



if [ -f *_uniq.ref.unconv ];then

	sort -k 2,2 merged_dupl.ref | uniq -f 1 > merged_uniq.temp
	cat <(cat *_uniq.ref.unconv| uniq -f 1) <(cat merged_uniq.temp)| sort -k 3,3 | uniq -f 2 -D |grep $unconv_tag >> merged_uniq.temp

	sort -k 3,3 merged_uniq.temp > merged_uniq.ref
	rm merged_uniq.temp

else

	sort -k 2,2 merged_dupl.ref | uniq -f 1| sort -k 3,3 > merged_uniq.ref

fi



cd -

touch $working_dir/noiseReduction.done
