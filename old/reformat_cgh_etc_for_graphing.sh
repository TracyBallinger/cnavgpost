#!/bin/sh 

if [ $# -ne 1 ]; then 
	echo "usage: reformat_cgh_etc_for_graphing.sh tracks_directory"
	exit 65; 
fi 

mydir=$1
if [ ! -e $mydir ]; then 
	echo "directory $mydir doesn't exist!" 
	exit 65; 
fi 

level2dir=$mydir/HMS__HG-CGH-244A/Level_2
for file in $level2dir/TCGA*.tsv; do 
	echo $file
	if [ ! -e $level2dir/tmp ]; then 
		cp $file $level2dir/tmp
	else
		cut -f2 $file | paste $level2dir/tmp - > $level2dir/tmp2
		mv $level2dir/tmp2 $level2dir/tmp 
	fi 
done 
mv $level2dir/tmp $level2dir/combined_sets.tsv
cat $mydir/HMS__HG-CGH-244A/Level_3/hms.harvard*.txt | grep -v barcode > $mydir/HMS__HG-CGH-244A/Level_3/all_copy_number_analysis.txt 

level2dir=$mydir/MSKCC__HG-CGH-244A/Level_2 
for file in $level2dir/MSK*.mat ; do 
    echo $file
    if [ ! -e $level2dir/tmp ]; then                         
        cp $file $level2dir/tmp
    else
        cut -f2 $file | paste $level2dir/tmp - > $level2dir/tmp2
        mv $level2dir/tmp2 $level2dir/tmp
    fi
done
mv $level2dir/tmp $level2dir/combined_sets.tsv
cat $mydir/MSKCC__HG-CGH-244A/Level_3/mskcc.org*.txt | grep -v barcode > $mydir/MSKCC__HG-CGH-244A/Level_3/all_copy_number_analysis.txt 

awk '{printf "%1.3f\n", $4}' $mydir/bambam.cov.bg | sort | uniq -c > $mydir/bambam.covhist.txt 
