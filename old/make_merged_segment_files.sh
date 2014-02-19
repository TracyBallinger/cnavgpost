#!/bin/sh 

# make_merged_segment_files.sh written by Tracy Ballinger, last edited 6/15/2013

if [ $# -ne 2 ]; then 
	echo "usage: make_merged_segment_files.sh cnavg_dir segmented_histories.txt" 
	echo "EX: make_merged_segment_files.sh /inside/home/dzerbino/cn-avg/runs/gbm/TCGA-06-0145 > segment_histories.txt "
	exit 65; 
fi 

workingdir=/inside/home/tballing/cnavg-study
cnavg_dir=$1
if [ ! -e $cnavg_dir ]; then 
	echo "$cnavg_dir doesn't exist"
	exit 65; 
fi
output=$2 
outputdir=/scratch/tjbtmp
if [ ! -e $outputdir ]; then 
	mkdir $outputdir
	exit 65; 
fi 

for treefile in `ls $cnavg_dir/HISTORY_TREES_?`; do
	echo treefile is $treefile
	x=`echo $treefile | awk '{split($1, a, "HISTORY_TREES_"); print a[2]}'`
	braneyfile=$cnavg_dir/HISTORIES_$x.braney
	$workingdir/get_history_tree_data.py $treefile $braneyfile --history_id=2500 | awk '{print $0"\t"x}' x=$x > $outputdir/history_$x.dat
done 

cat $outputdir/history_?.dat | sort -k1,1 -k2,2n -k3,3n -k8,8n > $outputdir/tmps.dat
$workingdir/get_segment_values_dist.py $outputdir/tmps.dat > $output 

rm -r $outputdir
