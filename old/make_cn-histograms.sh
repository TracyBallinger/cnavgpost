#!/bin/sh 

# make_cn-histograms.sh written by Tracy Ballinger, last edited 4/19/2014

if [ $# -ne 1 ]; then 
	echo "usage: make_cn-histograms.sh cactus_tree_outputdir graphs_outputdir"
	echo "This will create graphs_outputdir if it doesn't exist and put graphs"
	echo "of cactus tree output in there. "
	exit 65; 
fi 

sampleid=$1
cactusout=$1
if [ ! -e $cactusout ]; then 
	echo "$cactusout doesn't exist"
	exit 65; 
fi 
output=$2
if [ ! -e $output ]; then 
	mkdir $output
fi 

# preprep for some of the input
awk '{printf "%1.3f\n", $4}' $cactusout/tracks/bambam.cov.bg | sort | uniq -c > $output/bambam.covhist.txt
newbb=/inside/depot4/bambam/five3/gbm/$sampleid"_whole_cnv.txt"

R --vanilla --args $cactusout $output/bambam.covhist.txt $output/cn_histograms.pdf $newbb < cn_histograms.v2.R
