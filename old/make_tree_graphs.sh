#!/bin/sh 

# make_tree_graphs.sh written by Tracy Ballinger, last edited 4/14/2013

if [ $# -ne 3 ]; then 
	echo "usage: make_tree_graphs.sh HISTORY_TREES_0 HISTORIES_0.braney output.pdf" 
	exit 65; 
fi 

history_trees=$1
if [ ! -e $history_trees ]; then 
	echo "$history_trees doesn't exist"
	exit 65; 
fi 

history_braney=$2
if [ ! -e $history_braney ]; then 
	echo "$history_braney doesn't exist"
	exit 65; 
fi 
output=$3

treenum=1
head -$treenum $history_trees > tmp.tre
gzip -dc $history_braney | grep ^A | cut -f8-11,13-15 > tmp
gzip -dc $history_braney | grep -v ^A | cut -f4-7,9-11 | cat tmp - | sed 's/-//' | sort -u > tmpdata.txt
awk '$2==(treenum-1)' treenum=$treenum tmp | sed 's/-//' | sort -u > tmpdata.txt 

R --vanilla --args tmp.tre tmpdata.txt $output < tree_plots.R

# rm tmp tmpdata.txt tmp.tre 
	
