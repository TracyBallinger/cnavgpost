#!/bin/sh 

# run_analysis_on_cnavg_output.sh written by Tracy Ballinger, last edited 10/20/2013

if [ $# -ne 1 ]; then 
	echo "usage: run_analysis_on_cnavg_output.sh <list of cnavg_output directories>"
	exit 65
fi 
bindir=/inside/home/tballing/cnavg-study/cnavgmbin
tabixfile=/inside/home/tballing/cnavg-study/cancer-gene.unique.coordinates.bed.bgz
bedfile=/inside/home/tballing/cnavg-study/cancer-gene.merged.bed
outpost=/inside/home/tballing/grotto/cnavg_outpost
samplelist=$1

cat $samplelist | while read line; do
	set `echo $line`
	cnavgout=$1
	samplename=$2
	mkdir -p $outpost/$samplename
	cd $outpost/$samplename
	if [ ! -e $samplename.evnts ]; then 
		($bindir/score_and_link_cycles.py --cnavg $cnavgout --outpickle $samplename.pevnts --events $samplename.evnts --links $samplename.links ) &> event_link_error  
	fi 
	if [ ! -e $samplename.evntswgns ]; then 
		($bindir/paint_cycles_with_genes.py $samplename.evnts $tabixfile > $samplename.evntswgns) &> errtbx 	
	fi
	cd $outpost
	#$bindir/score_gene_pairs.py $bedfile $cnavgout > $samplename.gnprs.dat
done 
