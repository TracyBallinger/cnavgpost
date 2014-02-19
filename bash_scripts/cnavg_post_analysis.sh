#!/bin/sh 

# cnavg_post_analysis.sh written by Tracy Ballinger, last edited 10/20/2013

if [ $# -ne 2 ]; then 
	echo "usage: cnavg_post_analysis.sh <cnavg_output directory><samplename>"
	exit 65
fi 
bindir=/inside/home/tballing/cnavg-study/cnavgmbin
tabixfile=/inside/home/tballing/cnavg-study/cancer-gene.unique.coordinates.bed.bgz
bedfile=/inside/home/tballing/cnavg-study/cancer-gene.merged.bed
outpost=/inside/home/tballing/grotto/cnavg_outpost

cnavgout=$1
samplename=$2
mkdir -p $outpost/$samplename
cd $outpost/$samplename

if [ ! -e $samplename.pevnts ]; then 
	echo "($bindir/score_and_link_cycles.py --cnavg $cnavgout --outpickle $samplename.pevnts --links $samplename.links ) &> event_link_error" 
	($bindir/score_and_link_cycles.py --cnavg $cnavgout --outpickle $samplename.pevnts --links $samplename.links ) &> event_link_error  
fi 

if [ ! -e $samplename.pedgs ]; then 
	echo "$bindir/score_edges_within_pevents.py --outpickle $samplename.pedgs $samplename.pevnts"
	$bindir/score_edges_within_pevents.py --outpickle $samplename.pedgs $samplename.pevnts
fi 

if [ ! -e $samplename.ann ]; then 
	echo "$bindir/annotate_events.py $samplename.pevnts $tabixfile > $samplename.ann"	
	$bindir/annotate_events.py $samplename.pevnts $tabixfile > $samplename.ann 	
fi

if [ ! -e $samplename.gnrank ]; then 
	echo "$bindir/histories_to_gene_orders.py $samplename.pevnts $samplename.ann > $samplename.gnrank"	
	$bindir/histories_to_gene_orders.py $samplename.pevnts $samplename.ann > $samplename.gnrank 	
fi

if [ ! -e $samplename.gnpr ]; then 
	echo "skipping $bindir/likelihood_score_gene_pairs.py $samplename.pevnts $samplename.ann > $samplename.gnpr"	
#	$bindir/likelihood_score_gene_pairs.py $samplename.pevnts $samplename.ann > $samplename.gnpr 	
fi

