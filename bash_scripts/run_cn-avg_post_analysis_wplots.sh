#!/bin/sh 

# run_cn-avg_post_analysis_wplots.sh written by Tracy Ballinger, last edited 8/28/2013

if [ $# -ne 2 ]; then 
	echo "usage: run_cn-avg_post_analysis_wplots.sh cnavg_dir sample_name"
	echo "a directory called sample_name.out2 will be created that contains analysis files."
	exit 65; 
fi 

bindir=/inside/home/tballing/cnavg-study/cnavgmbin
# need to change this bedfile to consolidate the coordinates...
bedfile=/inside/home/tballing/cnavg-study/cancer-gene.merged.bed
tabixfile=/inside/home/tballing/cnavg-study/cancer-gene.merged.bed.bgz
cnavg_dir=$1
if [ ! -e $cnavg_dir ]; then
    echo "$cnavg_dir doesn't exist"
    exit 65;
fi

samplename=$2

outputdir=$samplename.out2
if [ ! -e $outputdir ]; then
    mkdir $outputdir
    exit 65;
fi
odir=$outputdir
cd $odir

# Get universal segments across the different history simulations
$bindir/create_segment_history_file.py $cnavg_dir > segment_histories.txt 
$bindir/get_best_segment_history.py segment_histories.txt > segment_historiesb.txt 
rm segment_histories.txt

# Do a per gene analysis to see where the first change was 
$bindir/get_gene_firstCNchange_stats.py $bedfile $cnavg_dir --outdir datfiles > gene_cnprev.stats 

# Make a file of the rearrangement cycles
$bindir/run_braney_to_score_cycles.py $cnavg_dir > $samplename.evnts
$bindir/paint_cycles_with_genes.py $samplename.evnts $tabixfile > $samplename.evntswgenes

# Score gene pairs 
$bindir/score_gene_pairs.py $bedfile $cnavg_dir > $samplename.gene_pairs.dat
