#!/inside/home/common/bin/python2.7 

import sys, os
import argparse
from Bio import Phylo
import gzip
from cStringIO import StringIO
import numpy as np
import history_segments_module as histseg

bindir="/inside/home/tballing/cnavg-study/cnavgmbin"

def get_data_for_region(chr, start, end, cnavgdir, sim):
	mysegs=[]
	myhistories={} #keys: a historyid, values: a list of event IDs
	braneyfn="%s/HISTORIES_%d.braney" % (cnavgdir, sim)
	braneyf=gzip.open(braneyfn, 'rb')
	# get information from the trees
	treefn="%s/HISTORY_TREES_%d" % (cnavgdir, sim)
	treef=open(treefn, 'r')
	treenum=0
	treeline=treef.readline()
	mytree=Phylo.read(StringIO(treeline), 'newick')
	mytmpsegs=[]  # a holder for overlapping segments all from the same simulated history step
	prevhistoryid=0
	for braneyline in braneyf:
		data=braneyline.strip().split('\t')
		if len(data) > 1 and data[0] != 'A': # skip all adjacencies
 			# check if the segment overlaps our region of interest
			braneyseg=histseg.Braney_seg(braneyline)
			if (braneyseg.chr == chr and braneyseg.start < int(end) and braneyseg.end>int(start)):
				while treenum < braneyseg.historyid:
					treeline=treef.readline()
					mytree=Phylo.read(StringIO(treeline), 'newick')
					treenum = treenum + 1
				if treenum == braneyseg.historyid:
					get_prevalence_from_tree(mytree, braneyseg)
				else: # should only get here when treenum> historyid, which shouldn't happen
					sys.stderr.write("ERROR: didn't find tree for %d\n" % (braneyseg.historyid))
				if braneyseg.historyid > prevhistoryid:
					firstseg=get_first_event(mytmpsegs)
					mysegs.append(firstseg)
					mytmpsegs=[]
					prevhistoryid=braneyseg.historyid
#				sys.stderr.write("afmytmpsegs: %d, braneyseg: %s, phid %d\n" % (len(mytmpsegs), braneyseg, prevhistoryid))
				mytmpsegs.append(braneyseg)
	firstseg=get_first_event(mytmpsegs)
	mysegs.append(firstseg)	
	braneyf.close()
	treef.close()
	return(mysegs)
	
def get_first_event(mytmpsegs):
	mytmpsegs.sort(key=lambda x: x.preval, reverse=True)
	keepseg=None
	# want the event that has a cn change > 1
	for seg in mytmpsegs: 
		if abs(seg.cnval/seg.preval) > 1.5:
			keepseg= seg
			break
	if not keepseg: 
		keepseg=mytmpsegs[0]
	return keepseg

def get_prevalence_from_tree(mytree, seg):
	for clade in mytree.get_terminals():
		if clade.name:
			treeptr=clade.name
		else:
			treeptr=str("%d" % clade.confidence)
		if treeptr == seg.ptrid:
			seg.preval=clade.branch_length

def get_prev_and_cn_stats(seglist):
	# create arrays of cn changes and prevalences so you can get stats on them 
	cnchanges= np.empty(len(seglist), dtype=float)
	prevals = np.empty(len(seglist), dtype=float)
	weights = np.empty(len(seglist), dtype=float)
	for i in xrange(len(seglist)): 
		seg=seglist[i]
		cnchanges[i]=seg.cnval/seg.preval
		prevals[i]=seg.preval
		weights[i]=1/seg.complexity
	return (np.average(cnchanges, weights=weights), np.std(cnchanges), np.average(prevals, weights=weights), np.std(prevals), cnchanges.size)

def read_in_dat(datfn):
	mysegs=[]
	datf=open(datfn, 'r')
	for line in datf:
		data=line.strip().split('\t')
		braneyseg=histseg.Braney_seg("%s\t%s\t%s\t%f\t%s\t0\t0\t0\t0\t%s\t0\n" % (data[0], data[1], data[2], -1* float(data[3]), data[6], data[5]))
		braneyseg.preval=float(data[4])
		mysegs.append(braneyseg)
	datf.close()
	return mysegs

parser = argparse.ArgumentParser(description='Will output a distribution of CN change and prevalence for the earliest CN change in a genomic region')
parser.add_argument('bed', help='a bedfile of regions to gets distributions for')
parser.add_argument('cnavg', help='the CN-AVG output directory of a sample')
parser.add_argument('--outdir', '-o', help='a directory to put .dat files with all of the datapoints for each region in the bedfile. If the .dat file for a region already exists in the directory, it will be read in unless the --overwrite option is used.')
parser.add_argument('--overwrite', help='Will overwrite .dat files')
args=parser.parse_args()
cnavgdir=args.cnavg


if args.outdir and not os.path.exists(args.outdir):
	os.makedirs(args.outdir)
	
for bedline in open(args.bed, 'r'):
	(chr, start, end)=bedline.strip().split('\t')
	datfn="%s/%s_%s-%s.dat" % (args.outdir, chr, start, end)	
	if args.outdir:
		if not os.path.exists(args.outdir):
			os.makedirs(args.outdir)
		if not os.path.exists(datfn) or args.overwrite:
			outdat=open(datfn, 'w')
			for sim in xrange(10): # get data from each of the simulations 
				sys.stderr.write("working on sim: %d\n" % (sim))
				braneysegs=get_data_for_region(chr, start, end, args.cnavg, sim)
				for seg in braneysegs:
					outdat.write("%s\t%d\n" % (str(seg), sim))
			outdat.close()
		else: 
			braneysegs=read_in_dat(datfn)
			stats=get_prev_and_cn_stats(braneysegs)
			sys.stdout.write("%s\t%s\t%s\t%s\n" % (chr, start, end, "\t".join(map(str, stats))))
		 
