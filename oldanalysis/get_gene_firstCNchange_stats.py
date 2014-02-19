#!/inside/home/common/bin/python2.7 

import sys, os
import numpy as np
import history_segments_module as histseg
import pysam
bindir="/inside/home/tballing/cnavg-study/cnavgmbin"

def get_data_for_region(chr, start, end, tabixfn):
	mysegs=[]
	myhistories={} #keys: a historyid, values: a list of event IDs
	mytmpsegs=[]  # a holder for overlapping segments all from the same simulated history step
	prevhistoryid=-1
	tabixfile=pysam.Tabixfile(tabixfn)
	braneylines=[]
	try:
		braneylines=tabixfile.fetch(reference=chr, start=start, end=end)
	except ValueError:
		sys.stderr.write("Error in tabix fetch for %s, %s:%d-%d, %s\n" % (tabixfn, chr, start, end, str(ValueError)))
	bsegs=[]
	for bline in braneylines: 
		bseg=histseg.Braney_seg(bline) 
		bsegs.append(histseg.Braney_seg(bline))
	# sort bsegs by historyid
	if len(bsegs) >0 : 
		bsegs.sort(key=lambda x: x.historyid) 
		for braneyseg in bsegs:
			if braneyseg.historyid > prevhistoryid:
				if len(mytmpsegs)>0:
					firstseg=get_first_event(mytmpsegs)
					mysegs.append(firstseg)
				prevhistoryid=braneyseg.historyid
				mytmpsegs=[braneyseg]
			else: 
				mytmpsegs.append(braneyseg)			
		firstseg=get_first_event(mytmpsegs)
		mysegs.append(firstseg)	
	return(mysegs)
	
def get_first_event(mytmpsegs):
	mytmpsegs.sort(key=lambda x: x.preval, reverse=True)
	keepseg=None
	# want the event that has a cn change > 1
	for seg in mytmpsegs: 
		if abs(seg.cnval/seg.preval) >= 1:
			keepseg= seg
			break
	if not keepseg: 
		keepseg=mytmpsegs[0]
	return keepseg

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

def read_in_dat(datfh):
	mysegs=[]
	for line in datfh:
		data=line.strip().split('\t')
		braneyseg=histseg.Braney_seg("%s\t%s\t%s\t%f\t%s\t%s\t0\t0\t0\t0\t%s\t0\n" % (data[0], data[1], data[2], -1* float(data[3]), data[4], data[6], data[5]))
		mysegs.append(braneyseg)
	return mysegs

if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description='Will output a distribution of CN change and prevalence for the earliest CN change in a genomic region')
	parser.add_argument('bed', help='a bedfile of regions to gets distributions for')
	parser.add_argument('cnavg', help='the CN-AVG output directory of a sample')
	parser.add_argument('--outdir', '-o', help='a directory to put .dat files with all of the datapoints for each region in the bedfile. If the .dat file for a region already exists in the directory, it will be read in unless the --overwrite option is used.')
	parser.add_argument('--overwrite', help='Will overwrite .dat files', action='store_true', default=False)
	parser.add_argument('--tabixdir', '-t', help='a directory to put tabix files created from the .braney files.', default="./") 
	args=parser.parse_args()
	cnavgdir=args.cnavg

	if args.outdir and not os.path.exists(args.outdir):
		os.makedirs(args.outdir)
	
	# make a tabix file of the .braney so that it's easy to fetch overlapping regions 
	numsims=10
	tabixfiles=[]
	for sim in xrange(numsims): # get data from each of the simulations 
		braneyfn="%s/HISTORIES_%d.braney" % (args.cnavg, sim)
		(tabixsegs, tabixadj)=histseg.make_tabix_from_braney(braneyfn, args.tabixdir)
		tabixfiles.append((tabixsegs, tabixadj))
	
	#print out a header line to label the output columns 
	sys.stdout.write("chr\tstart\tend\tave_cnval\tstd_cnval\tave_prevalence\tstd_prevalence\tnum_histories\n")
	for bedline in open(args.bed, 'r'):
		data=bedline.strip().split('\t')
		chr=data[0]
		start=int(data[1])
		end=int(data[2])
		if args.outdir: 
			datfn="%s/%s_%d-%d.dat" % (args.outdir, chr, start, end)	
			if (not os.path.exists(datfn)) or (os.path.getsize(datfn) ==0) or args.overwrite:
				outdat=open(datfn, 'r+')
		else: 
			outdat=os.tmpfile()
		for sim in xrange(numsims): # get data from each of the simulations 
			(tabixsegs, tabixadj) = tabixfiles[sim]	
			braneysegs=get_data_for_region(chr, start, end, tabixsegs)
			for seg in braneysegs:
					outdat.write("%s\t%d\n" % (str(seg), sim))
		outdat.seek(0)
		braneysegs=read_in_dat(outdat)
		outdat.close()
		if len(braneysegs) >0:
			stats=get_prev_and_cn_stats(braneysegs)
			sys.stdout.write("%s\t%s\t%s\t%s\n" % (chr, start, end, "\t".join(map(str, stats))))
		else: 		
			stats=("NA", "NA", "NA", "NA", 0)
			sys.stdout.write("%s\t%s\t%s\t%s\n" % (chr, start, end, "\t".join(map(str, stats))))
		 
