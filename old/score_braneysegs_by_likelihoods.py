#!/inside/home/common/bin/python2.7
import sys, os
import gzip 
import glob 
import re
import combined_segments_module as histseg
bindir="/inside/home/tballing/cnavg-study/cnavgmbin"

def score_braneysegs_by_likelihood(sampledir):
	allsegments=[]
	totalprob=0
	braneyfiles=glob.glob(sampledir+"/"+"HISTORIES_?.braney")
	for braneyfn in braneyfiles:
		hi=int(re.match(".*HISTORIES_(\d+)\.braney", braneyfn).group(1))
		totalprob+=histseg.get_total_history_prob(braneyfn)
		sys.stderr.write("totalprob is %s\n" % (str(totalprob)))
		unique_segs=combine_segs_across_histories(braneyfn)
		sys.stderr.write("number of segs: %d\n" % (len(unique_segs)))
		# change the history id to keep track of which simulation they came from 
		for seg in unique_segs: 
			newhids= [x+(hi*10000) for x in seg.histories]
			seg.histories = newhids
		allsegments+=unique_segs
		tmpunisegs= get_unique_segments(sorted(allsegments, key=lambda x: (x.chr, x.start, x.end)))
		allsegments=tmpunisegs
	for seg in allsegments: 
		seg.likelihood=histseg.compute_likelihood(seg.costs, totalprob)
		sys.stdout.write(str(seg))

def combine_segs_across_histories(braneyfile): 
	numhistory=0
	unique_segs=[]
	seglines=[]
	braneyf=gzip.open(braneyfile, 'rb')
	bline=braneyf.readline()
	while len(bline) >0: 
		braneyseg=histseg.Braney_seg(bline)
#		braneyseg.order_ends()
		myseg=histseg.Combined_braneyseg(braneyseg)
		if myseg.histories[0] == numhistory:
			seglines.append(myseg)
		else: # moved on to a new history and process the segments of the previous history
			if len(seglines) >0: 
				unique_segs += seglines
				if len(unique_segs)>5000:	
					tmp_unisegs=get_unique_segments(sorted(unique_segs, key=lambda x: (x.chr, x.chr2, x.start, x.end)))
					unique_segs=tmp_unisegs
			numhistory=myseg.histories[0]
#			sys.stderr.write("working on history %d, there are %d unique_segs\n" % (numhistory, len(unique_segs)))
			seglines=[]
		bline=braneyf.readline()
	# process the last history 
	segslist=get_unique_segments(sorted(seglines, key=lambda x: (x.chr, x.chr2, x.start, x.end)))
	unique_segs += segslist
	tmp_unisegs=get_unique_segments(sorted(unique_segs, key=lambda x: (x.chr, x.chr2, x.start, x.end)))
	unique_segs=tmp_unisegs
	return(unique_segs)

def get_unique_segments(seglist): 
	unique_segments=[]
	current_seg=seglist[0]
	working_segs=[current_seg]
	for myseg in seglist[1:]:
		beenadded=False
		working_tmpsegs=[]
		for workseg in working_segs:
			if workseg == myseg and (abs(workseg.prevals[0]-myseg.prevals[0]) < prevalence_error):
				if myseg.histories[0] not in workseg.histories:
					workseg.addin(myseg)
				working_tmpsegs.append(workseg)
				beenadded=True
			elif workseg.comes_before(myseg): 
				unique_segments.append(workseg)
			elif myseg.comes_before(workseg):
				working_tmpsegs.append(workseg)
			else: # the segments overlap or myseg comes before workseg 
				working_tmpsegs.append(workseg)	
		if not beenadded:
			working_tmpsegs.append(myseg)	
		working_segs=working_tmpsegs
	# print out all of the rest of the segments at the end
	unique_segments+=working_segs 
	return unique_segments

	
if __name__ == "__main__": 
	import argparse
	parser = argparse.ArgumentParser(description='Given a cn-avg output directory, it will merge the .braney files into one, combining equivalent adjacencies and segments and scoring them by their likelihood.') 
	parser.add_argument('sample', help='a cn-avg output directory.')
	parser.add_argument('--prevalence_error', help='the difference in prevalences to be considered the same.', type=float, default=0.05)
	args=parser.parse_args()
	global prevalence_error
	prevalence_error=args.prevalence_error
	score_braneysegs_by_likelihood(args.sample)
#	braneyfile="/inside/home/dzerbino/cn-avg/runs/gbm2/TEST/HISTORIES_0.braney"
#	braneyfile="/inside/home/dzerbino/cn-avg/runs/gbm/TCGA-06-0145/HISTORIES_0.braney"
#	unisegments=combine_segs_across_histories(braneyfile)
#	sys.stderr.write("end: len unisegments: %d\n" % (len(unisegments)))
#	for seg in unisegments: 
#		sys.stdout.write(str(seg))

