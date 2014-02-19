#!/inside/home/common/bin/python2.7
import sys, os
import gzip 
import glob 
import re
import history_segments_module as histseg
bindir="/inside/home/tballing/cnavg-study/cnavgmbin"

def score_segment_and_break_histories(sampledir):
	allsegment_hists=[]
	totalprob=0
	braneyfiles=glob.glob(sampledir+"/"+"HISTORIES_?.braney")
	for braneyfn in braneyfiles:
		hi=int(re.match(".*HISTORIES_(\d+)\.braney", braneyfn).group(1))
		totalprob+=histseg.get_total_history_prob(braneyfn)
		unisegments=combine_seghists_across_histories(braneyfile)
		sys.stderr.write("unisegments is %d\n" % (len(unisegments)))
		# change the history id to keep track of which simulation they came from 
		for seg in unisegments: 
			newhids= [x+(hi*10000) for x in seg.histories]
			seg.histories = newhids
#			sys.stderr.write("uniseg%d\t%s" % (hi, seg))
		allsegment_hists+=unisegments
		tmpunisegs= get_unisegments(sorted(allsegment_hists, key=lambda x: (x.chr, x.start, x.end)))
		allsegment_hists=tmpunisegs
		sys.stderr.write("allsegment_hists is %d\n" % (len(allsegment_hists)))
#		for seg in allsegment_hists: 
#			sys.stderr.write("allseg%d:\t%s" % (hi, seg))
	for seghist in allsegment_hists: 
		seghist.likelihood=histseg.compute_likelihood(seghist.costs, totalprob)
		sys.stdout.write(str(seghist))

def combine_seghists_across_histories(braneyfile): 
	numhistory=0
	universal_seghists=[]
	seglines=[]
	braneyf=gzip.open(braneyfile, 'rb')
	bline=braneyf.readline()
	while len(bline) >0: 
		if re.match("chr", bline): 
			braneyseg=histseg.Braney_seg(bline)
			myseghist=histseg.Segment_history(braneyseg)
			if myseghist.histories[0] == numhistory:
				seglines.append(myseghist)
			else: # moved on to a new history and process the segments of the previous history
				if len(seglines) >0: 
					segshist=get_unisegments(sorted(seglines, key=lambda x: (x.chr, x.start, x.end)))
#					sys.stderr.write("len of segshist: %d\n" % (len(segshist)))
#					for seg in segshist:
#						sys.stderr.write("segshist%d\t%s" % (numhistory, seg))
					universal_seghists += segshist
					tmp_uniseghists=get_unisegments(sorted(universal_seghists, key=lambda x: (x.chr, x.start, x.end)))
					universal_seghists=tmp_uniseghists
#					sys.stderr.write("len of universal_seghists: %d\n" % (len(universal_seghists)))
#					for seg in universal_seghists:
#						sys.stderr.write("uniseg%d\t%s" % (numhistory, seg))
				numhistory=myseghist.histories[0]
				seglines=[]
		bline=braneyf.readline()
	# process the last history 
	segshist=get_unisegments(sorted(seglines, key=lambda x: (x.chr, x.start, x.end)))
	universal_seghists += segshist
	tmp_uniseghists=get_unisegments(sorted(universal_seghists, key=lambda x: (x.chr, x.start, x.end)))
	universal_seghists=tmp_uniseghists
	return(universal_seghists)

def get_unisegments(seghists): 
	unisegments=[]
	current_seghist=seghists[0]
	working_seghists=[current_seghist]
	for myseg in seghists[1:]:
#		sys.stderr.write("myseg: %s" % (myseg))
		working_tmphists=[]
		for workseg in working_seghists:
#			sys.stderr.write("workseg: %s" % (workseg))
			if workseg.overlaps(myseg):
				new_segments=workseg.split_and_merge_in(myseg)
				working_tmphists+=new_segments
			elif myseg.comes_before(workseg):
				working_tmphists.append(workseg)
			elif workseg.comes_before(myseg): 
				unisegments.append(workseg)
		workseg=working_seghists[-1]
		if workseg.comes_before(myseg):
			working_tmphists.append(myseg)
		if len(working_tmphists) >1: 
			working_seghists=remove_duplicate_histsegs(working_tmphists)
		else: 
			working_seghists=working_tmphists
#		for seg in working_seghists:
#			sys.stderr.write("tmpseg: %s" % (seg))
	# print out all of the rest of the segments at the end
	unisegments+=working_seghists 
	return unisegments

def remove_duplicate_histsegs(histseglist):
	uniqsegs=[]
	tmpseglist=[]
	while len(histseglist) >1: 
		sega=histseglist[0]
		for j in xrange(1, len(histseglist)): 
			segb=histseglist[j]
#			sys.stderr.write("checking sega: %s against segb: %s" % (sega, segb))
			if sega == segb: 
#				sys.stderr.write("sega == segb\n") 
				if segb.numhists > sega.numhists: 
					sega=segb #keep the seghist with the most number of histories
			elif sega.mergedwith(segb):
#				sys.stderr.write("sega mergedwith segb\n") 
				if len(segb.cnvals) > len(sega.cnvals): 
					sega=segb
			else: 
				tmpseglist.append(segb)
		histseglist=tmpseglist
		tmpseglist=[]
		uniqsegs.append(sega)
	return uniqsegs + histseglist
				
	
if __name__ == "__main__": 
	import argparse
	parser = argparse.ArgumentParser(description='Given a cn-avg output directory, it will merge the last step of multiple simulations and create file of universal segments with lists of CNvalues and prevalences across the different histories')
	parser.add_argument('sample', help='a cn-avg output directory.')
	args=parser.parse_args()
	score_segment_and_break_histories(args.sample)
#	braneyfile="/inside/home/dzerbino/cn-avg/runs/gbm2/TEST/HISTORIES_9.braney"
#	unisegments=combine_seghists_across_histories(braneyfile)
#	sys.stderr.write("end: len unisegments: %d\n" % (len(unisegments)))
#	for seg in unisegments: 
#		sys.stdout.write(str(seg))

