#!/usr/bin/env python 

import argparse
import sys, os
import cPickle as pickle
import numpy as np
from cnavgpost.genehistory.segments_history_module import * 
import cnavgpost.mergehistories.event_cycles_module as ecycles

def create_seghists_from_eventsegs(esegs, histScores): #Note that the edges must be sorted by segstr.  
	totalp=ecycles.compute_totalp(histScores)
	myseghists=[]
	working_seghists=[SegmentHistory(esegs[0])]
	for e in esegs[1:]:
#		sys.stderr.write("%sworking_seghists: %d\tmyseghists: %d\n" % (str(e), len(working_seghists), len(myseghists)))
		tmpseghists=[]
		for seghist in working_seghists: 
			if seghist.overlaps(e): 
				new_seghists=merge_in_edge(e, seghist)
#				sys.stderr.write("overlaps, %d\n" % len(new_seghists))
				tmpseghists+=new_seghists
			elif seghist.comes_after(e): 
			 	tmpseghists.append(seghist)
			elif seghist.comes_before(e): 
				myseghists.append(seghist)
		lastseg=working_seghists[-1]
		if lastseg.overlaps(e) and e.end > lastseg.end: 
			newseg=SegmentHistory(e)
			newseg.start=working_seghists[-1].end+1
			tmpseghists.append(newseg)
		elif lastseg.comes_before(e): 
			tmpseghists.append(SegmentHistory(e))
		working_seghists=tmpseghists
	myseghists+=working_seghists
	return myseghists

#	sys.stderr.write("done making seghists, now creating the CN profiles.\n")
#	pickle.dump(myseghists, open("seghists.pck", 'w'), pickle.HIGHEST_PROTOCOL)
#	for sgh in myseghists: 
#		create_CNprofiles_from_Edges(sgh, histScores, totalp)
#	sys.stderr.write("done making CN profiles\n")
#	return(myseghists)

def make_esegs_from_events(events): 
	breakpointcn=-10
	myesegs=[]
	for e in events:
		e.unpack()
		for s in e.segs: 
			if s.seg: # if it's a genomic segment, not an adjacency
				myesegs.append(EventSegment(e, s.chr, s.start, s.end))
			else: # if it's an adjacency, treat each breakpoint like a mini deletion
				e1=EventSegment(e, s.chr, s.start, s.start)
				e1.cnval=breakpointcn
				e2=EventSegment(e, s.chr2, s.end, s.end)
				e2.cnval=breakpointcn
				myesegs+=[e1, e2]
	sortedesegs=sorted(myesegs, key=lambda x: (x.chr, x.start, x.end))
	return sortedesegs	

if __name__ == "__main__": 
	import argparse
	parser = argparse.ArgumentParser(description='Will take the edges from events and split the genome at breakpoints, with histories for each segment.')
	parser.add_argument('pedges', help='a pickled files of event edges.')
	parser.add_argument('hstats', help='historystats.txt file')
	parser.add_argument('outf', help='name of the output file')
	args=parser.parse_args()
	histScores=np.loadtxt(args.hstats, dtype=int)
	events=pickle.load(open(args.pedges, 'rb'))
	myesegs=make_esegs_from_events(events)
	sys.stderr.write("loaded and sorted events\n")
	seghists=create_seghists_from_eventsegs(myesegs, histScores)
	fh=open(args.outf, 'w')
	for sgh in seghists: 
		fh.write("sgh\t"+str(sgh))
	fh.close()
	pickle.dump(seghists, open("seghists.pk", 'w'), pickle.HIGHEST_PROTOCOL)

	
