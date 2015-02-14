#!/usr/bin/env python 

import sys, os
import cnavgpost.mergehistories.event_cycles_module as histseg
import argparse
import cPickle as pickle 
import subprocess, re
import numpy as np

def filter_events_by_cn_overlap(events, bedfn, perc, cbs=True):
	if cbs: 
		mystr=subprocess.check_output("sed 1d %s | cut -f2-4" % bedfn, shell=True)
		bedlines=mystr.split("\n")
	else: 
		bedlines=open(bedfn, 'r').readlines()
	bedentries=[]
	for line in bedlines: 
		if line!="": 
			bedentries.append(BedEntry(line))
	bedentries=sorted(bedentries, key=lambda x: (x.chr, x.start, x.end))
	myevents=[]
	for event in events: 
		if event_overlaps_CN(event, bedentries, perc): 
			myevents.append(event) 
	return(myevents)

def event_overlaps_CN(event, bedentries, perc):
	if event.segs==[]:
		event.make_segs_from_str()
	keep=True
	for seg in event.segs: 
		if seg.seg:
			region=BedEntry("%s\t%d\t%d\n" % (seg.chr, seg.start, seg.end)) 
			poverlaps=percent_bed_overlaps(bedentries, region, perc)
#			sys.stderr.write("seg %s:%d-%d overlaps: %s\n" % (seg.chr, seg.start, seg.end, str(poverlaps)))
			if len(poverlaps) >0 and (max(poverlaps) < perc): 
				keep=False
	return keep

# find what the overlap is for the given region in the bedentries		
def percent_bed_overlaps(bedentries, region, perc):
	i=len(bedentries)/2
	bedentry=bedentries[i]
	overlap=bedentry.overlap(region)
	myvals=[]
	if overlap>0: 
		poverlap=min(float(overlap)/(bedentry.end-bedentry.start+1), 
					float(overlap)/(region.end-region.start+1))
		myvals.append(poverlap) 
	if bedentry.comes_before(region) and (i+1)< len(bedentries):
		myvals+=percent_bed_overlaps(bedentries[(i+1):], region, perc)
	elif i > 0: 
		myvals+=percent_bed_overlaps(bedentries[:i], region, perc)
	return myvals


class BedEntry: 
	def __init__(self, bedline):
		(chr, start, end) = bedline.strip().split('\t')[:4]
		self.chr=re.sub("chr", "", chr)
		self.start=int(start)
		self.end=int(end)
	
	def __str__(self): 
		return "%s:%d-%d" % (self.chr, self.start, self.end)

	def comes_before(self, other): 
		return ((self.chr < other.chr) or (self.chr == other.chr and self.start < other.start))		
		
	def comes_after(self, other): 
		return ((self.chr > other.chr) or (self.chr==other.chr and self.start>other.start))

	def overlap(self, region): 
		overlap=0
		if (self.chr == region.chr): 
#			sys.stderr.write("%s: %d-%d, %s: %d-%d\n" % (self.chr, self.start, self.end, region.chr, region.start, region.end))
			overlap=min(region.end, self.end) - max(region.start, self.start)
		return overlap 

## This will filter out events that cancel eachother out (ie, an amplification immediately followed by a deletion). 
def filter_fleeting_events(events, fcutoff=0): 
	for e in events: 
		(newstr, sign) = histseg.remove_signs_from_segstr(e.segstr)
		e.cnval=sign*e.cnval 
	sevents=sorted(events, key=lambda x: (x.segstr, x.prevalmean))
	e1=sevents[0]
	sameEs=[e1]
	filtered=[]
	for e2 in sevents[1:]: 
		if e1.segstr == e2.segstr:
			sameEs.append(e2)
		else: 
			filtered+=remove_canceling_histories(sameEs)	
			sameEs=[e2]
	return filtered

def remove_canceling_histories(events): 
	e1=events[0]
	for i in xrange(len(events)): 
		e1=events[i]
		for e2 in events[i:]: 
			if e1.segstr == e2.segstr and ((e1.cnval + e2.cnval) == 0): 
				# only keep the histories that are unique to each event. 
				histseg.cancel_Event_data(e1,e2)
	filtered=[]
	for e in events: 
		if len(e.histories) >0: 
			filtered.append(e)
	return filtered 


def main(args): 
	allevents=pickle.load(open(args.pevnt, 'rb'))
	sys.stderr.write("Begin, %d events\n" % (len(allevents)))
	for e in allevents:
		e.histories=histseg.listout_ranges(e.histRanges)
	if args.CNoverlap: 
		allevents=filter_events_by_cn_overlap(allevents, args.bed, args.perc, args.cbs)
		sys.stderr.write("After CNoverlap, %d events\n" % (len(allevents)))
	if args.fleeting: 
		allevents=filter_fleeting_events(allevents, fcutoff=0) 
		sys.stderr.write("After fleeting, %d events\n" % (len(allevents)))
	if args.zeroCost: 
		allevents=filter_zero_cost_events(allevents)
		sys.stderr.write("After zeroCost, %d events\n" % (len(allevents)))
	historyScores=np.loadtxt(open(args.histstats, 'r'), dtype=int)
	for e in allevents:
		totalp=histseg.compute_totalp(historyScores) 
		e.update(historyScores, totalp)
	pickle.dump(allevents, open(args.out, 'w'), pickle.HIGHEST_PROTOCOL)  


if __name__ == '__main__': 
	parser=argparse.ArgumentParser("This will filter events based on the options selected.")
	parser.add_argument('pevnt', help='the pickled event file (.pevnt)')
	parser.add_argument('histstats', help='historystats.txt file')
	parser.add_argument('out', help='The name of the output file.')
	parser.add_argument('--fleeting', action='store_true', help='Filter events based on being fleeting.')
	parser.add_argument('--zeroCost', action='store_true', help='Filter events based on having a cost of zero.')
	parser.add_argument('--CNoverlap', action='store_true', help='Filter events based on CN overlap')
	parser.add_argument('--bed', help='A bed file of the CN changes')
	parser.add_argument('--cbs', help='The file is a CBS_OUT file of the CN changes', action='store_true')
	parser.add_argument('--perc', help='A percentage overlap', default=0.9, type=float)
	args=parser.parse_args()
	main(args)

