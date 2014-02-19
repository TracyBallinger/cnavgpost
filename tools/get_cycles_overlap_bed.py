#!/inside/home/common/bin/python2.7

import sys, os
import event_cycles_module as histseg
import subprocess, StringIO

def get_cycles_overlapping_bed(chr, start, end, evntfn): 
	evntfile=open(evntfn, 'r')
	myevents=[]
	for evntline in evntfile: 
		if bed_overlap_event(chr, start, end, evntline): 
			myevents.append(evntline)
	return myevents

def bed_overlap_event(chr, start, end, evntline):
	myevent=histseg.Event(evntline)
	myevent.make_segs_from_str()
	overlap=False
	for seg in myevent.segs: 
		if seg.adj:
			if ((seg.chr == chr and seg.start >= start and seg.start<=end) or (seg.chr2==chr and seg.end >= start and seg.start <= end)):
				overlap=True
		else: # seg is a segment, not an adjacency
			if (seg.chr ==chr and seg.start <= end and seg.end >= start): 
				overlap=True
	return overlap


if __name__ == "__main__":	
	import argparse
	parser = argparse.ArgumentParser(description='Will output the subset of events from .braney file that overlap the specified genomic region.  Will filter out duplicate identical events and indicate the number of histories that contain each version of an event.')
	parser.add_argument('bed', help='a bedfile of regions to select events for')
	parser.add_argument('evnts', help='a CN-AVG evnts file.')
	args=parser.parse_args()
	for bedline in open(args.bed, 'r'):
		beddata=bedline.strip().split('\t')
		(chr, start, end)=beddata[0:3]
		myevents=get_cycles_overlapping_bed(chr, int(start), int(end), args.evnts)
		for evnt in myevents: 
			sys.stdout.write(evnt)		

