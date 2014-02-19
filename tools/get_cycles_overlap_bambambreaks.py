#!/inside/home/common/bin/python2.7

import sys, os
import argparse
import event_cycles_module as histseg
import pysam
import re

def get_cycles_overlapping_break(breakline, evntfn): 
	evntfile=open(evntfn, 'r')
	myevents=[]
	for evntline in evntfile: 
		if break_overlap_event(breakline, evntline): 
			myevents.append(evntline)
	return myevents

class bambamBreak: 
	def __init__(self, line): 
		data=line.strip().split('\t')
		loc1=data[0]
		m=re.search('(\w+):(\d+)-(\d+)', loc1)
		(chr, start, end) = m.group(1,2,3)
		(self.chr1, self.start1, self.end1)=(chr, int(start), int(end))
		loc2=data[1]
		m=re.search('(\w+):(\d+)-(\d+)', loc2)
		(chr, start, end) = m.group(1,2,3)
		(self.chr2, self.start2, self.end2)=(chr, int(start), int(end))
		self.str1=data[2]
		self.str2=data[3]	

def break_overlap_event(breakline, evntline): 
	mybreak=bambamBreak(breakline)
	myevent=histseg.Event(evntline)
	myevent.make_segs_from_str()
	overlap=False
	for seg in myevent.segs: 
		if seg.adj:
			if seg.chr.startswith("chr"): seg.chr=seg.chr[3:] 
			if seg.chr2.startswith("chr"): seg.chr2=seg.chr2[3:] 
			if (((seg.chr==mybreak.chr1 and seg.start>=mybreak.start1 and seg.start <= mybreak.end1) 
				and (seg.chr2==mybreak.chr2 and seg.end >= mybreak.start2 and seg.end <= mybreak.end2))
				or ((seg.chr==mybreak.chr2 and seg.start >= mybreak.start2 and seg.start <= mybreak.end2)
				and (seg.chr2==mybreak.chr1 and seg.end >= mybreak.start1 and seg.end <= mybreak.end1))):
				overlap=True
	return overlap


if __name__ == '__main__': 
	import argparse
	parser = argparse.ArgumentParser(description='Given a .breaks file from bambam and a .evnts file containing cycles from cnavg, will find the cycles that contain adjencies overlapping the breaks in the breaks file.')
	parser.add_argument('breaks', help='bambam breaks file')
	parser.add_argument('evnts', help='CNAVG evnts file')
	args=parser.parse_args()
	breakf=open(args.breaks, 'r')
	for breakline in breakf: 
		myevents=get_cycles_overlapping_break(breakline, args.evnts)
		sys.stderr.write("# events: %d\n" % (len(myevents)))
		for evnt in myevents:
			sys.stdout.write(evnt)

