#!/inside/home/common/bin/python2.7

import sys, os
import argparse
import numpy
import copy 
from sets import Set 

parser = argparse.ArgumentParser(description='Given a file of segments from different histories, it will output the regions that overlap and the history information in that segment.')
parser.add_argument('history_segments', help='a concatenated file of segments from multiple histories, ordered by genomic location.  (sort -k1,1 -k2,2n -k3,3n -k8,8n history_?.txt)')
parser.add_argument('-m', '--merge_distance', type=float, help='The minimum CN distance two events in a segment have to be in order to be called the same event')
args=parser.parse_args()
merge_distance=args.merge_distance
infile=open(args.history_segments, 'r')

class Segment_history: 
	def __init__(self, myline):
		data=myline.strip().split('\t')
		self.chr=data[0]
		self.start=int(data[1])
		self.end=int(data[2])
		self.cnvals=[float(data[3])]
		self.timings=[int(data[4])]
		self.prevals=[float(data[5])]
		self.histories=[int(data[7])]
		self.id=data[0] + ":" + data[1] + "-" + data[2]
	def addline(self, line):
		data=line.strip().split('\t')
		self.cnvals.append(float(data[3]))
		self.timings.append(int(data[4]))
		self.prevals.append(float(data[5]))
		self.histories.append(int(data[7]))
	def __str__(self):
		mystr="\t".join([self.chr, str(self.start), str(self.end)])
		return "\t".join([mystr, ",".join(map(str, self.cnvals)), ",".join(map(str, self.timings)), ",".join(map(str, self.prevals)), ",".join(map(str, self.histories))])

class Consensus_event: 
	def __init__(self, myseg_hist):
		self.chr=myseg_hist.chr
		self.start=myseg_hist.start
		self.end=myseg_hist.end
		self.cnvals=[]
		self.prevals=[]
	def __str__(self):
		mystr="\t".join([self.chr, str(self.start), str(self.end)])
		return "\t".join([mystr, str(numpy.mean(self.cnvals)),  str(numpy.mean(self.prevals)), str(len(self.cnvals)), ",".join(map(str, self.cnvals)), ",".join(map(str, self.prevals))])

def split_segment_by_CN(segment_hist, merge_distance): 
	values=segment_hist.cnvals
	# get a list of clusters of events
	setlist=[]
	for iA in xrange(0, len(values)): 
		valA=values[iA]
		myset=Set([iA])
		for iB in xrange(0, len(values)):
			valB=values[iB]
			if (iB != iA) and abs(valA - valB) < merge_distance and (valA * valB >= 0) :  #check they are the same sign
				myset.add(iB)
		setlist.append(myset)
	# merge identical clusters. 
	finalsetlist=[]
	while 0 < len(setlist):
		set1=setlist[0]
		finalsetlist.append(set1)
		tmpsetlist=[]
		for set2 in setlist[1:]:
			if (set1 < set2): # if set1 is just a subset of set2, then replace it with the bigger set
				finalsetlist[-1]=set2
				set1=set2
			elif not (set2 <= set1): #if they aren't equivalent or subsets, then add set2 as a unique set
				tmpsetlist.append(set2) 
		setlist=tmpsetlist
	# create events from the sets
	myevents=[]
	for myset in finalsetlist:
		myevent=Consensus_event(segment_hist)
		myevent.cnvals=[segment_hist.cnvals[i] for i in myset]
		myevent.prevals=[segment_hist.prevals[i] for i in myset]
		myevents.append(myevent)
	return myevents 

line=infile.readline()
current_seghist=Segment_history(line)
working_seghists=[current_seghist]
#Should put in a check somewhere that the input is properly ordered by chr location. 

for line in infile: 
	myseg=Segment_history(line)
	working_tmphists=[]
	for workseg in working_seghists: 
		# Check if this segment overlaps the current one
		if (myseg.chr == workseg.chr) & (myseg.start <= workseg.end) & (myseg.end >= workseg.start): 
			# if this segment completely overlaps the current one, then add it. 	
			if (myseg.start <= workseg.start) & (myseg.end >= workseg.end): 
				workseg.addline(line)
				working_tmphists.append(workseg)
			# if this segment partially overlaps the current one, split the current one into or 2 or 3 non-overlapping ones
			else: 
				if (myseg.start > workseg.start):
					tmpseg=copy.deepcopy(workseg)
					tmpseg.end=myseg.start-1
					working_tmphists.append(tmpseg)
				else: # myseg.start <= workseg.start
					myseg.start=workseg.start
				tmpseg=copy.deepcopy(workseg)
				tmpseg.start=myseg.start
				tmpseg.addline(line)
				tmpseg.end = min(myseg.end, workseg.end)
				working_tmphists.append(tmpseg)
				if (myseg.end < workseg.end):
					tmpseg=copy.deepcopy(workseg)
					tmpseg.start=myseg.end+1
					working_tmphists.append(tmpseg)
		elif (workseg.chr == myseg.chr) & (myseg.end < workseg.start): 
			working_tmphists.append(workseg)
		elif (myseg.chr > workseg.chr) | (myseg.start > workseg.end):
			# print out the current working segment - we're done with it.
			consensus_events=split_segment_by_CN(workseg, merge_distance)
			for event in consensus_events: 
				sys.stdout.write(str(event) + "\n")
			
	# after checking for overlap with all of the current working segments
	workseg=working_seghists[-1]
	if (myseg.chr == workseg.chr) & (myseg.start <= workseg.end) & (myseg.end > workseg.end): 
		myseg.start=workseg.end+1
		working_tmphists.append(myseg)
	elif (myseg.chr > workseg.chr) | (myseg.start > workseg.end):
		working_tmphists.append(myseg)
	working_seghists=working_tmphists

# print out all of the rest of the segment at the end
for seg in working_seghists:
	consensus_events=split_segment_by_CN(seg, merge_distance)
	for event in consensus_events: 
		sys.stdout.write(str(event) + "\n")

