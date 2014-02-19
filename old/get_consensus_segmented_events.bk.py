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
		return "\t".join([mystr, str(numpy.mean(self.cnvals)),  str(numpy.mean(self.prevals)), ",".join(map(str, self.cnvals)), ",".join(map(str, self.prevals))])

def split_segment_by_CN(segment_hist, merge_distance): 
	myevents=[]
	values=segment_hist.cnvals
	eventA=Consensus_event(segment_hist)
	eventA.cnvals.append(segment_hist.cnvals[0])
	eventA.prevals.append(segment_hist.prevals[0])
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
	# if any of the clusters are identical, then just keep one of them 
	tmpsetlist=[]
	prev_numsets=len(tmpsetlist)
	while len(setlist) != prev_numsets:
		prev_numsets=len(tmpsetlist)
		for set1 in setlist: 
			tmpsetlist.append(set1)
			for set2 in setlist[1:]:
				if set2 != set1:
					tmpsetlist.append(set2) 
		setlist=tmpsetlist		

	sys.stderr.write("eventA is %s\n" % (str(eventA)))
	unmerged=range(1, len(values))
	while len(unmerged) > 0:
		sys.stderr.write("unmerged is: " + str(unmerged) + "\n")
		valA=eventA.cnvals[0]
		tmpunmerged=[] 
		for i in unmerged :  
			valB=values[i]
			sys.stderr.write("valA: %f, valB: %f, \n" % (valA, valB))
			if abs(valA - valB) < merge_distance and (valA * valB >= 0) :  #check they are the same sign
				eventA.cnvals.append(valB)
				eventA.prevals.append(segment_hist.prevals[i])
				sys.stderr.write("updated eventA: %s\n" % (str(eventA)))
			else: 
				tmpunmerged.append(i)
		sys.stderr.write("tmpunmerged is " + str(tmpunmerged) + "\n")
		if len(tmpunmerged) == len(unmerged): # if no more values were merged into this event, make a new one
			myevents.append(copy.deepcopy(eventA))
			eventA=Consensus_event(segment_hist)
			eventA.cnvals.append(segment_hist.cnvals[tmpunmerged[0]])
			eventA.prevals.append(segment_hist.prevals[tmpunmerged[0]])
			sys.stderr.write("new eventA: %s\n" % (str(eventA)))
			unmerged=tmpunmerged[1:]
		else:
			unmerged=tmpunmerged
	# add the last event on
	myevents.append(eventA)
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
			sys.stderr.write("workseg: %s\n" % (str(workseg))) 
			consensus_events=split_segment_by_CN(workseg, merge_distance)
			for event in consensus_events: 
				sys.stdout.write(str(event) + "\n")
				sys.stderr.write(str(event) + "\n")
			
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

