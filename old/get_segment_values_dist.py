#!/inside/home/common/bin/python2.7

import sys, os
from numpy import array 
import copy 
import history_segments_module as histseg
import argparse

parser = argparse.ArgumentParser(description='Given a file of segments from different histories, it will output the regions that overlap and the history information in that segment.')
parser.add_argument('history_segments', type=argparse.FileType('r'), default=sys.stdin, help='a concatenated file of segments from multiple histories, ordered by genomic location.  (sort -k1,1 -k2,2n -k3,3n -k8,8n history_?.txt)')

args=parser.parse_args()
infile=args.history_segments

line=infile.readline()
current_seghist=histseg.Segment_history(line)
working_seghists=[current_seghist]
#Should put in a check somewhere that the input is properly ordered by chr location. 

for line in infile: 
	myseg=histseg.Segment_history(line)
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
			sys.stdout.write(str(workseg) + "\n")
			
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
	sys.stdout.write(str(seg) + "\n")

