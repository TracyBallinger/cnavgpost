#!/inside/home/common/bin/python2.7 

import sys, os
import re
from numpy import * 
import history_segments_module as histseg

def calculate_history_stats(values, histids):
	myhistids={}
	histid_crossindex=0
	for id in histids: 
		if id not in myhistids.keys(): 
			myhistids[id] = histid_crossindex
			histid_crossindex += 1
	numhistories=len(myhistids.keys())
	maxpts=max(myhistids, key=lambda x: myhistids[x])
	mydistances=ones((len(values), numhistories), dtype=float) * 100 
	# find the minimum difference in CN values between two points for every history combination
	for i in xrange(0, len(values)): 
		valA=values[i]
		histA=histids[i]
		for j in xrange(0, len(values)): 
			valB=values[j]
			histB=histids[j]
			histBi=myhistids[histB]
			dist=abs(valA-valB)
			if dist < mydistances[i][histBi]: 
				mydistances[i][histBi]=dist
	mystats=array([mydistances.mean(axis=1), mydistances.std(axis=1)])
	return(mystats)

if __name__ == "__main__": 
	import argparse
	parser = argparse.ArgumentParser(description='Given a file of universal segments with histories (output from get_segment_values_dist.py), this will calculate the ave and sd for events across the histories')
	parser.add_argument('segment_histories', help='a file with a comma separated list of CN values, timings, prevalences, and history IDs') 
	args=parser.parse_args()
	infile=open(args.segment_histories, 'r')
	sys.stdout.write("length\tn\tcn_value\tcn_mean\tcn_sd\tprevalence\tp_mean\tp_std\n")
	for line in infile: 
		myseghist=histseg.Segment_history(line)
		myhistids={}
		for id in myseghist.histories: 
			if id in myhistids: 
				myhistids[id]+= 1
			else: 
				myhistids[id]=1 
		numhistories=len(myhistids.keys())
		cnstats=calculate_history_stats(myseghist.cnvals, myseghist.histories)
		prevalstats=calculate_history_stats(myseghist.prevals, myseghist.histories)
		for i in xrange(0,len(myseghist.cnvals)): 
			sys.stdout.write("%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n" % (myseghist.end-myseghist.start, numhistories, myseghist.cnvals[i], cnstats[0,i], cnstats[1,i], myseghist.prevals[i], prevalstats[0,i], prevalstats[1,i]))


