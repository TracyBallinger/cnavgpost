#!/usr/bin/env python 

import argparse
import cPickle as pickle
import sys, os
import cnavgpost.mergehistories.event_cycles_module as histseg
import numpy as np
import re

def cluster_events(events): 
	mydat=np.zeros((len(events), len(events)), dtype=int)
	for i in xrange(len(events)):
		e1=events[i] 
		for j in xrange((i+1),len(events)): 
			e2=events[j]
			mydat[i,j]=sum(np.in1d(np.array(e1.histories), np.array(e2.histories)))
	return mydat

def main(args): 
	events=pickle.load(open(args.pevnts, 'rb'))
	for e in events: 
		e.histories=histseg.listout_ranges(e.histRanges)
	mymat=cluster_events(events)	
	np.savetxt(args.out, mymat, fmt='%d', delimiter="\t")

if __name__=='__main__': 
	parser=argparse.ArgumentParser(description='given a .pevnts file, will create a square matrix of the number of histories shared between every two events')
	parser.add_argument('pevnts', help='a .pvents file')
	parser.add_argument('out', help='name of the output file')
	args=parser.parse_args()
	main(args)	
