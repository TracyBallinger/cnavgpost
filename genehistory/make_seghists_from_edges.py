#!/usr/bin/env python 
import sys, os
import pickle
from cnavgpost.mergehistories.history_segments_module import * 
import cnavgpost.mergehistories.event_cycles_module as ecycles

def create_seghists_from_edges(sortededges): 
	myseghists=[]
	working_seghists=[SegmentHistory(sortededges[0])]
	for edge in sortededges[1:]:
		tmpseghists=[]
		for seghist in working_seghists: 
			if seghist.overlaps(edge): 
				new_seghists=merge_in_edge(edge, seghist)
			 	tmpseghists+=new_seghists
			elif seghist.comes_after(edge): 
			 	tmpseghists.append(seghist)
			elif seghist.comes_before(edge): 
				myseghists.append(seghist)
		working_seghists=tmpseghists
	myseghists.append(working_seghists)
	return(myseghists)
	

if __name__ == "__main__": 
	import argparse
	parser = argparse.ArgumentParser(description='Will take the edges from events and split the genome at breakpoints, with histories for each segment.')
	parser.add_argument('pegdes', help='a pickled files of event edges.')
	args=parser.parse_args()
	edges=pickle.load(open(args.pedges, 'rb'))
	sortededges=sort(edges, key=lambda x: (x.segstr, x.cnval))
	seghists=create_seghists_from_edges(sortededges)
	pickle.dump(seghists, stdout, pickle.HIGHEST_PROTOCOL)

	
