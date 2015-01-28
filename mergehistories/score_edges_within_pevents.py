#!/usr/bin/env python 
#This is used to look at individual edges (segments or adjacencies) across histories as opposed to the cyclic events.   

import sys, os
import copy
import cPickle as pickle 
import cnavgpost.mergehistories.event_cycles_module as histseg
import argparse
import numpy as np
bindir="/inside/home/tballing/cnavg-study/cnavgmbin"

def unique_c_edges(sorted_edges): 
	unique_edges=[]
	prevedge=sorted_edges[0]
	for edge in sorted_edges[1:]: 
		if edge == prevedge: 
			consolidate_values(prevedge, edge) #This will modify prevedge
		else: 
			unique_edges.append(prevedge)
			prevedge=edge
	unique_edges.append(prevedge)
	return unique_edges

def consolidate_values(edge1, edge2): 
	for i in xrange(len(edge2.histories)):
		hid = edge2.histories[i] 
		if hid not in edge1.histories: 
			edge1.histories.append(hid)
			edge1.prevals.append(edge2.prevals[i])
			edge1.orders.append(edge2.orders[i])
			edge1.uppercosts.append(edge2.orders[i])
			edge1.lowercosts.append(edge2.orders[i])

def unique_loc_edges(sorted_edges): 
	unique_edges=[]
	prevedge=sorted_edges[0]
	prevsegstr=prevedge.segstr
#	(prevsegstr, sign)=histseg.remove_signs_from_segstr(prevedge.segstr)
	for edge in sorted_edges[1:]:
		#(mysegstr, sign)=histseg.remove_signs_from_segstr(edge.segstr) 
		mysegstr=edge.segstr
		if mysegstr== prevsegstr: 
			consolidate_values_sumCN(prevedge, edge) #This will modify prevedge
		else: 
			unique_edges.append(prevedge)
			prevedge=edge
			prevsegstr=mysegstr
	unique_edges.append(prevedge)
	return unique_edges
	
def consolidate_values_sumCN(edge1, edge2):
	edge1.cnval+=edge2.cnval 
	for i in xrange(len(edge2.histories)):
		hid = edge2.histories[i] 
		if hid not in edge1.histories: 
			edge1.histories.append(hid)
			edge1.prevals.append(edge2.prevals[i])
			edge1.orders.append(edge2.orders[i])
			edge1.uppercosts.append(edge2.orders[i])
			edge1.lowercosts.append(edge2.orders[i])


def score_edges_within_pevents(allevents, historyScores, totalp, prev_error=0.05, ignore_cn=True): 
	prevalence_error=prev_error
	sys.stderr.write("number of events: %d\n" % (len(allevents)))
	sys.stderr.write("ignore_cn: %s\n" % ignore_cn) 
	alledges=[]
	for event in allevents:
		event.unpack() 
		for seg in event.segs: 
			edge=copy.deepcopy(event)
			edge.segs=[seg]
			edge.make_segstr()
			alledges.append(edge)
	sortededges=sorted(alledges, key=lambda x: (x.segstr, x.cnval))
	if ignore_cn: 
		unique_edges=histseg.merge_events_by_type(sortededges)  #unique_loc_edges(sortededges)
	else: 
		unique_edges=histseg.unique_c_events_sorted_list(sortededges)
		splitoffs=histseg.get_split_offs(unique_edges)
		unique_edges+=splitoffs
	sys.stderr.write("totalp: %s\n" % (str(totalp)))
	sys.stderr.write("number of edges is: %d\n" % (len(unique_edges)))
	for edge in unique_edges: 
		edge.update(historyScores)
		edge.likelihood=histseg.compute_likelihood_histories(edge.histories, historyScores, totalp)
		edge.trim()
	return unique_edges 

def add_score_edges_options(parser): 
	parser.add_argument('pevnts', help='a .pevnts file.')
	parser.add_argument('historystats', help='The file with historystats')
	parser.add_argument('--outpickle', help='pickle the edges and write them to this file.')
	parser.add_argument('--edges', help='write edges to this file in text (not pickled).')
	parser.add_argument('--prevalence_error', help='the difference in prevalences to be considered the same.', type=float, default=0.05)
	parser.add_argument('--ignore_cn', help='merge together edges with different CN values.', default=False, action='store_true')
	parser.add_argument('--totalp', help='total probability of the histories', type=float)
	parser.add_argument('--binwidth', help='the multiplier between history ids of independent runs', default=histseg.Global_BINWIDTH, type=int)	

if __name__ == "__main__": 
	parser = argparse.ArgumentParser(description='Given an .pevnts file, it will split events into edges, combine equivalent edges (segments or adjacencies), and score them by likelihood.')
	add_score_edges_options(parser)
	args=parser.parse_args()
	histseg.Global_BINWIDTH=args.binwidth
	allevents=pickle.load(open(args.pevnts, 'rb'))
	historyScores=np.loadtxt(args.historystats, dtype=int)
	totalp=0
	if args.totalp: 
		totalp=args.totalp
	else: 
		totalp = histseg.compute_likelihood_histories(historyScores[:,0], historyScores) 
	alledges = score_edges_within_pevents(allevents, historyScores, totalp, args.prevalence_error, args.ignore_cn)
	if args.edges:
		outfile=open(args.edges, 'w') 
		for edge in alledges: 
			outfile.write(str(edge))
	if args.outpickle: 
		pickle.dump(alledges, open(args.outpickle, 'wb'), pickle.HIGHEST_PROTOCOL)


