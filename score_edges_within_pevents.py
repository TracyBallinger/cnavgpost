#!/inside/home/common/bin/python2.7
import sys, os
import copy, pickle 
import event_cycles_module as histseg
import argparse
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
			edge1.costs.append(edge2.costs[i])
			edge1.prevals.append(edge2.prevals[i])
			edge1.orders.append(edge2.orders[i])

def score_edges_within_pevents(args): 
	global prevalence_error
	prevalence_error=args.prevalence_error
	allevents=pickle.load(open(args.pevnts, 'rb'))
	sys.stderr.write("number of events: %d\n" % (len(allevents)))
	alledges=[]
	allcosts=[]
	allhistoryids=[]
	for event in allevents: 
		if not args.totalp: 
			for i in xrange(len(event.histories)): 
				hid = event.histories[i]
				if hid not in allhistoryids: 
					allhistoryids.append(hid)
					allcosts.append(event.costs[i])
		event.make_segs_from_str()
		edges=event.segs
		for seg in event.segs: 
			edge=histseg.Event([seg]) # copy.deepcopy(event)
			edge.histories=event.histories
			edge.prevals=event.prevals
			edge.orders=event.orders
			edge.costs=event.costs
			edge.make_segstr()
			alledges.append(edge)
	#	sys.stderr.write("number of edges is: %d\n" % (len(alledges)))
	sortededges=sorted(alledges, key=lambda x: (x.segstr, x.cnval))
	unique_edges=unique_c_edges(sortededges)
	totalp= histseg.compute_likelihood(allcosts, 1)
	sys.stderr.write("totalp: %s\n" % (str(totalp)))
	for edge in unique_edges: 
		edge.numhists=len(edge.histories)
		edge.compute_timing_wmeansd()
		edge.get_hranges()
		edge.likelihood=histseg.compute_likelihood(edge.costs, totalp)
		edge.trim()
	if args.edges:
		outfile=open(args.edges, 'w') 
		for edge in unique_edges: 
			outfile.write(str(edge))
	if args.outpickle: 
		outfile=open(args.outpickle, 'wb')
		pickle.dump(unique_edges, outfile, pickle.HIGHEST_PROTOCOL)

def add_score_edges_options(parser): 
	parser.add_argument('pevnts', help='a .pevnts file.')
	parser.add_argument('--outpickle', help='pickle the edges and write them to this file.')
	parser.add_argument('--edges', help='write edges to this file.')
	parser.add_argument('--prevalence_error', help='the difference in prevalences to be considered the same.', type=float, default=0.05)
	parser.add_argument('--totalp', help='total probability of the histories', type=float)

if __name__ == "__main__": 
	parser = argparse.ArgumentParser(description='Given an .pevnts file, it will split events into edges, combine equivalent edges (segments or adjacencies), and score them by likelihood.')
	add_score_edges_options(parser)
	args=parser.parse_args()
	score_edges_within_pevents(args) 
