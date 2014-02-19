#!/inside/home/common/bin/python2.7

import os, sys
import event_cycles_module as histseg
import pickle
import argparse

def analyze_simulation(edges, args): 
	totalp=subtract_from_likelihood(edges, args.trueID)
	refhistoryid=args.trueID
	TP=[0,0,0,0]
	FP=[0,0,0,0]
	TN=[0,0,0,0]
	FN=[0,0,0,0]
	FNedges=[]
	header="edge_type\tLscore\tCNval\ttrue\tlength\tprevals\torders\tnumhists\n"
	args.datout.write(header)
	types=histseg.Global_EVENTTYPES
	for edge in edges: 
		type=edge.determineEventType()
		length=edge.get_Event_length()
		isTrue=0
		refindex=0
		refpreval=-1
		reforder=-1
		if refhistoryid in edge.histories: 
			refindex=edge.histories.index(refhistoryid)
			if len(edge.histories)>1:
				edge.histories.pop(refindex)
				refpreval=edge.prevals.pop(refindex)
				reforder=edge.orders.pop(refindex)
				refcosts=edge.costs.pop(refindex)
				edge.likelihood = histseg.compute_likelihood(edge.costs, totalp)	
				edge.compute_timing_wmeansd()
				TP[0]+=1
				TP[type]+=1
				isTrue=1
			else: 
				FN[0]+=1
				FN[type]+=1
				isTrue=2
				FNedges.append(edge)
				edge.likelihood=1
				refpreval=edge.prevals[0]
				reforder=edge.orders[0]
				refcosts=edge.costs[0]
		else: 
			FP[0]+=1
			FP[type]+=1
			edge.likelihood = histseg.compute_likelihood(edge.costs, totalp)
		prevals=",".join(map(str, [refpreval, edge.prevalmean, edge.prevalsd]))
		orders=",".join(map(str, [reforder, edge.ordermean, edge.ordersd]))
		mystr="%s\t%s\t%f\t%d\t%d\t%s\t%s\t%d\n" % (types[type], str(edge.likelihood), edge.cnval, isTrue, length, prevals, orders, len(edge.histories))	
		args.datout.write(mystr)
	if len(FNedges) >0:
		TN[0]=checkForCancellingEdges(FNedges)
	args.stats.write("type\ttotal\tAmp\tDel\tAdj\n")
	args.stats.write("TP\t%s\nFP\t%s\nFN\t%s\nTN\t%s\n" % ("\t".join(map(str, TP)), "\t".join(map(str, FP)), "\t".join(map(str, FN)), "\t".join(map(str, TN)) ))
	f1score = float(2*TP[0])/float(2*TP[0]+FN[0]+FP[0])
	args.stats.write("F1Score:\t%s\n" % (str(f1score)))

def checkForCancellingEdges(edges):
	mysegs=[]
	for edge in edges: 
		edge.make_segs_from_str()
		mysegs += edge.segs
	sortedsegs=sorted(mysegs, key=lambda x: (x.chr, x.start, x.chr2, x.end, x.st1, x.st2))
	preseg=sortedsegs[0]
	numpairs=0
	for seg in sortedsegs[1:]: 
		if seg.same_coords(preseg) and (seg.cnval == -1 * preseg.cnval):
			numpairs+=1
		preseg=seg
	return numpairs

def subtract_from_likelihood(events, refid):
	allhistoryids=[]
	allcosts=[]
	for event in events:
		for i in xrange(len(event.histories)):
			hid = event.histories[i]
			if hid != refid and hid not in allhistoryids:
				allhistoryids.append(hid)
				allcosts.append(event.costs[i])
	totalp = histseg.compute_likelihood(allcosts, 1)
	return totalp

def add_options(parser):
	parser.add_argument('pedgs', help='a .pedgs file')
	parser.add_argument('trueID', help='the ID of the true history', type=int)
	parser.add_argument('--datout', help='the file to write data to', type=argparse.FileType('w')) #, default=sys.stdout) 
	parser.add_argument('--stats', help='the file to write stats to.', type=argparse.FileType('w')) #, default=sys.stderr) 
	

if __name__== '__main__': 
	parser = argparse.ArgumentParser(description='does analysis for simulation tests')
	add_options(parser)
	args=parser.parse_args()
	allevents=pickle.load(open(args.pedgs, 'rb'))
	analyze_simulation(allevents, args)

