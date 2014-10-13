#!/usr/bin/env python 

# For use in determining the number of TP and FP, etc, of events across a CNAVG sampling of histories given that one of the histories is the truth.  This history is what all other will be compared to.  
import os, sys
import cnavgpost.mergehistories.event_cycles_module as histseg
import cPickle as pickle
import argparse
import numpy as np
import re

class EdgeSimulationData:
	def __init__(self, edge): 
		self.edge=edge
		self.isTrue=0  # this will be 0 if edge is FP, 1 if TP, 2 if TN, -1 if FN 
		self.refindex=-1
		self.refpreval=-1
		self.reforder=-1
		self.avecost=-1
		self.type=edge.determineEventType()
		(self.segstr, self.sign)=histseg.remove_signs_from_segstr(edge.segstr)
		self.cnval=edge.cnval*self.sign		

def analyze_simulation(edges, refhistoryid, historyScores, datout_fh, stats_fh, breaks_fh):
	#make the cost of the refhistoryid 0 so that is doesn't get included in the likelihood calculation 
	myhistScores=np.copy(historyScores)
	myhistScores[np.where(historyScores[:,0] == refhistoryid),:]=0	 
	totalp=histseg.compute_likelihood_histories(myhistScores[:,0], myhistScores)
	TP=[0,0,0,0]
	FP=[0,0,0,0]
	TN=[0,0,0,0]
	FN=[0,0,0,0]
	FNedges=[]
	types=histseg.Global_EVENTTYPES
	myEdgeSimData=[] # a list of tuples (edge, isTrue, refpreval, reforder)
	for edge in edges:
		if not edge.histories: edge.histories=histseg.listout_ranges(edge.histRanges)
		myedgesim=EdgeSimulationData(edge)
		type=myedgesim.type
		if refhistoryid in edge.histories: 
			refindex=edge.histories.index(refhistoryid)
			myedgesim.refindex=refindex
			if len(edge.histories)>1:
				TP[0]+=1
				TP[type]+=1
				myedgesim.isTrue=1
				edge.histories.pop(refindex)
				myedgesim.refpreval=edge.prevals.pop(refindex)
				myedgesim.reforder=edge.orders.pop(refindex)
				edge.likelihood = histseg.compute_likelihood_histories(edge.histories, myhistScores, totalp)	
				edge.compute_timing_wmeansd(myhistScores)
				edge.histories.insert(refindex, refhistoryid)	
				edge.prevals.insert(refindex, myedgesim.refpreval)	
				edge.orders.insert(refindex, myedgesim.reforder)	
				upperc=edge.uppercosts.pop(refindex)
				lowerc=edge.lowercosts.pop(refindex)
				myedgesim.avecost=np.mean(np.array(edge.uppercosts+edge.lowercosts))
				edge.uppercosts.insert(refindex, upperc)	
				edge.lowercosts.insert(refindex, lowerc)	
			else: 
				FN[0]+=1
				FN[type]+=1
				FNedges.append(myedgesim)
				myedgesim.isTrue=-1
				edge.likelihood=1
				myedgesim.avecost=np.mean(np.array(edge.uppercosts+edge.lowercosts))
				myedgesim.refpreval=edge.prevals[refindex]
				myedgesim.reforder=edge.orders[refindex]
		else: 
			FP[0]+=1
			FP[type]+=1
			edge.likelihood = histseg.compute_likelihood_histories(edge.histories, myhistScores, totalp)	
			if edge.likelihood >1: 
				sys.stderr.write("bad lscore: %s\t%s\t%d\n" % (str(edge.likelihood), str(totalp), len(edge.costs)))
			myedgesim.isTrue=0
			myedgesim.avecost=np.mean(np.array(edge.uppercosts+edge.lowercosts))
		myEdgeSimData.append(myedgesim)
	if len(FNedges) >0: 
		TN=checkForCancellingEdges(FNedges) #this will also modify the isTrue value of FNedges
		for i in xrange(len(TN)): 
			FN[i]=FN[i]-TN[i]	
	
	if datout_fh: 
		header="event_id\tevent_type\tavecost\tLscore\tCNval\ttrue\tlength\tprevals\torders\tnumhists\n"
		datout_fh.write(header)
		for edgesim in myEdgeSimData:
			edge=edgesim.edge 
			prevals=",".join(map(str, [edgesim.refpreval, edge.prevalmean, edge.prevalsd]))
			orders=",".join(map(str, [edgesim.reforder, edge.ordermean, edge.ordersd]))
			type=edge.determineEventType()
			length=edge.get_Event_length()
			mystr="\t".join(map(str, [edge.id, types[type], edgesim.avecost, edge.likelihood, edge.cnval, edgesim.isTrue, length, prevals, orders, len(edge.histories)])) + "\n"	
			datout_fh.write(mystr)
	
	if stats_fh: 
		stats_fh.write("type\ttotal\tAmp\tDel\tAdj\n")
		stats_fh.write("TP\t%s\nFP\t%s\nFN\t%s\nTN\t%s\n" % ("\t".join(map(str, TP)), "\t".join(map(str, FP)), "\t".join(map(str, FN)), "\t".join(map(str, TN)) ))
		f1score = float(2*TP[0])/float(2*TP[0]+FN[0]+FP[0])
		stats_fh.write("F1Score:\t%s\n" % (str(f1score)))
	
	if breaks_fh: 
		breakpoints=histseg.get_breakpoints(edges, refhistoryid)
		for loc in breakpoints.keys(): 
			(n, t) = breakpoints[loc]
			breaks_fh.write("%s\t%d\t%d\n" % (loc, n, t))
		breaks_fh.write("Breakpoints: %d\n" % len(breakpoints))


def checkForCancellingEdges(edgesims):
	sortededgesims=sorted(edgesims, key=lambda x: (x.segstr, x.edge.cnval))  
	preseg=sortededgesims[0]
	numpairs=0
	TN=[0,0,0,0]
	types=histseg.Global_EVENTTYPES
	for edgesim in sortededgesims[1:]: 
		if edgesim.segstr == preseg.segstr and (edgesim.cnval == -1 * preseg.cnval):
			numpairs+=1
			edgesim.isTrue=2
			preseg.isTrue=2
			TN[0]+=2
			TN[preseg.type]+=1
			TN[edgesim.type]+=1		
		preseg=edgesim
	return TN

def add_options(parser):
	parser.add_argument('pedgs', help='a .pedgs file')
	parser.add_argument('trueID', help='the ID of the true history', type=int)
	parser.add_argument('historystats', help='The file of history stats')
	parser.add_argument('--datout', help='the file to write data to', type=argparse.FileType('w')) #, default=sys.stdout) 
	parser.add_argument('--stats', help='the file to write stats to.', type=argparse.FileType('w')) #, default=sys.stderr) 
	parser.add_argument('--breaks', help='the file to write breaklocations to.', type=argparse.FileType('w')) #, default=sys.stderr) 
	parser.add_argument('--binwidth', dest='binwidth', help='the multiplier between history ids of independent runs', default=histseg.Global_BINWIDTH, type=int)	

if __name__== '__main__': 
	parser = argparse.ArgumentParser(description='does analysis for simulation tests')
	add_options(parser)
	args=parser.parse_args()
	histseg.Global_K=0  #this used to be 1 initially.  
	histseg.Global_BINWIDTH=args.binwidth
	allevents=pickle.load(open(args.pedgs, 'rb'))
	historyScores=np.loadtxt(args.historystats, dtype=int)
	analyze_simulation(allevents, args.trueID, historyScores, args.datout, args.stats, args.breaks)

