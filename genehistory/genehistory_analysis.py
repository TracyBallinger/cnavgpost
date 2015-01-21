#!/usr/bin/env python 
import os, sys
import argparse
import cPickle as pickle 
import cnavgpost.mergehistories.event_cycles_module as histseg 
from cnavgpost.mergehistories.history_node_links_module import *
import cnavgpost.genehistory.annotate_events as annmod
import numpy as np

def wholehistory_analysis(evnts, outdir, historyScores, cutoff): 
	datfn=os.path.join(outdir, "events.dat")
	datfile=open(datfn, 'w')
	header="event_id\tevent_type\tavecost\tLscore\tCNval\ttrue\tlength\tprevalmean\tprevalsd\tordermean\tordersd\tnumhists\n"
	datfile.write(header)
	goodevents=[]	
	for e in evnts: 
		if e.likelihood > cutoff: 
			e.unpack()
			goodevents.append(e)
			datfile.write(write_event_data(e))
	datfile.close() 
	linkfn=os.path.join(outdir, "links.dat")
	linkfile=open(linkfn, 'w')
	links=link_events_by_order_within_histories(goodevents)
	for l in links: 
		l.likelihood=histseg.compute_likelihood_histories(l.histories, historyScores)
		linkfile.write(str(l))
	linkfile.close()	

def write_event_data(event):
	#header="event_id\tevent_type\tavecost\tLscore\tCNval\ttrue\tlength\tprevals\torders\tnumhists\n"
	avecost=np.mean(np.array(event.uppercosts+event.lowercosts))
	types=histseg.Global_EVENTTYPES
	return "\t".join(map(str, [event.id, types[event.determineEventType()], avecost, event.likelihood, event.cnval, 0, event.get_Event_length(), event.prevalmean, event.prevalsd, event.ordermean, event.ordersd, event.numhists])) + "\n"

def single_gene_history(chr, start, end, evnts, outdir, historyScores, cutoff):
	myevents=[]
	datfn=os.path.join(outdir, "%s_%d-%d_event.dat" % (chr, start, end))
	datfile=open(datfn, 'w')
	header="chr\tstart\tend\tevent_id\tevent_type\tavecost\tLscore\tCNval\ttrue\tlength\tprevalmean\tprevalsd\tordermean\tordersd\tnumhists\n"
	datfile.write(header)
	for event in evnts: 
		if event.likelihood > cutoff:
			if event.check_overlap(chr, start, end): 
				datfile.write(print_overlapping_segs(event, chr, start, end))
	datfile.close()

def print_overlapping_segs(event, chr, start, end):
	mystr=""
	for seg in event.segs:
		if (	seg.seg and 
				(seg.chr == chr and 
				seg.start < end and seg.end>start)):
			mystr=mystr+"%s\t%d\t%d\t%s" % (chr, seg.start, seg.end, write_event_data(event))
			#mystr=mystr+"%s\t%s" % (event.segstr, write_event_data(event))
		else:
			if (seg.chr == chr and seg.start <= end and seg.start >= start):  
				mystr=mystr+"%s\t%d\t%d\t%s" % (chr, seg.start, seg.start, write_event_data(event))
			elif (seg.chr2 ==chr and seg.end <= end and seg.end >= start):
				mystr=mystr+"%s\t%d\t%d\t%s" % (chr, seg.end, seg.end, write_event_data(event))
	return mystr 
				

if __name__=='__main__': 
	parser = argparse.ArgumentParser(description='filters events based on p-value, annotates them, links them, and plots them')
	parser.add_argument('pevnts', help='a .pevnts file')
	parser.add_argument('outdir', help='the directory where you want the output to go')
	parser.add_argument('historystats', help='historystats.txt file')
	parser.add_argument('--genes', help='This is a text file that has a list of gene names or a bed file with gene coordinates.')
	parser.add_argument('--cutoff', help='the pvalue cutoff to use', default=0.5, type=float)
	parser.add_argument('--plotpdf', help='generate a scatter plot of the events')
	args=parser.parse_args()
	if not os.path.exists(args.outdir): 
		os.mkdir(args.outdir)
	events=pickle.load(open(args.pevnts, 'rb'))
	for e in events: 
		e.unpack()
	historyScores=np.loadtxt(args.historystats, dtype='int')
	if args.genes: 
		bedlines=open(args.genes, 'r').readlines()
		for bedline in bedlines:
			beddat=bedline.strip().split('\t') 
			chr=beddat[0]
			start=int(beddat[1])
			end=int(beddat[2])
			single_gene_history(chr, start, end, events, args.outdir, historyScores, args.cutoff)	
	else: 
		wholehistory_analysis(events, args.outdir, historyScores, args.cutoff)


	
