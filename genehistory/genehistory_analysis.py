#!/usr/bin/env python 
import os, sys
import argparse
import cPickle as pickle 
import cnavgpost.mergehistories.event_cycles_module as histseg 
from cnavgpost.mergehistories.history_node_links_module import *
import cnavgpost.genehistory.annotate_events as annmod
import cnavgpost.genehistory.bedFileModule as bedmod
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

def single_gene_history(chr, start, end, evnts, outdir):
	myevents=[]
	datfn=os.path.join(outdir, "%s_%d-%d_event.dat" % (chr, start, end))
	datfile=open(datfn, 'w')
	header="chr\tstart\tend\tevent_id\tevent_type\tavecost\tLscore\tCNval\ttrue\tlength\tprevalmean\tprevalsd\tordermean\tordersd\tnumhists\n"
	datfile.write(header)
	for event in evnts: 
		if event.check_overlap(chr, start, end): 
			datfile.write(print_overlapping_segs(event, chr, start, end))
	datfile.close()

def make_genes_history(bedfn, evnts, outdir):
	outfn=os.path.join(outdir, "genes_history.dat")
	outfh=open(outfn, 'w')
	bedlines=open(bedfn, 'r').readlines()
	for bedline in bedlines:
		bede=bedmod.BedEntry(bedline)
		bede.chr="chr"+bede.chr  #need to do this because bedmod removes leading chr
		for event in evnts: 
			if event.check_overlap(bede.chr, bede.start, bede.end): 
				geffect=determine_gene_effect(bede, event)
				datout="\t".join([bede.name, geffect, str(event.prevalmean), str(event.prevalsd), str(event.likelihood)])	
				outfh.write(datout+"\n")


# In this, it is important to note that for segs, the CN value has the inverse sign, so amplificiations will have a negative cn and deletions will have a positive cn. 
def determine_gene_effect(bede, event):
	effect="none"
	for seg in event.segs: 
		if seg.seg: 
			seg.order_ends()
			if (seg.chr == bede.chr and seg.start < bede.end and seg.end > bede.start): 
				if seg.cnval <0 and effect is not "del": 
					effect= "amp"
				elif seg.cnval > 0 and effect is not "amp": 
					effect= "del"
				else: 
					effect = "amdl"
		elif seg.adj: 
			if ((seg.chr == bede.chr and (seg.start <= bede.end and seg.start >= bede.start)) 
				or (seg.chr2 == bede.chr and seg.start <= bede.end and seg.end>=bede.start)): 
				if effect is "none": 
					effect="brk"
	return effect

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
	parser.add_argument('--historystats', help='historystats.txt file')
	parser.add_argument('--genes', help='This is a bed file with gene coordinates.')
	parser.add_argument('--cutoff', help='the pvalue cutoff to use', default=0.5, type=float)
#	parser.add_argument('--plotpdf', help='generate a scatter plot of the events')
	args=parser.parse_args()
	if not os.path.exists(args.outdir): 
		os.mkdir(args.outdir)
	events=pickle.load(open(args.pevnts, 'rb'))
	myevents=[]
	for e in events: 
		e.unpack()
		if e.likelihood > args.cutoff: 
			myevents.append(e)
	sys.stderr.write("There are %d events\n" % len(myevents))
	if args.genes: 
		make_genes_history(args.genes, myevents, args.outdir)
	elif args.historyScores:
		historyScores=np.loadtxt(args.historystats, dtype='int')
		wholehistory_analysis(myevents, args.outdir, args.historyScores)
	else: 
		sys.stderr.write("Need to give either genes or historyScores.\n")


	
