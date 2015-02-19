#!/usr/bin/env python 

import sys, os
import cnavgpost.mergehistories.event_cycles_module as histseg
from cnavgpost.genehistory.bedFileModule import *  
import argparse
import cPickle as pickle 
import subprocess, re
import numpy as np

def filter_events_by_cn_overlap(events, bedfn, perc, cbs=True, datout=False):
	if cbs: 
		mystr=subprocess.check_output("sed 1d %s | cut -f2-4" % bedfn, shell=True)
		bedlines=mystr.split("\n")
	else: 
		bedlines=open(bedfn, 'r').readlines()
	bedentries=[]
	for line in bedlines: 
		if line!="": 
			bedentries.append(BedEntry(line))
	bedentries=sorted(bedentries, key=lambda x: (x.chr, x.start, x.end))
	myevents=[]
	for event in events:
		if datout:
			flag=0 
			if event_overlaps_CN(event, bedentries, perc):
				flag=1
			myevents.append(printEventData(event, flag))
		elif event_overlaps_CN(event, bedentries, perc):
			myevents.append(event)
	return(myevents)

def printEventData(event, info): 
	return "%s\t%d\t%f\t%d\t%s\n" % (event.id, event.determineEventType(), event.likelihood, event.numhists, str(info))

def event_overlaps_CN(event, bedentries, perc):
	if event.segs==[]:
		event.make_segs_from_str()
	keep=True
	for seg in event.segs: 
		if seg.seg:
			region=BedEntry("%s\t%d\t%d\n" % (seg.chr, seg.start, seg.end)) 
			poverlaps=percent_bed_overlaps(bedentries, region, perc)
#			sys.stderr.write("seg %s:%d-%d overlaps: %s\n" % (seg.chr, seg.start, seg.end, str(poverlaps)))
			if len(poverlaps) >0 and (max(poverlaps) < perc): 
				keep=False
	return keep

# find what the overlap is for the given region in the bedentries		
def percent_bed_overlaps(bedentries, region, perc):
	i=len(bedentries)/2
	bedentry=bedentries[i]
	overlap=bedentry.overlap(region)
	myvals=[]
	if overlap>0: 
		poverlap=min(float(overlap)/(bedentry.end-bedentry.start+1), 
					float(overlap)/(region.end-region.start+1))
		myvals.append(poverlap) 
	if bedentry.comes_before(region) and (i+1)< len(bedentries):
		myvals+=percent_bed_overlaps(bedentries[(i+1):], region, perc)
	elif i > 0: 
		myvals+=percent_bed_overlaps(bedentries[:i], region, perc)
	return myvals


## This will filter out events that cancel eachother out (ie, an amplification immediately followed by a deletion). 
def filter_fleeting_events(events, histScores, totalp, fcutoff=0, datout=False): 
	for e in events: 
		(newstr, sign) = histseg.remove_signs_from_segstr(e.segstr)
		e.cnval=sign*e.cnval 
		e.segstr=newstr
	sevents=sorted(events, key=lambda x: (x.segstr, x.prevalmean))
	e1=sevents[0]
	sameEs=[e1]
	filtered=[]
	for e2 in sevents[1:]: 
		if e1.segstr == e2.segstr:
			sameEs.append(e2)
		else: 
#			filtered+=remove_canceling_histories(sameEs, histScores, totalp, datout)	
			filtered+=add_canceling_events(sameEs, datout)	
			sameEs=[e2]
			e1=e2
#	filtered+=remove_canceling_histories(sameEs, histScores, totalp, datout)	
	filtered+=add_canceling_events(sameEs, datout)	
	return filtered

def add_canceling_events(events, datout=False): 
	for i in xrange(len(events)): 
		e1=events[i]
		for e2 in events[i:]: 
			if (e1.segstr == e2.segstr) and ((e1.cnval + e2.cnval) == 0): 
		#		sys.stderr.write("%s is cancelling %s, %f, %f\n" % (e1.id, e2.id, e1.likelihood, e2.likelihood))
				if e1.likelihood>e2.likelihood: 
					e1.likelihood=e1.likelihood-e2.likelihood
					e2.likelihood=0
				else: 
					e2.likelihood=e2.likelihood-e1.likelihood
					e1.likelihood=0
		#		sys.stderr.write("after %s is cancelling %s, %f, %f\n" % (e1.id, e2.id, e1.likelihood, e2.likelihood))
	filtered=[]
	if datout: 
		for e in events: 
			filtered.append(printEventData(e, str(e.likelihood)))
		#	sys.stderr.write("event %s %f\n" % (e.id, e.likelihood))
	else: 
		for e in events: 
			filtered.append(e)
	return filtered 

def remove_canceling_histories(events, histScores, totalp, datout=False): 
	myinfo=[]
	if datout: 
		for e in events: 
			myinfo.append("%f\t%d" % (e.likelihood, len(e.histories)))
	for i in xrange(len(events)): 
		e1=events[i]
		for e2 in events[i:]:
			if (e1.segstr == e2.segstr) and ((e1.cnval + e2.cnval) == 0): 
				# only keep the histories that are unique to each event. 
				histseg.cancel_Event_data(e1,e2)
			#	sys.stderr.write("%s is cancelling %s\n" % (e1.id, e2.id))
	filtered=[]
	if datout: 
		for (e, info) in zip(events, myinfo):
			#e.compute_timing_wmeansd(histScores)
			#e.compute_likelihood(histScores, totalp)
			filtered.append(printEventData(e, info)) 
	else: 
		for e in events:
			if  len(e.histories) >0:
				e.compute_timing_wmeansd(histScores)
				e.compute_likelihood(histScores, totalp)
				filtered.append(e)
	return filtered 

def filter_zero_cost_events(events, histScores, totalp, dat): 
	myinfo=[]
	for e in events: 
		numzerol=sum(np.array(e.lowercosts)==0)
		numzeroh=sum(np.array(e.uppercosts)==0)
		myinfo.append(printEventData(e, "%f\t%d\t%d\t%d" % (e.likelihood, len(e.histories), numzerol, numzeroh)))
	return(myinfo)	

def main(args): 
	allevents=pickle.load(open(args.pevnt, 'rb'))
	historyScores=np.loadtxt(open(args.histstats, 'r'), dtype=int)
	totalp=histseg.compute_totalp(historyScores) 
	sys.stderr.write("Begin, %d events\n" % (len(allevents)))
	for e in allevents:
		e.histories=histseg.listout_ranges(e.histRanges)
	if args.CNoverlap: 
		allevents=filter_events_by_cn_overlap(allevents, args.bed, args.perc, args.cbs, args.dat)
		sys.stderr.write("After CNoverlap, %d events\n" % (len(allevents)))
	if args.fleeting: 
		allevents=filter_fleeting_events(allevents, historyScores, totalp, fcutoff=0, datout=args.dat) 
		sys.stderr.write("After fleeting, %d events\n" % (len(allevents)))
	if args.zeroCost: 
		allevents=filter_zero_cost_events(allevents, historyScores, totalp, args.dat)
		sys.stderr.write("After zeroCost, %d events\n" % (len(allevents)))
	fout=open(args.out, 'w')
	if args.dat:
		for e in allevents: 
			fout.write(e) 
	else: 
		for e in allevents:
			e.compute_timing_wmeansd(historyScores)
			e.compute_likelihood(historyScores, totalp)
		pickle.dump(allevents, fout, pickle.HIGHEST_PROTOCOL)  


if __name__ == '__main__': 
	parser=argparse.ArgumentParser("This will filter events based on the options selected.")
	parser.add_argument('pevnt', help='the pickled event file (.pevnt)')
	parser.add_argument('histstats', help='historystats.txt file')
	parser.add_argument('out', help='The name of the output file.')
	parser.add_argument('--dat', action='store_true', help='print data for the events rather than the filtered events.')
	parser.add_argument('--fleeting', action='store_true', help='Filter events based on being fleeting.')
	parser.add_argument('--zeroCost', action='store_true', help='Filter events based on having a cost of zero.')
	parser.add_argument('--CNoverlap', action='store_true', help='Filter events based on CN overlap')
	parser.add_argument('--bed', help='A bed file of the CN changes')
	parser.add_argument('--cbs', help='The file is a CBS_OUT file of the CN changes', action='store_true')
	parser.add_argument('--perc', help='A percentage overlap', default=0.9, type=float)
	args=parser.parse_args()
	main(args)

