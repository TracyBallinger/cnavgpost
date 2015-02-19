#!/usr/bin/env python 
import os, sys
import cnavgpost.mergehistories.event_cycles_module as histseg
import cPickle as pickle
import argparse
import numpy as np

def writedat(event, dat=False):
	if dat:  
		mystr="\t".join(map(str, [
			event.id, 
			event.determineEventType(), 
			np.mean(np.array(event.uppercosts + event.lowercosts)),
			event.likelihood, 
			event.cnval, 
			event.get_Event_length(), 
			event.numsegs,
			event.prevalmean, event.prevalsd,
			event.ordermean, event.ordersd, 
			event.numhists]
		)) + "\n"
	else: 
		mystr=str(event)
	return mystr



if __name__ == '__main__': 
	parser= argparse.ArgumentParser(description='Will unpickle a pickled events file.') 
	parser.add_argument('pickle', help='The pickled file to view')
	parser.add_argument('--historyid', type=int, help='The history you want to look at')
	parser.add_argument('--eventid', help='The id of the event you want to look at.')
	parser.add_argument('--bedfile', help='Get events that overlap these regions in the bedfile.')
	parser.add_argument('--dat', action='store_true', help='print event data.')
	args=parser.parse_args()
	if args.dat: 
		sys.stdout.write("event_id\tevent_type\tavecost\tLscore\tCNval\tlength\tnumsegs\tpvlmean\tpvlsd\tordmean\tordsd\tnumhists\n")
	events=pickle.load(open(args.pickle, 'rb'))
	if args.eventid != None: 
		eid=args.eventid
		sys.stderr.write("event id is %s\n" % eid)
		for event in events:
			if event.id == eid: 
				sys.stdout.write(writedat(event, args.dat))
	elif args.bedfile != None: 
		for bedline in open(args.bedfile, 'r'): 
			beddata=bedline.strip().split('\t')
			(chr, start, end) = beddata[0:3]
			for event in events: 
				event.unpack()
				if event.check_overlap(chr, int(start), int(end)): 
					sys.stdout.write(writedat(event, args.dat))
	else: 
		for e in events: 
			e.histories=histseg.listout_ranges(e.histRanges)
		if args.historyid != None: 
			hid=args.historyid
			sys.stderr.write("hid is %d\n" % hid)
			for event in events:
				if hid in event.histories: 
					sys.stdout.write(writedat(event, args.dat))
		else: 
			for event in events:
				sys.stdout.write(writedat(event, args.dat))
				#sys.stdout.write(str(event.prevals) + "\n")
