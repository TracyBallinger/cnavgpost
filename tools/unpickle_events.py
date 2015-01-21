#!/usr/bin/env python 
import os, sys
import cnavgpost.mergehistories.event_cycles_module as histseg
import cPickle as pickle
import argparse


if __name__ == '__main__': 
	parser= argparse.ArgumentParser(description='Will unpickle a pickled events file.') 
	parser.add_argument('pickle', help='The pickled file to view')
	parser.add_argument('--historyid', type=int, help='The history you want to look at')
	parser.add_argument('--eventid', help='The id of the event you want to look at.')
	parser.add_argument('--bedfile', help='Get events that overlap these regions in the bedfile.')
	args=parser.parse_args()
	events=pickle.load(open(args.pickle, 'rb'))
	if args.eventid != None: 
		eid=args.eventid
		sys.stderr.write("event id is %s\n" % eid)
		for event in events:
			if event.id == eid: 
				sys.stdout.write(str(event))
	elif args.bedfile != None: 
		for bedline in open(args.bedfile, 'r'): 
			beddata=bedline.strip().split('\t')
			(chr, start, end) = beddata[0:3]
			for event in events: 
				event.unpack()
				if event.check_overlap(chr, int(start), int(end)): 
					sys.stdout.write(str(event))
	else: 
		for e in events: 
			e.histories=histseg.listout_ranges(e.histRanges)
		if args.historyid != None: 
			hid=args.historyid
			sys.stderr.write("hid is %d\n" % hid)
			for event in events:
				if hid in event.histories: 
					sys.stdout.write(str(event))
		else: 
			for event in events:
				sys.stdout.write(str(event))
				#sys.stdout.write(str(event.prevals) + "\n")
