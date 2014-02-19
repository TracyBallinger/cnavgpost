#!/inside/home/common/bin/python2.7

import sys, os
import event_cycles_module as histseg 
from history_node_links_module import * 
import pickle

def add_event_link_options(parser): 
	parser.add_argument('--cnavg', help='the CN-AVG output directory for a sample')
	parser.add_argument('--events', help='the file to write events to')  
	parser.add_argument('--links', help='the file to write links to')
	parser.add_argument('--inpickle', help='read the evnts from a pickled file.') 
	parser.add_argument('--outpickle', help='pickle the evnts and write them to this file.') 
	parser.add_argument('--totalp', help='total probability of histories', type=float)
	parser.add_argument('--binwidth', help='the multiplier for each history id to distinguish independent simulations.', type=float)

def score_and_link_cycles(args):
	if args.binwidth: 
		histseg.Global_BINWIDTH=args.binwidth
	totalprob=1
	if args.totalp: 
		totalprob=args.totalp
	allevents=[]	
	if args.cnavg: 
		sys.stderr.write("using cnavg dir: %s\n" % (args.cnavg))
		allevents=histseg.get_events_from_cnavgdir(args.cnavg)
		if not args.totalp: 
			totalprob = histseg.get_total_likelihood(allevents)
	elif args.inpickle and os.path.isfile(args.inpickle): 
		sys.stderr.write("using pickle file\n")
		allevents=pickle.load(open(args.inpickle, 'rb'))
	sys.stderr.write("there are %d events\n" % (len(allevents)))
	if args.outpickle: 
		eventfile= open(args.outpickle, 'wb')
		pickle.dump(allevents, eventfile, pickle.HIGHEST_PROTOCOL)
	if args.events: 
		eventfile=open(args.events, 'w')
		for evnt in allevents: 
			eventfile.write("%s" % (str(evnt)))
	if args.links: 
		eventlinks = link_events_by_order_within_histories(allevents)
		linkfile=open(args.links, 'w')
		for link in eventlinks: 
			link.likelihood=histseg.compute_likelihood(link.costs, totalprob)
			linkfile.write("%s" % (str(link)))


if __name__ == '__main__': 
	import argparse
	parser = argparse.ArgumentParser(description='Will essentially reformat braney files to have one event (rearrangment cycle) per line, and will assign each event a likelihood score based on the history complexity values.  Will also create links between cycles that co-occur in the same histories.')
	add_event_link_options(parser) 
	args=parser.parse_args()
	score_and_link_cycles(args)

