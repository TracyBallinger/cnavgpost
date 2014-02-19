#!/inside/home/common/bin/python2.7

import sys, os
import argparse
import gzip
import event_cycles_module as histseg

# Input: the filename of a .braney file
# Output: a list of events, where each event is a list of braney segments and adjacencies 
def write_events(braneyfn):
	braneyf=gzip.open(braneyfn, 'rb')
	eventid=""
	prevhistid=-1
	allcosts=[]
	eventlines=[]
	for braneyline in braneyf: 
		braneyseg=histseg.Braney_seg(braneyline)
		if braneyseg.historyid != prevhistid: 
			allcosts.append(braneyseg.complexity)
			prevhistid=braneyseg.historyid
		if eventid=="": 
			eventid=braneyseg.ptrid  #the pointer id is unique between histories, even if everything else is the same
			eventlines.append(braneyseg)
		elif braneyseg.ptrid == eventid:
			eventlines.append(braneyseg)
		elif braneyseg.ptrid != eventid: 
			myevent=histseg.Event(eventlines)
			myevent.remove_dup_adj()
			myevent.get_hranges()
			myevent.compute_likelihood(1)
			sys.stdout.write("%s" % (str(myevent)))
			eventlines=[braneyseg]
			eventid=braneyseg.ptrid
	return (allcosts)

if __name__ == '__main__': 	
	parser = argparse.ArgumentParser(description='Will essentially reformat braney files to have one event (rearrangment cycle) per line.') 
	parser.add_argument('braney', help='a .braney file', default='-')
	parser.add_argument('--costs', help='a text file to write all of the costs to')
	args=parser.parse_args()
	braneyfn=args.braney
	allcosts=write_events(braneyfn)
	if args.costs:  
		if args.costs == "stderr":
			costf=sys.stderr
		else:
			costf=open(args.costs, 'w')
		for cost in allcosts: 
			costf.write("%d\n" % (cost))

