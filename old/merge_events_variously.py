#!/inside/home/common/bin/python2.7

import sys, os
import event_cycles_module as histseg 
import copy 

def merge_by_cnsign_or_genes(events, bygene, byCNsign): 
	eventa=events[0]
	myevent=copy.deepcopy(eventa)
	cnvals=[eventa.cnval]
	pvals=[eventa.prevalmean]
	scores=[eventa.likelihood]
	myevent=copy.deepcopy(eventa)
	for eventb in events[1:]:
		if (bygene and same_genes(myevent, eventb) or (byCNsign and same_cnsign(myevent, eventb): 
			cnvals.append(eventb.cnval)
			pvals.append(eventb.prevalmean)
			scores.append(eventb.likelihood)
		else: 
			myevent.cnval=mean(cnvals)
			myevent.prevalmean=mean(pvals)
			myevent.likelihood=sum(scores)
			fevents.append(myevent)
			myevent=copy.deepcopy(eventb)
	# add the last event 
	myevent.cnval=mean(cnvals)
	myevent.prevalmean=mean(pvals)
	myevent.likelihood=sum(scores)
	fevents.append(myevent)
	return fevents

def same_genes(eventa, eventb): 
	return if (eventa.genes != "None" and eventa.genes == eventb.genes)

def same_cnsign(eventa, eventb): 
	return if (eventa.segstr == eventb.segstr and 
		(eventa.cnval>0 and eventb.cnval>0) or (eventa.cnval <0 and eventb.cnval<0))
 

if __name__== '__main__': 
	import argparse
	parser=argparse.ArgumentParser(description='Will combine events in an .evnts file based on various criteria')
	parser.add_argument('events', help='the .evnts file')
	parser.add_argument('--genes', help='combine by the genelist names (This will be the last column in the evnts file.)', action='store_true')
	parser.add_argument('--CNsign', help='combine same cycles if they have the same CN sign', action='store_true')
	args=parser.parse_args()
	#read in all of the events and sort them by prevalence (timing)
	eventfile = open(args.events, 'r')
	allevents=[]
	for eline in eventfile: 
		evnt=histseg.Event(eline)
		allevents.append(evnt)
	sortedevents=sorted(allevents, key=lambda x: (x.prevalmean))
	
	if args.genes: 	
		tmpevents=merge_by_genes(sortedevents)
		sortedevents=tmpevents
	if args.CNsign:
		tmpevents=merge_by_cnsign(sortedevents)
		sortedevents=tmpevents
	
	for evnt in sortedevents: 
		sys.stdout.write(str(evnt))
	


	  
