#!/inside/home/common/bin/python2.7

import sys, os
import event_cycles_module as histseg
import math
import subprocess
import tempfile

bindir="/inside/home/tballing/cnavg-study/cnavgmbin"

def print_event_details(evnt, datfn): 
	if os.path.exists(datfn):
		datfile=open(datfn, 'a')
	else: 
		datfile=open(datfn, 'w')
	for i in xrange(len(evnt.prevals)): 
		datfile.write("%f\t%f\t%f\n" % (evnt.prevals[i], evnt.histories[i], evnt.costs[i]))
	datfile.close()

if __name__ == '__main__': 
	import argparse
	parser = argparse.ArgumentParser(description='Print out the prevalence values, historyids, and complexity costs for the given event.')
	parser.add_argument('cnavg', help='the CN-AVG output directory for a sample')
	parser.add_argument('evnts', help='A list of a events to print all of the data for.')
	parser.add_argument('datdir', help='directory to output all the values for prevalence to get distribution information. A .dat file will be created for each event, with the name of the file being the id of the event.')
	args=parser.parse_args()
	cnavgdir=args.cnavg
	if not os.path.exists(args.datdir): 
		os.makedirs(args.outdir)
	eventstokeep=[]
	datfnhash={}
	eventf=open(args.evnts, 'r')
	for eline in eventf: 
		evnt=histseg.Event(eline)
		eventstokeep.append(evnt)
		hashkey="%f+%s" % (evnt.cnval, evnt.segstr) 
		datfnhash[hashkey]="%s.dat" % (evnt.id)

	sys.stderr.write("number of events to keep: %d" % (len(eventstokeep)))
	for sim in xrange(10):
		sys.stderr.write("working on sim: %d\n" % (sim))
		braneyfn="%s/HISTORIES_%d.braney" % (cnavgdir, sim)
		events=histseg.make_events_from_braneyfn(braneyfn)
		for evnt in events: 
			if evnt in eventstokeep:
				hashkey="%f+%s" % (evnt.cnval, evnt.segstr) 
				datfn=args.datdir + "/" + datfnhash[hashkey]
				print_event_details(evnt, datfn) 	

