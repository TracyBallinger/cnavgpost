#!/inside/home/common/bin/python2.7

import sys, os
import event_cycles_module as histseg
import glob
import re
bindir="/inside/home/tballing/cnavg-study/cnavgmbin"

if __name__ == '__main__': 
	import argparse
	parser = argparse.ArgumentParser(description='Will essentially reformat braney files to have one event (rearrangment cycle) per line, and will assign each event a likelihood score based on the history complexity values')
	parser.add_argument('cnavg', help='the CN-AVG output directory for a sample')
	args=parser.parse_args()
	cnavgdir=args.cnavg
	totalprob=0
	allevents=[]
	braneyfiles=glob.glob(cnavgdir+"/"+"HISTORIES_*.braney")
	for braneyfn in braneyfiles:
		sim=int(re.match(".*HISTORIES_(\d+)\.braney", braneyfn).group(1))
#		sys.stderr.write("braneyfn: %s, sim: %s\n" % (braneyfn, sim))
		totalprob+=histseg.get_total_history_prob(braneyfn)
		events=histseg.make_events_from_braneyfn(braneyfn)
		dbevent=events[0]
		sortedevents=sorted(events, key=lambda x: (x.segstr, x.cnval))
		uniqueevents=histseg.unique_c_events_sorted_list(sortedevents)
		dbevent=events[0]
		dbevent=uniqueevents[0]
		for evnt in uniqueevents:
			for i in xrange(evnt.numhists):
				evnt.histories[i] = evnt.histories[i] + sim*10000 
		allevents+=uniqueevents
		sys.stderr.write("totalprob: %s, num_events: %d, filtered: %d\n" % (str(totalprob), len(events), len(allevents)))
	# do a final merge of events across different simulations
	sortedevents=sorted(allevents, key=lambda x: (x.segstr, x.cnval))
	finalevents=histseg.unique_c_events_sorted_list(sortedevents)
	for evnt in finalevents: 
		evnt.likelihood=histseg.compute_likelihood(evnt.costs, totalprob)
		sys.stdout.write("%s" % (str(evnt)))

