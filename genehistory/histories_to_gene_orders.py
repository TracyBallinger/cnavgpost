#!/usr/bin/env python

import sys, os
import cpickle
import cnavgpost.genehistory.likelihood_score_gene_pairs as genemod
import cnavgpost.mergehistories.event_cycles_module as histseg

def main(allevents, annotationfn, totalp, outputfh): 
	annotations=genemod.read_in_annotations(annotationfn)
	sys.stderr.write("Read in all the annotations: %d\n" % (len(annotations)))
	geneEvents=genemod.create_gene_events_hash(allevents, annotations)  #key: a gene name, value: a list of event indexes
	sys.stderr.write("geneEvents: %d\n" % (len(geneEvents)))
	if not totalp: 
		totalp = histseg.get_total_likelihood(allevents)
	sys.stderr.write("totalp: %s\n" % (str(totalp)))
	for gene in geneEvents.keys():
		myevents=[allevents[i] for i in geneEvents[gene]]
		sys.stderr.write("Working on %s with %d events\n" % (gene, len(myevents)))
		mergedEvent=histseg.merge_events(myevents)
		mergedEvent.compute_timing_wmeansd()
		mergedEvent.likelihood = histseg.compute_likelihood(mergedEvent.costs, totalp)
		mergedEvent.numhists=len(mergedEvent.histories)
		self = mergedEvent
		mystr=("%s\t%f\t%s\t%f\t%s\t%s\t%d\n" % (gene, self.ordermean, str(self.ordersd), self.prevalmean, str(self.prevalsd), str(self.likelihood), self.numhists))
		outputfh.write(mystr)	
		sys.stderr.write(mystr)	
		

if __name__ == "__main__": 
	import argparse
	parser=argparse.ArgumentParser(description='')
	parser.add_argument('pevnts', help='a .pevnts file for a sample')
	parser.add_argument('annotation', help='an annotation file of the events (see annotate_event.py)')
	
	parser.add_argument('--totalp', type=float)
	args = parser.parse_args()
	allevents=pickle.load(open(args.pevnts, 'rb'))
	sys.stderr.write("loaded all the events: %d\n" % (len(allevents)))
	main(allevents, args.annotation, args.totalp, sys.stdout)

