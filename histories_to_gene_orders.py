#!/inside/home/common/bin/python2.7

import sys, os
import pickle
import likelihood_score_gene_pairs as genemod
import event_cycles_module as histseg


if __name__ == "__main__": 
	import argparse
	parser=argparse.ArgumentParser(description='')
	parser.add_argument('pevnts', help='a .pevnts file for a sample')
	parser.add_argument('annotation', help='an annotation file of the events (see annotate_event.py)')
	
	parser.add_argument('--totalp', type=float)
	args = parser.parse_args()
	allevents=pickle.load(open(args.pevnts, 'rb'))
	sys.stderr.write("loaded all the events: %d\n" % (len(allevents)))
	annotations=genemod.read_in_annotations(args.annotation)
	sys.stderr.write("Read in all the annotations: %d\n" % (len(annotations)))
	geneEvents=genemod.create_gene_events_hash(allevents, annotations)  #key: a gene name, value: a list of event indexes
	sys.stderr.write("geneEvents: %d\n" % (len(geneEvents)))
	if not args.totalp: 
		totalp = histseg.get_total_likelihood(allevents)
	else: 
		totalp = args.totalp
	sys.stderr.write("totalp: %s\n" % (str(totalp)))

	for gene in geneEvents.keys():
		myevents=[allevents[i] for i in geneEvents[gene]]
		mergedEvent=histseg.merge_events(myevents)
		mergedEvent.compute_timing_wmeansd()
		mergedEvent.likelihood = histseg.compute_likelihood(mergedEvent.costs, totalp)
		mergedEvent.numhists=len(mergedEvent.histories)
		self = mergedEvent
		mystr=("%s\t%f\t%s\t%f\t%s\t%s\t%d\n" % (gene, self.ordermean, str(self.ordersd), self.prevalmean, str(self.prevalsd), str(self.likelihood), self.numhists))
		sys.stdout.write(mystr)	
		#sys.stderr.write(mystr)	
		
