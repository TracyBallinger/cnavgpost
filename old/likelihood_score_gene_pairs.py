#!/inside/home/common/bin/python2.7
import sys, os
import event_cycles_module as histseg

class GenePair_Event: 
	def __init__(self, event, genepair): 
		if isinstance(event, histseg.Event): 
			geneA=genepair[0]
			geneB=genepair[1]
			histories=event.histories
			costs=event.costs
			count=1
			countsA=1
			countsB=1
			countsABplus=0
			likelihood=0
			distance=0
	def __str__(self): 
		mystr=("%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\n" % (self.geneA, self.geneB, str(self.likelihood), str(self.geneAscore), str(self.geneBscore), self.count, self.countA, self.countB, self.count-self.countsABplus))
		return mystr

def compute_likelihood_sum(events):
	l=0
	N=0
	for event in events:
		l+=event.likelihood
		N+=1
	return(l, N)

def unique_d_events(eventsa, eventsb): 
	bothevents=[]
	for eventa in eventsa: 
		for eventb in eventsb: 
			if eventa==eventb:
				bothevents.append[eventa]
	return bothevents

if __name__ == "__main__":
	import argparse
	parser=argparse.ArgumentParser(description='Given a list of genes and an .evnts file for a sample, will assign a likelihood score L=sum(l(Events with geneA and geneB))/sum(l(Events with geneA))+sum(l(events with geneB)).', epilog="The output file has the following columns: \n<1.geneA name>\n<2.geneB name>\n<3.total likelihood of geneA and geneB being in the same event.>\n<4. Likelihood of events with geneA>\n<5. Likelihood of events with geneB>\n<6.Number of events containing geneA and geneB>\n<7.number of events with geneA>\n<8. number of events with geneB>\n<9. distance between geneA and geneB>\n<10. Number of events with just geneA and geneB>")
	parser.add_argument('bed', help='a bedfile of genes to compare')
	parser.add_argument('pevnts', help='an .pevnts file for a sample')
	args=parser.parse_args()
	gene_namehash={} #keys: a genomic locus, values: a gene name
	for bedline in open(args.bed, 'r'):
		beddata=bedline.strip().split('\t')
		(chr, start, end)=beddata[0:3]
		locusid="%s_%s_%s" % (chr, start, end)
		gene_namehash[locusid]=beddata[3]
		(start, end) = (int(start), int(end))
		sys.stderr.write("args.evnts: %s\n" % (args.evnts))
		events=histseg.get_overlapping_events(chr, start, end, args.evnts)
		if len(events)>0: 
			genehash[locusid]=overlapping_events

	#After getting all of the events that overlap each gene, compare events between genes
	genes=genehash.keys()
	for i in xrange(len(genes)):
		geneA=genes[i]
		for geneB in genes[i+1:]:
			eventsA=genehash[geneA]	
			eventsB=genehash[geneB]
			(lEa, Na) = compute_likelihood_sum(eventsA)	
			(lEb, Nb) = compute_likelihood_sum(eventsB)	
			eventsAB=histseg.unique_d_events(eventsA, eventsB)
			(lEab, Nab) = compute_likelihood(eventsAB)
			finalscore=lEab/(lEa+lEb)
			# find the distance between the genes
			(chra, starta, enda)=geneA.split('_')
			(chrb, startb, endb)=geneB.split('_')
			(starta, enda, startb, endb) = map(int, (starta, enda, startb, endb))
			distance="NA"
			if chra == chrb:
				if enda < startb:
					distance=startb-enda+1
				elif endb < starta:
					distance=starta-endb+1
				else:
					distance=0 # the genes' coordinates overlap 
			geneAname=gene_namehash[geneA]
			geneBname=gene_namehash[geneB]
			sys.stdout.write("%s\t%s\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%s\n" % (geneAname, geneBname, finalscore, lEab, lEa, lEb, Nab, Na, Nb, str(distance)))
