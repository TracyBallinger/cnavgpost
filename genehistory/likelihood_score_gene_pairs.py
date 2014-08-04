#!/inside/home/common/bin/python2.7
import sys, os
import event_cycles_module as histseg
import pickle, pysam 

class GenePair: 
	def __init__(self, geneA, geneB, eventsA, eventsB, annotations):
		if 1==1: 
			self.geneA = geneA
			self.geneB = geneB
			self.histories=[]
			self.costs=[]
			self.count=0
			self.countsA=0
			self.countsB=0
			self.countsABplus=0
			self.likelihood=0
			self.distance=0
			self.scoreA=0
			self.scoreB=0
			score_gene_pair(self, geneA, geneB, eventsA, eventsB, annotations)

	def __str__(self): 
		mystr=("%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\n" % (self.geneA, self.geneB, str(self.likelihood), str(self.scoreA), str(self.scoreB), self.count, self.countsA, self.countsB, self.count-self.countsABplus))
		return mystr

	def adjust_likelihoods(self, totalprob): 
		self.likelihood = self.likelihood/totalprob
		self.scoreA = self.scoreA/totalprob
		self.scoreB = self.scoreB/totalprob

def merge_histories_costs(events):
	histories=events[0].histories
	costs=events[0].costs
	for event in events[1:]: 
		addedi = [i for i in range(len(event.histories)) if event.histories[i] not in histories]
		addedhistories=[event.histories[i] for i in addedi]
		histories += addedhistories
		addedcosts=[event.costs[i] for i in addedi]
		costs += addedcosts
	return((histories, costs))
				
def	score_gene_pair(self, geneA, geneB, eventsA, eventsB, annotations):
	(historiesA, costsA) =merge_histories_costs(eventsA)
	(historiesB, costsB) =merge_histories_costs(eventsB)
	eventsAB=[]
	for event in eventsA:
		genelist = annotations[event.id]
		if geneA in genelist and geneB in genelist: 
			eventsAB.append(event)
			if len(genelist) >2: # there are other genes besides A and B
				self.countsABplus+=1
	(self.hists, self.costs) = merge_histories_costs(eventsAB)
	self.count = len(self.hists)
	self.countsA =len(historiesA)
	self.countsB = len(historiesB)
	self.scoreA = histseg.compute_likelihood(costsA, 1)
	self.scoreB = histseg.compute_likelihood(costsB, 1)
	self.likelihood= histseg.compute_likelihood(self.costs, 1)
	
def read_in_annotations(annotationsfile): 
	myannotations={}
	for line in open(annotationsfile): 
		(eventid, genes) = line.strip().split('\t')
		myannotations[eventid]=genes.split(',')
	return myannotations

def create_gene_events_hash(events, annotations): 
	geneEvents={}
	for i in xrange(len(events)):
		event=events[i] 
		genelist=annotations[event.id]
		if genelist[0] != "None":
			for gene in genelist: 
				if gene in geneEvents.keys(): 
					geneEvents[gene].append(i)
				else: 
					geneEvents[gene] = [i]
	return geneEvents

def likelihood_score_gene_pairs(allevent, annotations, tabixfn): 
	geneEvents=create_gene_events_hash(allevents, annotations) # key: a gene name, value: a list of event indexes for the events with geneX
	sys.stderr.write("Hashed in all the annotations: %d\n" % (len(geneEvents)))
	eventi=0
	myGenepairs=[]
	mygeneunpairs=[]
	pairIDs=[]
	allhistoryids=[]
	allcosts=[]
	while eventi < len(allevents):
		sys.stderr.write("working on event %d\n" % (eventi))
		myevent=allevents[eventi]
		for i in xrange(len(myevent.histories)):
			hid=myevent.histories[i]
			if hid not in allhistoryids: 
				allhistoryids.append(hid)
				allcosts.append(myevent.costs[i])
		genes=annotations[myevent.id]
		if genes[0] != "None":
			for ia in xrange(len(genes)): 
				geneA = genes[ia]
				for ib in xrange(ia+1, len(genes)): 
					geneB = genes[ib]
					genepairID="%s,%s" % (geneA, geneB)
				#	sys.stderr.write("working on %s\n" % (genepairID))
					if genepairID not in pairIDs:
						eventsA=[allevents[i] for i in geneEvents[geneA]]
						eventsB=[allevents[i] for i in geneEvents[geneB]]
				 		mypair = GenePair(geneA, geneB, eventsA, eventsB, annotations)
						pairIDs.append(genepairID)
						myGenepairs.append(mypair)
		eventi+=1
	totalp = histseg.compute_likelihood(allcosts, 1)
	sys.stderr.write("totalp: %s\n" % (str(totalp)))
	if tabixfn: 
		pair_distances = find_distance_between_genes(pairIDs, args.tabix) 
	for pair in myGenepairs:
		if tabixfn: pair.distance = pair_distances[pair.geneA+pair.geneB] 
		pair.adjust_likelihoods(totalp)
		sys.stdout.write(str(pair))

def find_distance_between_genes(myGenepairs, tabixfn): 
	mytabix = pysam.Tabixfile(tabixfn, 'r')
	locs = {}
	for loc in mytabix.fetch(): 
		(chr, start, end, name) = loc.split('\t')
		loc[name] = (chr, start, end)
	distances={}
	for pair in myGenepairs: 
		(geneA, geneB) = pair.split(',')
		(chra, starta, enda)=loc[geneA]
		(chrb, startb, endb)=loc[geneB]
		(starta, enda, startb, endb) = map(int, (starta, enda, startb, endb))
		distance="NA"
		if chra == chrb:
			if enda < startb:
				distance=startb-enda+1
			elif endb < starta:
				distance=starta-endb+1
			else:
				distance=0 # the genes' coordinates overlap 
		distances[pair]=distance
	return distances
		

if __name__ == "__main__":
	import argparse
	parser=argparse.ArgumentParser(description='Given a list of genes and an .evnts file for a sample, will assign a likelihood score L=sum(l(Events with geneA and geneB))/sum(l(Events with geneA))+sum(l(events with geneB)).', epilog="The output file has the following columns: \n<1.geneA name>\n<2.geneB name>\n<3.total likelihood of geneA and geneB being in the same event.>\n<4. Likelihood of events with geneA>\n<5. Likelihood of events with geneB>\n<6.Number of events containing geneA and geneB>\n<7.number of events with geneA>\n<8. number of events with geneB>\n<9. distance between geneA and geneB>\n<10. Number of events with just geneA and geneB>")
	parser.add_argument('pevnts', help='an .pevnts file for a sample')
	parser.add_argument('annotation', help='a annotation file of the events (see annotate_events.py)')
	parser.add_argument('--tabix', help='a tabix file with the gene annotations that were used, for calculated distances between genes')
	args=parser.parse_args()
	allevents=pickle.load(open(args.pevnts, 'rb'))
	sys.stderr.write("loaded all the events: %d\n" % (len(allevents)))
	annotations=read_in_annotations(args.annotation) #annotations will be a hash: key: eventid, value: a list of gene names
	sys.stderr.write("Read in all the annotations: %d\n" % (len(annotations)))
	myGenepairs=likelihood_score_gene_pairs(allevents, annotations, args.tabix)

