#!/inside/home/common/bin/python2.7

import sys, os
import event_cycles_module as histseg
import math
import subprocess

def score_intersection_of_eventlists(eventsA, eventsB):
	tf=os.tmpfile()
	for evnt in eventsA:
		tf.write("%s\t%s\t%f\n" % (evnt.id, evnt.histlist[0], evnt.complexity_costs[0]))
	for evnt in eventsB: 
		tf.write("%s\t%s\t%f\n" % (evnt.id, evnt.histlist[0], evnt.complexity_costs[0]))
	tf.seek(0)
	tfout=os.tmpfile()
	returncode=subprocess.check_call(["sort | uniq -d | cut -f2,3 | sort -u "], stdin=tf, stdout=tfout, shell=True)			
	tf.close()
	tfout.seek(0)
	myscores=[] # a list of sets (eventid, historyid, complexity_score)
	myscore=0
	n=0
	for line in tfout: 
		(h, cost)=line.strip().split('\t')
		c=float(cost)
		myscore += math.exp(-c)
		n+=1
	return (myscore, n)

if __name__ == "__main__":	
	import argparse
	parser = argparse.ArgumentParser(description='Given a list of genes, it will cassign each pair of genes a score based on the likelihood of them both being in the same event across multiple histories.', epilog="The output file has the following columns: \n<1.geneA name>\n<2.geneB name>\n<3.Likelihood of geneA and geneB being in the same event.>\n<4.Number of events containing geneA and geneB>\n<5.number of events with geneA>\n<6. number of events with geneB>\n<7. distance between geneA and geneB>")
	parser.add_argument('bed', help='a bedfile of genes to compare')
	parser.add_argument('cnavg', help='the CN-AVG output directory for a sample')
	args=parser.parse_args()
	cnavgdir=args.cnavg
	numsims=10
	tabixfiles=range(numsims)
	totalprob=0
	for sim in xrange(numsims):
		braneyfn="%s/HISTORIES_%d.braney" % (cnavgdir, sim)
		(bsegtabixfn, badjtabixfn)=histseg.make_tabix_from_braney(braneyfn, "./") 
		totalprob+=histseg.get_total_history_prob(braneyfn)
		#sys.stderr.write("prob for %d is %s\n" % (sim, str(totalprob)))
		tabixfiles[sim]=(bsegtabixfn, badjtabixfn)

	genehash={} # keys are gene loci, values are lists of events
	gene_namehash={} # keys are gene loci, values are gene names
	for bedline in open(args.bed, 'r'):
		beddata=bedline.strip().split('\t')
		(chr, start, end)=beddata[0:3]
		locusid="%s_%s_%s" % (chr, start, end)
		gene_namehash[locusid]=beddata[3]
		overlapping_events=[]
		# compute the total likelihood
		for sim in xrange(numsims):
			(bsegtabixfn, badjtabixfn)=tabixfiles[sim]
			overevents=histseg.get_overlapping_events(chr, int(start), int(end), bsegtabixfn, badjtabixfn)
			for evnt in overevents: 
				evnt.histlist=[evnt.histlist[0]+10000*sim] 
			overlapping_events+=overevents
		if len(overlapping_events)>0:
			genehash[locusid]=overlapping_events

	# After getting all of the events that overlap each gene, compare events between genes		
	genes=genehash.keys()
	for i in xrange(len(genes)): 
		geneA=genes[i]
		for geneB in genes[i+1:]: 
			eventsA=genehash[geneA]
			eventsB=genehash[geneB]
			sys.stderr.write("working on gene %s with %d events and gene %s with %d events\n" % (geneA, len(eventsA), geneB, len(eventsB)))
			(score,N)=score_intersection_of_eventlists(eventsA, eventsB)
			likelihood=score/totalprob
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
			sys.stdout.write("%s\t%s\t%f\t%d\t%d\t%d\t%s\n" % (geneAname, geneBname, likelihood, N, len(eventsA), len(eventsB), str(distance)))


