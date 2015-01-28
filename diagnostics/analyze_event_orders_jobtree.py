#!/usr/bin/env python 

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from optparse import OptionGroup
from optparse import OptionParser
from sonLib.bioio import logger
from sonLib.bioio import system
import os, sys, glob, re
import subprocess
import cPickle as pickle
import cnavgpost.mergehistories.event_cycles_module as histseg
import cnavgpost.diagnostics.get_event_order_counts as ordcnts
import cnavgpost.simulations.correct_simulation_events as correctevents 
import cnavgpost.cnavg_post_analysis_jobtree as cpaj
import numpy as np

class Setup(Target):
	def __init__(self, options):
		Target.__init__(self)
		self.options=options

	def run(self):
		opts=self.options
		self.logToMaster("setting up...")
		mylines=open(opts.samplelist, 'r').readlines()
		for line in mylines:
			(cnavgpost, sampleid) = line.strip().split('\t')[:2]
			self.logToMaster("working on %s in %s" % (sampleid, cnavgpost))
			statsfn=os.path.join(cnavgpost, "historystats.txt")
			outputdir=cnavgpost
			pevntsfile=os.path.join(outputdir, sampleid + ".pevnts")	
			self.addChildTarget(RunEventOrderAnalysis(sampleid, cnavgpost, statsfn, opts.edges, opts.cutoff, opts.simulation))	

class RunEventOrderAnalysis(Target): 
	def __init__(self, sampleid, cnavgpost, statsfn, useEdges, cutoff, sim): 
		Target.__init__(self)
		self.sampleid=sampleid
		self.cnavgpost=cnavgpost
		self.statsfn=statsfn
		self.useEdges=useEdges
		self.cutoff=cutoff
		self.simulation=sim
	
	def run(self):
		pvalcutoff=self.cutoff 
		outdir=self.cnavgpost
		useEdges=self.useEdges
		if useEdges: 
			pevntsfile=os.path.join(outdir, "%s.pedgs" % self.sampleid)
			outfh=open(os.path.join(outdir, "edges_ordcnts.dat"), 'w')
		else: 
			pevntsfile=os.path.join(outdir, "%s.pevnts" % self.sampleid)
			outfh=open(os.path.join(outdir, "events_ordcnts.dat"), 'w')
		# filter the events
		events=pickle.load(open(pevntsfile, 'rb'))
		myevents=[]
		if pvalcutoff >0:
			for e in events:
				if e.likelihood > pvalcutoff:
					myevents.append(e)
		else:
			myevents=events
		self.addChildTarget(ordcnts.get_events_order_counts(myevents, outfh, self.simulation))
		outfh.close()
		if False: #I'm still debugging this so it's not included
			if useEdges:
				outfh1=open(os.path.join(outdir, "edges_earlycnts.dat"), 'w')
			else:
				outfh1=open(os.path.join(outdir, "event_earlycnts.dat"), 'w')
			outfn2=os.path.join(outdir, "histlengths.dat")
			self.addChildTarget(ordcnts.count_earlylate_with_correction(myevents, historyScores, outfh1, outfn2))
			outfh1.close()


class PrintEventCountsPerHistory(Target): 
	def __init__(self, sampleid, cnavgpost, opts): 
		Target.__init__(self)
		self.opts=opts
		self.sampleid=sampleid
		self.cnavgpost=cnavgpost

	def run(self): 
		opts=self.opts
		braneyfiles=glob.glob(os.path.join(opts.cnavgout, self.sampleid, "HISTORIES_*.braney"))
		self.logToMaster("%s, braneyfiles: %s\n" % (opts.cnavgout, braneyfiles))
		sys.stderr.write("%s, braneyfiles: %s\n" % (opts.cnavgout, braneyfiles))
		runlens=[]
		myarrays=[]
		mysims=[]
		for bfile in braneyfiles: 
			sim=int(re.match(".*/HISTORIES_(\d+).braney", bfile).group(1))
			mysims.append(sim)
			histlengths=subprocess.check_output("zcat %s | grep ^A | cut -f10,14 | sort -u | cut -f1 | sort -n | uniq -c | awk '{print $2\"\\n\"$1}'" % (bfile), shell=True)
			myarray=np.fromstring(histlengths, dtype=int, sep='\n')
			runlen=myarray.size/2
			myarray=np.reshape(myarray, (runlen, 2))
			myarray[:,0] = myarray[:,0] + sim*histseg.Global_BINWIDTH
			myarrays.append(myarray)
			runlens.append(runlen)
			self.logToMaster("runlens is %d, and runlen is %d\n" % (len(runlens), runlen))
		runlen=max(runlens)
		mystats=np.array([], dtype=int)
		for s in xrange(len(mysims)):
			sim=mysims[s]  
			histstats=myarrays[s]
			i=sim*runlen
			if mystats.size==0: 
				mystats=np.zeros(((max(mysims)+1)*runlen, histstats.shape[1]), dtype=int)
			mystats[i:(i+histstats.shape[0]),:]=histstats
		np.savetxt(os.path.join(self.cnavgpost, "historylengths.txt"), mystats, fmt='%d', delimiter="\t")
	

def main():
	#parser = OptionParser(usage = "Right now the do_something is to count the number of events per history from the .braney files and print it to a history_eventcounts.txt file.") 
	parser = OptionParser(usage = "Right now the do_something is to run get_event_order_counts on the list of samples in samplelist.") 
	parser.add_option("--samplelist", dest="samplelist", help="The list of CNAVG outputs and sample ids. Should have the form <directory><ID>.")
	parser.add_option("--edges", dest="edges", action='store_true', help="Do the analysis using the .pedges instead of .pevnts files.")
	parser.add_option('--simulation', dest="simulation", action='store_true', help='whether the history is a simulation')
	parser.add_option('--cutoff', dest="cutoff", help='only look at events with a likelihood above this cutoff', default=0, type=float)
	Stack.addJobTreeOptions(parser)
	options, args = parser.parse_args()
	i = Stack(Setup(options)).startJobTree(options)
	if i:
		raise RuntimeError("The jobtree contains %d failed jobs.\n" % i)


if __name__=="__main__":
	from analyze_event_orders_jobtree import *
	main()

