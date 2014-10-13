#!/usr/bin/env python

# This uses jobTree (https://github.com/benedictpaten/jobTree) to run analysis on the CNAVG output.  

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from optparse import OptionGroup
from optparse import OptionParser
from sonLib.bioio import logger
from sonLib.bioio import system

import subprocess, os
import argparse
import numpy as np
import pysam 
import re, glob
import cPickle as pickle 

from cnavgpost.mergehistories.score_edges_within_pevents import *
from cnavgpost.mergehistories.score_and_link_cycles import *
import cnavgpost.mergehistories.event_cycles_module as histseg
import cnavgpost.mergehistories.create_pevnts_file_jobtree as pevntsjobtree 
import cnavgpost.simulations.analyze_simulation as analyze_simulation
import cnavgpost.genehistory.annotate_events as annotate_events
import cnavgpost.genehistory.histories_to_gene_orders as histories_to_gene_orders
import cnavgpost.diagnostics.mcmc_mixing_analysis_jobtree as mcmcjobtree

#======= Main Setup ===================
class Setup(Target):
	def __init__(self, options): 
		Target.__init__(self)
		self.options=options
		self.events=[]
		self.totalp=0
		self.edges=[]
		self.historyScores=[]	
	
	def run(self):
		self.logToMaster("Setting up...") 
		opts=self.options
		histseg.Global_BINWIDTH=opts.binwidth
		sampleid=opts.sampleid
		outputdir=opts.outputdir
		subprocess.call("mkdir -p %s" % outputdir, shell=True)

		historystatsfile=os.path.join(outputdir, "historystats.txt")
		if not os.path.exists(historystatsfile): 
			logger.info("historystatsfile: %s" % historystatsfile)
			CombineHistoryStatsfiles(opts, historystatsfile).run()
		
		self.historyScores=np.loadtxt(historystatsfile, dtype=int)
		self.totalp=histseg.compute_likelihood_histories(self.historyScores[:,0], self.historyScores)
		logger.info("Global_BINWIDTH: %d" % histseg.Global_BINWIDTH)
		logger.info("totalp is %s" % str(self.totalp))	

		pevntsfile=os.path.join(outputdir, opts.sampleid + ".pevnts")
		if opts.pevnts or not os.path.exists(pevntsfile): 
			logger.info("pevntsfile: %s" % pevntsfile)
			pevntsjobtree.CreatePevntsFile(pevntsfile, self.historyScores, self.totalp, opts).run()

		pedgesfile=os.path.join(outputdir, sampleid + ".pedgs")
		if opts.pedges or not os.path.exists(pedgesfile):
			logger.info("pedgesfile: %s" % pedgesfile)
			CreatePedgesFile(pickle.load(open(pevntsfile, 'rb')), pedgesfile, self.historyScores, self.totalp, ignore_cn=False).run()
		
		mrgpedgesfile=os.path.join(outputdir, sampleid + ".mrgpedgs")
		if not os.path.exists(mrgpedgesfile):
			logger.info("mrgpedgesfile: %s" % mrgpedgesfile)
			CreatePedgesFile(pickle.load(open(pevntsfile, 'rb')), mrgpedgesfile, self.historyScores, self.totalp, ignore_cn=True).run()
		
		linksfile =os.path.join(outputdir, sampleid +".links")
		if opts.links and not os.path.exists(linksfile): 
			logger.info("linksfile: %s" % linksfile)
			self.addChildTarget(CreateLinksFile(pevntsfile, linksfile, self.totalp))		
	
		breaksfile=os.path.join(outputdir, "breakpoints.txt")
		if not os.path.exists(breaksfile): 
			breaklocs=histseg.get_breakpoints(pickle.load(open(pedgesfile, 'rb')), opts.trueID)
			breaklocs2=histseg.get_breakpoints(pickle.load(open(mrgpedgesfile, 'rb')), opts.trueID)
			breaksfh=open(breaksfile, 'w')
			for loc in sorted(breaklocs.keys()):
				(n, t) = breaklocs[loc]
				(n2, t2) = breaklocs2[loc]
				breaksfh.write("%s\t%d\t%d\t%d\t%d\n" % (loc, n, t, n2, t2))	
		annotationfile=os.path.join(outputdir, sampleid + ".ann")
		generankfile=os.path.join(outputdir, sampleid + ".gnrank")
		if opts.generank and not os.path.exists(generankfile): 
			logger.info("generankfile: %s" % generankfile)
			if not self.events: 
				self.events=pickle.load(open(pevntsfile, 'rb'))
			self.addChildTarget(CreateGeneRankFile(self.events, opts.tabixfile, self.totalp, annotationfile, generankfile))	
			logger.info("after creating generankfile")
		elif opts.ann and not opts.generank: 
			logger.info("annotationfile: %s" % annotationfile)
			if not self.events: 
				self.events=pickle.load(open(pevntsfile, 'rb'))
			self.addChildTarget(CreateAnnotationFile(self.events, opts.tabixfile, annotationfile))
		
		if opts.mcmcmix: 	
			mcmcdir=os.path.join(outputdir, "mcmcdata")
			mcmcdat=os.path.join(mcmcdir, "edge_counts.dat")
			mcmcdir=os.path.join(outputdir, "mcmcdata")
			mcmcdat=os.path.join(mcmcdir, "edge_counts.dat")
			if not os.path.exists(mcmcdir) or not os.path.exists(mcmcdat):
				subprocess.call("mkdir -p %s" % mcmcdir, shell=True)	
				opts.pevnts=pevntsfile 
				opts.pedges=pedgesfile
				self.addChildTarget(mcmcjobtree.SetupMCMC(opts, mcmcdir))

		if opts.simulation: 
			simoutput=os.path.join(outputdir, "events.stats")
			if ((not os.path.exists(simoutput)) or (os.path.getsize(simoutput) == 0)): 
				self.addChildTarget(SimAnalysisJob(pevntsfile, opts.trueID, self.historyScores, "events", outputdir, opts.binwidth))
			simoutput2=os.path.join(outputdir, "edges.stats")
			if ((not os.path.exists(simoutput2)) or (os.path.getsize(simoutput2) == 0)): 
				self.addChildTarget(SimAnalysisJob(pedgesfile, opts.trueID, self.historyScores, "edges", outputdir, opts.binwidth))
			simoutput3=os.path.join(outputdir, "mrgedges.stats")
			if ((not os.path.exists(simoutput3)) or (os.path.getsize(simoutput3) == 0)): 
				self.addChildTarget(SimAnalysisJob(mrgpedgesfile, opts.trueID, self.historyScores, "mrgedges", outputdir, opts.binwidth))


class CreateGeneRankFile(Target):
	def __init__(self, events, tabixfile, totalp, annotationfile, generankfile):
		Target.__init__(self)
		self.events=events
		self.tabixfile=tabixfile
		self.annotationfile=annotationfile
		self.generankfile=generankfile
		self.totalp=totalp
	
	def run(self):
		logger.info("Should be creating %s" % self.annotationfile)
		CreateAnnotationFile(self.events, self.tabixfile, self.annotationfile).run()
		logger.info("Should be creating %s" % self.generankfile)
		histories_to_gene_orders.main(self.events, self.annotationfile, self.totalp, open(self.generankfile, 'w')) 	

#  This can be gotten rid of.... 
class CombineHistoryStatsfiles(Target):
	def __init__(self, options, historystatsfile):
		Target.__init__(self)
		self.options=options
		self.historystatsfile=historystatsfile
		self.cnavgout=options.cnavgout

	def run(self):
		opts=self.options
		if opts.simulation:
			truefile=os.path.join(opts.cnavgout, "true.braney")
			truehist=os.path.join(opts.cnavgout, "HISTORIES_0.braney") 
			subprocess.call("grep -v ^$ %s | gzip > %s" % (truefile, truehist), shell=True)
			make_STATS_from_truebraney(truefile, os.path.join(opts.cnavgout, "HISTORY_STATS_0"))
		statsfiles=glob.glob(self.cnavgout+"/"+"HISTORY_STATS*")
		sys.stderr.write("statsfiles: %s\n" % (str(statsfiles)))
		historyScores=histseg.combine_history_statsfiles(self.cnavgout)
		np.savetxt(self.historystatsfile, historyScores, fmt='%d', delimiter='\t')	


class SimAnalysisJob(Target): 
	def __init__(self, peventsfile, trueID, historyScores, outname, outputdir, binwidth): 
		Target.__init__(self)
		self.events=pickle.load(open(peventsfile, 'rb'))
		self.outname=outname
		self.trueID = trueID
		self.historyScores=historyScores
		self.outputdir=outputdir
		self.binwidth=binwidth
	
	def run(self): 
		datout=open(os.path.join(self.outputdir, "%s.dat" % self.outname), 'w')
		statout=open(os.path.join(self.outputdir, "%s.stats" % self.outname), 'w')
		#breaksfile=open(os.path.join(self.outputdir, "breakpoints.txt"), 'w')
		breaksfile=""
		histseg.Global_BINWIDTH=self.binwidth
		analyze_simulation.analyze_simulation(self.events, self.trueID, self.historyScores, datout, statout, breaksfile)

def make_STATS_from_truebraney(truefile, statsfn):
	cost=subprocess.check_output("grep ^A %s | head -1 | cut -f15,16" % (truefile), shell=True)
	costs=cost.strip().split()
	sys.stderr.write
	errorcost=0
	length=subprocess.check_output("grep ^A %s | cut -f14 | sort -u | wc -l" % (truefile), shell=True).strip()
	medianLen=0
	maxLen=0
	fout=open(statsfn, 'w')
	fout.write("%s\t%s\t%d\t%s\t%d\t%d\n" % (costs[0], costs[1], errorcost, length, medianLen, maxLen))
	fout.close()

class CreateLinksFile(Target): 
	def __init__(self, pevntsfile, linksfile, totalp): 
		Target.__init__(self)
		self.pevntsfile=pevntsfile
		self.linksfile=linksfile
		self.totalp=totalp
	
	def run(self): 
		self.logToMaster("CreateLinksFile\n") 
		opts=self.options
		myargs=['--inpickle', self.pevntsfile, '--links', self.linksfile, '--totalp', str(self.totalp)]
		parser = argparse.ArgumentParser()
		add_event_link_options(parser)
		args=parser.parse_args(myargs)
		score_and_link_cycles(args)

class CreatePedgesFile(Target): 
	def __init__(self, events, pedgesfile, historyScores, totalp, ignore_cn): 
		Target.__init__(self)
		self.events=events
		self.pedgesfile=pedgesfile
		self.historyScores=historyScores
		self.totalp=totalp
		self.ignore_cn=ignore_cn	
	
	def run(self):
		self.logToMaster("CreatePedgesFile\n") 
		edges=score_edges_within_pevents(self.events, self.historyScores, self.totalp, self.ignore_cn)
		pickle.dump(edges, open(self.pedgesfile, 'wb'), pickle.HIGHEST_PROTOCOL) 
	
class CreateAnnotationFile(Target):
	def __init__(self, evnts, tabixfile, annotationfile): 
		Target.__init__(self)
		self.evnts=evnts
		self.tabixfile=tabixfile
		self.annfile=annotationfile
			
	def run(self): 
		logger.info("maybe creating %s" % self.annfile)
		if not os.path.exists(self.annfile):
			logger.info("am creating %s" % self.annfile)
			self.logToMaster("CreateAnnotationFile\n") 
			outputfh=open(self.annfile, 'w')
			allevents=self.evnts
			mytabix=pysam.Tabixfile(self.tabixfile, 'r')
			annotate_events.main(allevents, mytabix, outputfh)	

def add_analysis_options(parser): 
	group=OptionGroup(parser, "CNAVG Post Analysis Options")
	group.add_option("--sampleid", dest="sampleid", help="the sampleid")
	group.add_option("--cnavgout", dest="cnavgout", help="The CN-AVG output directory for the sample")
	group.add_option("--outputdir", dest="outputdir", help="The output directory for the analysis")
	group.add_option('--binwidth', dest='binwidth', help='the multiplier between history ids of independent runs', default=histseg.Global_BINWIDTH, type="int")
	group.add_option('--pevnts', dest="pevnts", default=False, action="store_true", help="create or rewrite .pevnts file")
	group.add_option('--pedges', dest="pedges", default=False, action="store_true", help="create or rewrite .pedges file")
	group.add_option('--links', dest="links", default=False, action="store_true", help="create or rewrite .links file")
	group.add_option('--mcmcmix', dest="mcmcmix", default=False, action="store_true", help="do analysis to look at the mixing across and within runs.")
	group.add_option('--generank', dest="generank", default=False, action="store_true", help="create or rewrite .gnrnk file")
	group.add_option('--ann', dest="ann", default=False, action="store_true", help="create or rewrite .annotation file")
	group.add_option('--tabixfile', dest="tabixfile", help="The tabix file containing gene annotations you are interested in.")
	group.add_option('--simulation', dest="simulation", default=False, action="store_true", help="do simulation analysis.")
	group.add_option('--trueID', dest="trueID", default=0, help="The history id of the true simulated history.", type="int")
	parser.add_option_group(group)

def main(): 
	parser = OptionParser(usage = "cn-avg_post_analysis_jobtree.py --cnavgout CN-AVG_outputdir --sampleid Sampleid --jobTree jobtreedir")
	
	add_analysis_options(parser)
	mcmcjobtree.add_mcmc_options(parser)
	Stack.addJobTreeOptions(parser)
	options, args = parser.parse_args()
	i = Stack(Setup(options)).startJobTree(options)
	if i: 
		raise RuntimeError("The jobtree contains %d failed jobs.\n" % i)

if __name__ == "__main__": 
	from cnavg_post_analysis_jobtree import *
	main() 
