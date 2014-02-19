#!/inside/home/common/bin/python2.7

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from optparse import OptionGroup
from optparse import OptionParser
from sonLib.bioio import logger
from sonLib.bioio import system

import subprocess

from score_edges_within_pevents import *
from score_and_link_cycles import *
import event_cycles_module as histseg
import analyze_simulation
import annotate_events
import argparse
import mcmc_mixing_analysis_jobtree as mcmcjobtree

#======= Main Setup ===================
class Setup(Target):
	def __init__(self, options): 
		Target.__init__(self)
		self.options=options
		self.events=[]
		self.totalp=1
		self.edges=[]	
	
	def run(self):
		self.logToMaster("Setting up...") 
		opts=self.options
		sampleid=opts.sampleid
		outputdir=opts.outputdir
		subprocess.call("mkdir -p %s" % outputdir, shell=True)
		pevntsfile=os.path.join(outputdir, opts.sampleid + ".pevnts")
		logger.info("pevntsfile: %s" % pevntsfile)
		if opts.pevnts or not os.path.exists(pevntsfile): 
			logger.info("pevntsfile: %s" % pevntsfile)
			CreatePevntsFile(pevntsfile, opts).run()
#		self.events=pickle.load(open(pevntsfile, 'rb'))
#		self.totalp=histseg.get_total_lieklihood(events)
		
		linksfile =os.path.join(outputdir, sampleid +".links")
		if opts.links and not os.path.exists(linksfile): 
			logger.info("linksfile: %s" % linksfile)
			self.addChildTarget(CreateLinksFile(pevntsfile, linksfile))		
		
		annotationfile=os.path.join(outputdir, sampleid + ".ann")
		if opts.ann and not os.path.exists(annotationfile): 
			logger.info("annotationfile: %s" % annotationfile)
			if not self.events: 
				self.events=pickle.load(open(pevntsfile, 'rb'))
			self.addChildTarget(CreateAnnotationFile(self.events, opts.tabixfile, annotationfile))
		
		pedgesfile=os.path.join(outputdir, sampleid + ".pedgs")
		if opts.pedges or not os.path.exists(pedgesfile):
			logger.info("pedgesfile: %s" % pedgesfile)
			CreatePedgesFile(pevntsfile, pedgesfile).run()
		
		mcmcdir=os.path.join(outputdir, "mcmcdata")
		if opts.mcmcmix and not os.path.exists(mcmcdir):
			opts.pevnts=pevntsfile 
			opts.pedges=pedgesfile
			opts.outputdir=mcmcdir
			self.addChildTarget(mcmcjobtree.SetupMCMC(opts))
		simoutput=os.path.join(outputdir, "edges.stats")
		if opts.simulation and ((not os.path.exists(simoutput)) or (os.path.getsize(simoutput) == 0)): 
			if not self.events: 
				logger.info("loading events now %s" % pevntsfile)
				self.events=pickle.load(open(pevntsfile, 'rb'))
			if not self.edges: 
				self.edges=pickle.load(open(pedgesfile, 'rb'))
			self.addChildTarget(SimAnalysisJob(self.events,"events", opts, pevntsfile))
			self.addChildTarget(SimAnalysisJob(self.edges, "edges", opts, pedgesfile))

class SimAnalysisJob(Target): 
	def __init__(self, events, outname, options, file): 
		Target.__init__(self)
		self.options=options
		self.events=events
		self.outname=outname
		self.pevntsfile=file

	def run(self): 
		opts=self.options
		datout=os.path.join(opts.outputdir, "%s.dat" % self.outname)
		statout=os.path.join(opts.outputdir, "%s.stats" % self.outname)
		parser=argparse.ArgumentParser()
		analyze_simulation.add_options(parser)
		myargs=['--datout', datout, '--stats', statout, self.pevntsfile, '0']
		args=parser.parse_args(myargs)
		analyze_simulation.analyze_simulation(self.events, args)

class CreatePevntsFile(Target): 
	def __init__(self, pevntsfile, options): 
		Target.__init__(self)
		self.options=options
		self.pevntsfile=pevntsfile
	
	def run(self):
		self.logToMaster("CreatePevntsFile\n") 
		opts=self.options
		if opts.simulation:
			truehist=os.path.join(opts.cnavgout, "HISTORIES_0.braney") 
			if not os.path.exists(truehist):
				truefile=os.path.join(opts.cnavgout, "true.braney")
				subprocess.call("grep -v ^$ %s | gzip > %s" % (truefile, truehist), shell=True)
		myargs=['--cnavg', opts.cnavgout, '--outpickle', self.pevntsfile]
		if opts.links: 
			linksfile=os.path.join(opts.outputdir, opts.sampleid + ".links") 
			myargs = myargs + ['--links', linksfile]
		parser = argparse.ArgumentParser()
		add_event_link_options(parser)
		args=parser.parse_args(myargs)
		score_and_link_cycles(args)

class CreateLinksFile(Target): 
	def __init__(self, pevntsfile, linksfile): 
		Target.__init__(self)
		self.pevntsfile=pevntsfile
		self.linksfile=linksfile
	
	def run(self): 
		self.logToMaster("CreateLinksFile\n") 
		opts=self.options
		myargs=['--inpickle', self.pevntsfile, '--links', self.linksfile]
		parser = argparse.ArgumentParser()
		add_event_link_options(parser)
		args=parser.parse_args(myargs)
		score_and_link_cycles(args)

class CreatePedgesFile(Target): 
	def __init__(self, pevntsfile, pedgesfile): 
		Target.__init__(self)
		self.pevntsfile=pevntsfile
		self.pedgesfile=pedgesfile
		
	def run(self):
		self.logToMaster("CreatePedgesFile\n") 
		myargs=['--outpickle', self.pedgesfile, self.pevntsfile]
		parser=argparse.ArgumentParser()
		add_score_edges_options(parser)
		args=parser.parse_args(myargs)
		score_edges_within_pevents(args)
			
class CreateAnnotationFile(Target):
	def __init__(self, evnts, tabixfile, annotationfile): 
		Target.__init__(self)
		self.evnts=evnts
		self.tabixfile=tabixfile
		self.annfile=annotationfile
			
	def run(self): 
		self.logToMaster("CreateAnnotationFile\n") 
		outputfh=open(self.annfile, 'w')
		allevents=self.evnts
		mytabix=pysam.Tabixfile(self.tabixfile, 'r')
		annotate_events.main(allevents, mytabix, outputfh)	

def add_analysis_options(parser): 
	group=OptionGroup(parser, "CNAVG Post Analysis Options")
	group.add_option('--pevnts', dest="pevnts", default=False, action="store_true", help="create or rewrite .pevnts file")
	group.add_option('--pedges', dest="pedges", default=False, action="store_true", help="create or rewrite .pedges file")
	group.add_option('--links', dest="links", default=False, action="store_true", help="create or rewrite .links file")
	group.add_option('--ann', dest="ann", default=False, action="store_true", help="create or rewrite .annotation file")
	group.add_option('--tabixfile', dest="tabixfile", help="The tabix file containing gene annotations you are interested in.")
	group.add_option('--mcmcmix', dest="mcmcmix", default=False, action="store_true", help="do mcmc analysis.")
	group.add_option('--simulation', dest="simulation", default=False, action="store_true", help="do simulation analysis.")
	group.add_option("--cnavgout", dest="cnavgout", help="The CN-AVG output directory for the sample")
	group.add_option("--outputdir", dest="outputdir", help="The output directory for the analysis")
	group.add_option("--sampleid", dest="sampleid", help="the sampleid")
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
