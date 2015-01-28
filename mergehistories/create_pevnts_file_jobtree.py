#!/usr/bin/env python 

''' 
Analysis on CN-AVG output to test how well histories are being sampled across runs 
''' 
import os, sys
import re
import subprocess, glob
import cPickle as pickle
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack 
from optparse import OptionGroup
from optparse import OptionParser
from sonLib.bioio import logger
import numpy as np
import cnavgpost.mergehistories.event_cycles_module as histseg 


#======== MAIN COMMAND ============
class SetupCreatePevnts(Target): 
	def __init__(self, options, outputdir): 
		Target.__init__(self)
		self.options=options
		self.outputdir=outputdir
		
	def run(self): 
		self.logToMaster("Merging .braney files...")
		opts=self.options
		outputdir=self.outputdir
		#global_dir=self.getGlobalTempDir()
		global_dir=self.outputdir
		cmd = "mkdir -p %s" % global_dir 
		logger.info(cmd)
		subprocess.call(cmd, shell=True)
		# get historystats
		historyScores=np.loadtxt(opts.historystatsfile, dtype=int)
		totalp=histseg.compute_likelihood_histories(historyScores[:,0], historyScores)
		pevntsfile=os.path.join(outputdir, opts.sampleid + ".pevnts")
		self.addChildTarget(CreatePevntsFile(pevntsfile, historyScores, totalp, opts))


class CreatePevntsFile(Target):
	def __init__(self, pevntsfile, historyScores, totalp, options):
		Target.__init__(self)
		self.pevntsfile=pevntsfile
		self.historyScores=historyScores
		self.totalp=totalp
		self.opts=options

	def run(self):
		self.logToMaster("CreatePevntsFile\n")
		opts=self.opts
		braneyfiles=glob.glob(opts.cnavgout+"/"+"HISTORIES_*.braney")
		sys.stderr.write("braneyfiles: %s\n" % (str(braneyfiles)))
		pevntsfiles=[]
		for braneyfn in braneyfiles:
			sim=int(re.match(".*HISTORIES_(\d+)\.braney", braneyfn).group(1))
			evntsfile=os.path.join(opts.outputdir, "sim.%d.pevnts" % sim)
			pevntsfiles.append(evntsfile)
			# This will create a pevnts file that is sorted and untrimmed. 
	#		self.addChildTarget(MergeSingleBraneyFile(braneyfn, evntsfile, sim))
			MergeSingleBraneyFile(braneyfn, evntsfile, sim).run()
			sys.stderr.write("finished Creating file %s\n" % (evntsfile))
	#	self.setFollowOnTarget(MergePevntsFiles(pevntsfiles, self.pevntsfile, self.historyScores, self.totalp))
		MergePevntsFiles(pevntsfiles, self.pevntsfile, self.historyScores, self.totalp).run()

class MergePevntsFiles(Target): 
	def __init__(self, pevntsfiles, pevntsout, historyScores, totalp): 
		Target.__init__(self)
		self.pevntsfiles=pevntsfiles
		self.pevntsout=pevntsout
		self.historyScores=historyScores
		self.totalp=totalp

	def run(self): 
		histseg.merge_pevnts_files(self.pevntsfiles, self.pevntsout, self.historyScores, self.totalp)
		# clean up 
		for f in self.pevntsfiles: 
			os.remove(f)

class MergeSingleBraneyFile(Target):
	def __init__(self, braneyfn, evntsfile, sim):
		Target.__init__(self)
		self.braneyfn=braneyfn
		self.evntsfile=evntsfile
		self.sim=sim

	def run(self):
		sim=self.sim
		sys.stderr.write("Creating file %s\n" % (self.evntsfile))
		events=histseg.make_events_from_braneyfn(self.braneyfn)
		# Need to change the id to distinguish from other braney files. 
		for evnt in events:
			for i in xrange(len(evnt.histories)):
				evnt.histories[i] = evnt.histories[i] + sim*histseg.Global_BINWIDTH
				(idhist, order) = map(int, evnt.id.split('.'))
			evnt.id = "%d.%d" % (idhist+(sim*histseg.Global_BINWIDTH), order)
		sortedevents=sorted(events, key=lambda x:(x.segstr, x.cnval, x.prevals[0]))	
		pickle.dump(sortedevents, open(self.evntsfile, 'wb'), pickle.HIGHEST_PROTOCOL)
		self.logToMaster("Created file %s\n" % (self.evntsfile))

def main(): 
	parser = OptionParser(usage = "create_pevnts_file_jobtree.py --cnavgout cnavgdir --outputdir outputdir ... jobtree_options\n")
	parser.add_option("--outputdir", dest="outputdir", help="where you want the created .pvents files to go", type="string")
	parser.add_option("--cnavgout", dest="cnavgout", help="where the .braney files are.", type="string")
	parser.add_option("--sampleid", dest="sampleid", help="The prefix of the .pevnts file.  ie sampleid.pevnts will be created.", type="string")
	parser.add_option("--histstats", dest="historystatsfile", help="The historystats.txt file.", type="string")
	parser.add_option("--binwidth", dest="binwidth", help="the multiplier for each history id to distinguish independent simulations.", type=int, default=histseg.Global_BINWIDTH)
	Stack.addJobTreeOptions(parser)
	options, args = parser.parse_args()
	histseg.Global_BINWIDTH=options.binwidth
	i = Stack(SetupCreatePevnts(options, options.outputdir)).startJobTree(options)
	if i: 
		raise RuntimeError("The jobtree contains %d failed jobs.\n" % i)


if __name__ == "__main__": 
	from create_pevnts_file_jobtree import *
	main()
	
