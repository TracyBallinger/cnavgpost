#!/inside/home/common/bin/python2.7 

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from optparse import OptionGroup
from optparse import OptionParser
from sonLib.bioio import logger

import subprocess, os 

import cnavg.simulator.simulator
import cnavg.cactus.oriented
import cnavg.cactus.graph
import cPickle as pickle


class SetupSim(Target):
	def __init__(self, options):
		Target.__init__(self)
		self.options=options
	
	def run(self):
		opts=self.options
		specfile=opts.specfile
		outputdir=opts.outputdir
		for specs in open(specfile, 'r'): 
			(blocks, events) = map(int, specs.strip().split())
			for i in xrange(opts.reps): 
				self.addChildTarget(RunCnavgSimulation(blocks, events, i, outputdir, opts))			
class Chdir: 
	def __init__(self, newPath):
		self.savedPath=os.getcwd()
		os.chdir(newPath)
	
	def __del__(self):
		os.chdir(self.savedPath)

class RunCnavgSimulation(Target): 
	def __init__(self, blocks, events, i, outputdir, options): 
		Target.__init__(self)
		self.options=options
		self.blocks=blocks
		self.events=events
		self.outputdir=outputdir
		self.i = i
	
	def run(self):
		opts=self.options 
		simoutdir=os.path.join(self.outputdir, "sim%d.%d.%d.%d" % (self.blocks, self.events, 1, self.i))
		subprocess.call("mkdir -p %s" % simoutdir, shell=True)
		cwd=os.getcwd()
		os.chdir(simoutdir)
		if opts.integer: 
			H = cnavg.simulator.simulator.RandomIntegerHistory(self.blocks, self.events, indel_length=10)
		else: 
			H = cnavg.simulator.simulator.RandomCancerHistory(self.blocks, self.events, indel_length=10)
		A = H.avg()
		C = cnavg.cactus.graph.Cactus(A)
		G = cnavg.cactus.oriented.OrientedCactus(C)

		file = open('true.braney', 'w')
		file.write(H.braneyText(G))
		file.close()

		file = open('CACTUS', 'w')
		pickle.dump(G, file)
		file.close()
		os.chdir(cwd)
		
		for j in xrange(1,10): 
			self.addChildTarget(RunCnavgForSim(simoutdir, j, opts))

class RunCnavgForSim(Target):
	def __init__(self, simoutdir, j, options): 
		Target.__init__(self)
		self.simoutdir=simoutdir
		self.j=j
		self.options=options

	def run(self):
		cwd=os.getcwd()
		os.chdir(self.simoutdir)
		if self.options.integer: 
			subprocess.call("/inside/home/dzerbino/cn-avg/bin/cn-avg.py -n -d . -i %d -s 1000" % self.j, shell=True)
		else: 
			subprocess.call("/inside/home/dzerbino/cn-avg/bin/cn-avg.py -d . -i %d -s 1000" % self.j, shell=True)
		os.chdir(cwd)

def add_simulation_options(parser): 
	group = OptionGroup(parser, "Make Simulation Options")
	group.add_option("--specfile", dest='specfile', help='A file specifying the number of blocks and events to put in each simulation', type='string')
	group.add_option("--reps", dest='reps', type="int", help="The number of reps to do for each simulation")
	group.add_option("--outputdir", dest='outputdir', help="The directory to write to")
	group.add_option("--integer", dest='integer', default=False, action="store_true", help="Use integer simulations rather than metagenomic.")
	parser.add_option_group(group)	

def main(): 
	parser= OptionParser(usage = "make_simulations_jobtree.py ....")
	add_simulation_options(parser)
	Stack.addJobTreeOptions(parser)
	options, args = parser.parse_args()
	i = Stack(SetupSim(options)).startJobTree(options)
	if i: 
		raise RuntimeError("The jobtree contains %d failed jobs.\n" % i)


if __name__== "__main__": 
	from make_simulations_jobtree import * 
	main()

 
