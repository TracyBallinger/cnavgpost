#!/inside/home/common/bin/python2.7

''' 
Analysis on CN-AVG output to test how well histories are being sampled across runs 
''' 
import os, sys
import re
import subprocess, pickle, glob
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack 
from optparse import OptionGroup
from optparse import OptionParser
from sonLib.bioio import logger
import numpy as np
import get_history_distances_between_mcmc_steps as mcmcdist 
import event_cycles_module as histseg 


#======== MAIN COMMAND ============
class SetupMCMC(Target): 
	def __init__(self, options): 
		Target.__init__(self)
		self.options=options
		
	def run(self): 
		args=self.options
		#global_dir=self.getGlobalTempDir()
		global_dir=args.outputdir
		cmd = "mkdir -p %s" % global_dir 
		logger.info(cmd)
		subprocess.call(cmd, shell=True)
		logger.info("loading events: %s\n" % args.pevnts)
		events=pickle.load(open(args.pevnts, 'rb'))
		evntdatadir=os.path.join(global_dir, "eventdata")
		subprocess.call("mkdir -p %s" % evntdatadir, shell=True)
		logger.info("loading edges: %s\n" % args.pedges)
		edges=pickle.load(open(args.pedges, 'rb'))
		edgedatadir=os.path.join(global_dir, "edgedata")
		subprocess.call("mkdir -p %s" % edgedatadir, shell=True)
		logger.info("numruns %d" % (args.numruns))
		for i in xrange(args.numruns):
			logger.info("Adding Child %d" % i)
			self.addChildTarget(DoMCMCforSingleRun(i, args, events, evntdatadir))
			self.addChildTarget(DoMCMCforIndepRuns(i, args, events, evntdatadir))
			self.addChildTarget(DoMCMCforSingleRun(i, args, edges, edgedatadir))
			self.addChildTarget(DoMCMCforIndepRuns(i, args, edges, edgedatadir))
		logger.info("error here %s\n" % ("after children"))
		self.setFollowOnTarget(CombineMCMCdataRuns(global_dir, evntdatadir, edgedatadir))

class DoMCMCforIndepRuns(Target): 
	def __init__(self, i, options, events, datadir): 
		Target.__init__(self)
		self.options=options
		self.i=i
		self.outdir=datadir
		self.events=events

	def run(self): 
		args=self.options
		refhistoryid=args.refhistoryid+ args.binwidth*self.i 
		i2 = self.i +1
		if self.i== args.numruns-1:
			i2=0 
		refhistoryid2=args.refhistoryid + args.binwidth * i2 
		outfile=os.path.join(self.outdir, "runs_%d.%d.dat" % (self.i, i2))
		outfh=open(outfile, 'w')
		logger.info("running %d: get_history_distances_between_mcmc_steps(events, %s, %s, %s, %s) > %s" % (self.i, args.refhistoryid, refhistoryid2, args.numsteps, args.stepsize, outfile))
		mcmcdist.get_history_distances_between_mcmc_steps(self.events, refhistoryid, refhistoryid2, args.numsteps, args.stepsize, outfh) 
		logger.info("error here %s\n" % ("hello"))

class DoMCMCforSingleRun(Target): 
	def __init__(self, i, options, events, datadir): 
		Target.__init__(self)
		self.options=options
		self.i=i
		self.outdir=datadir
		self.events=events

	def run(self): 
		args=self.options
		refhistoryid=args.refhistoryid + args.binwidth*self.i 
		outfile=os.path.join(self.outdir, "run_%d.dat" % self.i)
		outfh=open(outfile, 'w')
		logger.info("running %d: get_history_distances_between_mcmc_steps(events, %s, %s, %s, %s) > %s" % (self.i, refhistoryid, "" , args.numsteps, args.stepsize, outfile))
		mcmcdist.get_history_distances_between_mcmc_steps(self.events, refhistoryid, "", args.numsteps, args.stepsize, outfh) 
		logger.info("error here %s\n" % ("hello"))


class CombineMCMCdataRuns(Target): 
	def __init__(self, global_dir, evntdatadir, edgedatadir):
		Target.__init__(self)
		self.global_dir=global_dir
		self.evntdatadir=evntdatadir
		self.edgedatadir=edgedatadir

	def run(self):
		global_dir=self.global_dir 
		datafiles=glob.glob(self.evntdatadir+"/"+"run_*.dat")
		outfile=os.path.join(global_dir, "evntmixing_dep.dat")
		combine_datafiles(datafiles, outfile)
		datafiles=glob.glob(self.evntdatadir+"/"+"runs_*.*.dat")
		outfile=os.path.join(global_dir, "evntmixing_ind.dat")
		combine_datafiles(datafiles, outfile)
		datafiles=glob.glob(self.edgedatadir+"/"+"run_*.dat")
		outfile=os.path.join(global_dir, "edgemixing_dep.dat")
		combine_datafiles(datafiles, outfile)
		datafiles=glob.glob(self.edgedatadir+"/"+"runs_*.*.dat")
		outfile=os.path.join(global_dir, "edgemixing_ind.dat")
		combine_datafiles(datafiles, outfile)
	
	
def combine_datafiles(datafiles, outfile): 
	mydata=np.array([], dtype=float)
	myfmts=[]
	myheader=[]
	for i in xrange(len(datafiles)): 
		datafile=datafiles[i] 
		label=re.search('(run.*)\.dat', datafile).group(1)
		logger.info("datafile is " + datafile)
		data=np.loadtxt(datafile, skiprows=1)
		logger.info("data is %s" % str(data.shape))
		if i ==0:
			mydata=np.int_(data[:,0:2])
			myfmts.append("%d")	
			myfmts.append("%d")
			myheader.append("historyID_1")	
			myheader.append("historyID_2")	
			#p = subprocess.Popen(["sed 1d %s | awk '{print $5/($3+$4-$5)}'" % datafile], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			#out, err = p.communicate()
			#mydata.append(np.array(out.split('\n')), axis=0)
		x=data[:,4]/(data[:,2]+data[:,3]-data[:,4])
		x=x.reshape(data.shape[0], 1)
		myfmts.append("%f")
		myheader.append(label)	
		logger.info("x is %s, mydata is %s" % (str(x.shape), str(mydata.shape)))
		mydata=np.append(mydata, x, axis=1)
	np.savetxt(outfile, mydata, fmt="\t".join(myfmts), header="\t".join(myheader), delimiter="\t")	

def add_mcmc_options(parser): 
	group = OptionGroup(parser, "MCMC analysis options")
	group.add_option('--refhistoryid', dest='refhistoryid', help='the maximum number of steps in the sampling chain', default=2500, type="int")
	group.add_option('--binwidth', dest='binwidth', help='the multiplier between history ids of independent runs', default=10000, type="int")
	group.add_option('--numsteps', dest="numsteps", help='the number of steps to take.', default=0, type="int")
	group.add_option('--stepsize', dest="stepsize", help='the step size from the reference id(s) in the mcmc that you want to go.', default=-1, type="int")
	group.add_option('--numruns', dest="numruns", help='the number of independent histories.' , default=10, type="int")
	parser.add_option_group(group)

def main(): 
	parser = OptionParser(usage = "mcmc_mixing_analysis_jobtree.py --pevnts sample.pevnts --refhistoryid=2500 --numsteps=1000 --stepsize=1 --logInfo --jobTree=/outputdir --batchSystem=singleMachine")
	parser.add_option("--pevnts", dest="pevnts", help="a .pevnts file", type="string")
	parser.add_option("--pedges", dest="pedges", help="a .pedges file", type="string")
	parser.add_option("--outputdir", dest="outputdir", help="where you want the final output to go", type="string")
	add_mcmc_options(parser)
	Stack.addJobTreeOptions(parser)
	options, args = parser.parse_args()
	i = Stack(SetupMCMC(options)).startJobTree(options)
	if i: 
		raise RuntimeError("The jobtree contains %d failed jobs.\n" % i)

if __name__ == "__main__": 
	from mcmc_mixing_analysis_jobtree import *
	main()
	
