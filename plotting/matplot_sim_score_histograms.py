#!/inside/home/common/bin/python2.7 

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import argparse 
import os, sys
#from IPython.parallel.util import interactive 
#from ../event_cycles_module import Global_EVENTTYPES
Global_types=['any', 'amplification', 'deletion', 'adjacency']
class SimData: 
	def __init__(self, statsfile): 
		f=open(statsfile, 'r') 
		lines=f.readlines()
		header=lines[0]
		self.TP=map(int, (lines[1].strip().split('\t')[1:5]))
		self.FP=map(int, (lines[2].strip().split('\t')[1:5]))
		self.FN=map(int, (lines[3].strip().split('\t')[1:5]))
		self.TN=map(int, (lines[4].strip().split('\t')[1:5]))
		self.F1score=float(lines[5].strip().split('\t')[1])

def add_scores(datfile, truescores, fpscores):
	datfh=open(datfile, 'r')
	datfh.readline() # skip the header
	for line in datfh: 
		data=line.strip().split('\t')
		true=int(data[3])
		lscore=float(data[1])
		if lscore >1: 
			sys.stderr.write("Bad lscore: %f, %s\n%s" % (lscore, datfile, line))
		if true==1:
			truescores.append(lscore)
		elif true==0:
			fpscores.append(lscore)
			
def make_histograms(truescores, fpscores, title):
	plt.hist(truescores, normed=True, histtype='step', linewidth=3, label="True Events, %d" % len(truescores))	
	plt.hist(fpscores, normed=True, histtype='step', linewidth=3, label="FP Events, %d" % len(fpscores))
	plt.legend()
	plt.xlabel("Likelihood Score")
	plt.ylabel("Density of events")
	plt.title("Histogram of TP and FP Event Likelihoods")	

#@interactive 
def main(simulations, outname, edgeorevent, plot, savedata):
	mysimulations=open(simulations, 'r').readlines()
	truescores=[]
	fpscores=[]
	type=0
#	(min_count, max_count) = (5,12)
	(min_count, max_count) = (0,500)
	for i in xrange(len(mysimulations)): 
		line=mysimulations[i]
		if len(line.strip().split('\t')) ==3: 
			(dir, blocks, events)=line.strip().split('\t')
		else: 
			(dir, simid, blocks, events)=line.strip().split('\t')
		statsfile=os.path.join(dir, "%s.stats" % edgeorevent)
		data=SimData(statsfile)
		truecount=data.TP[type] + data.FN[type]
		if truecount >= min_count and truecount <= max_count:
			sys.stderr.write("truecount is %d\n" % (truecount))
			datfile=os.path.join(dir, "%s.dat" % edgeorevent)
			add_scores(datfile, truescores, fpscores)
	if savedata: 
		np.savetxt("%s.tp.txt" % outname, np.array(truescores))
		np.savetxt("%s.fp.txt" % outname, np.array(fpscores))
	if plot: 
		title=edgeorevent
		sys.stderr.write("making plot...\n")
		make_histograms(truescores, fpscores, title) 	
		plt.savefig(outname+".png")	

if __name__ == "__main__": 
	parser=argparse.ArgumentParser(description="makes plots")
	parser.add_argument('simulations', help='A list of simulation directories')
	parser.add_argument('outname', help='The basename for the files that are made')
	parser.add_argument('--plot', help='whether to plot a figure or not.', action='store_true')
	parser.add_argument('--savedata', help='whether to print that data out or not.', action='store_true')
	parser.add_argument('--edge', help='whether to analyze data for edges rather than events', action='store_true')
	args=parser.parse_args()
	if args.edge: 
		edgeorevent="edges"
	else: 
		edgeorevent="events"
	main(args.simulations, args.outname, edgeorevent, args.plot, args.savedata)


