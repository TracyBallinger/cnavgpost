#!/usr/bin/env python 

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

#@interactive 
def main(simulations, outname, edgeorevent, eventtypes, plot, data, xlimit):
	mysimulations=open(simulations, 'r').readlines()
	datalist=[]
	for i in xrange(len(mysimulations)): 
		line=mysimulations[i]
		(dir, blocks, events)=line.strip().split('\t')
		statsfile=os.path.join(dir, "%s.stats" % edgeorevent)
		datalist.append(SimData(statsfile))
	ntypes=4
	myaccuracy
	truecounts=range(ntypes)
	f1scores=range(ntypes)
	precisions=range(ntypes)
	recalls=range(ntypes)
	accuracies=range(ntypes)
	adjaccuracies=range(ntypes)
	types=eventtypes
	for type in types: 
		truecount=np.zeros(len(mysimulations))
		total=np.zeros(len(mysimulations))
		f1score=np.zeros(len(mysimulations))
		precision=np.zeros(len(mysimulations))
		recall=np.zeros(len(mysimulations))
		accuracy=np.zeros(len(mysimulations))
		adjaccuracy=np.zeros(len(mysimulations))
		for i in xrange(len(mysimulations)):
			data=datalist[i]	
			truecount[i]=data.TP[type]+data.FN[type]
			total[i]= truecount[i]+data.TN[type]+data.FP[type]
			if data.TP[type] >0: 
				f1score[i]=float(2*data.TP[type])/float(2*data.TP[type]+data.FP[type] + data.FN[type])
				precision[i]=float(data.TP[type])/float(data.TP[type]+data.FP[type])
				recall[i]=float(data.TP[type])/float(data.TP[type]+data.FN[type])
				accuracy[i]=float(data.TP[type])/float(total[i])
				adjaccuracy[i]=float(data.TP[type]+data.TN[type])/float(total[i])
			else: 
				f1score[i]=-1
				precision[i]=-1
				recall[i]=-1
				accuracy[i]=-1
				adjaccuracy[i]=-1
				
		truecounts[type]=truecount
		f1scores[type]=f1score
		precisions[type]=precision
		recalls[type]=recall
		accuracies[type]=accuracy
		adjaccuracies[type]=adjaccuracy
		if data:  
			np.savetxt("%s.%d.txt" % (outname, type), np.vstack((truecount, total,recall, precision, accuracy, f1score )).T, fmt='%d\t%d\t%f\t%f\t%f\t%f', header="\t".join(("True_events", "Total_events", "recall", "precision", "accuracy", "f1score")))
	
	if plot:
		fig=plt.figure(figsize=(8,16))
		ax1=fig.add_subplot(4,1,1)
		ax1.set_ylabel("Recall")
		ax1.set_xlabel("Number of true events")
		ax2=fig.add_subplot(4,1,2, sharex=ax1, sharey=ax1)
		ax2.set_ylabel("Precision")
		ax2.set_xlabel("Number of true events")
		ax3=fig.add_subplot(4,1,3, sharex=ax1, sharey=ax1)
		ax3.set_ylabel("F1score")
		ax3.set_xlabel("Number of true events")
		ax4=fig.add_subplot(4,1,4, sharex=ax1, sharey=ax1)
		ax4.set_ylabel("Accuracy (adjusted)")
		ax4.set_xlabel("Number of true events")
		ax1.set_xlim([-1,xlimit])
		ax1.set_ylim([-.1,1.1])
	
		typcolors=('black', 'red','orange', 'turquoise')
		for type in types:
			truecount=truecounts[type]
			precision=precisions[type]
			recall=recalls[type]
			accuracy = adjaccuracies[type]
			f1score=f1scores[type]
			ax1.scatter(truecount, recall, color=typcolors[type], label=Global_types[type])
			ax2.scatter(truecount, precision, color=typcolors[type], label=Global_types[type])
			ax3.scatter(truecount, f1score, color=typcolors[type], label=Global_types[type])
			ax4.scatter(truecount, accuracy, color=typcolors[type], label=Global_types[type])
		if len(types)>1: 
			ax1.legend()	
			ax2.legend()	
			ax3.legend()	
			ax4.legend()	
		plt.savefig(outname+".png")

if __name__ == "__main__": 
	parser=argparse.ArgumentParser(description="makes plots")
	parser.add_argument('simulations', help='A list of simulation directories')
	parser.add_argument('outname', help='The basename for the files that are made')
	parser.add_argument('--plot', help='whether to plot a figure or not.', action='store_true')
	parser.add_argument('--data', help='whether to print that data out or not.', action='store_true')
	parser.add_argument('--xlimit', help='The max value on the x axis', default=200, type=int)
	parser.add_argument('--eventtypes', help='The type of events to include (amp, del, ins)')
	parser.add_argument('--edge', help='whether to analyze data for edges rather than events', action='store_true')
	args=parser.parse_args()
	if args.edge: 
		edgeorevent="edges"
	else: 
		edgeorevent="events"
	if args.eventtypes: 
		eventtypes=map(int, args.eventtypes.strip().split())
	else:
		 eventtypes=[0]
#main(simulations, outname, edgeorevent, eventtypes, plot, data, xlimit):
	
	main(args.simulations, args.outname, edgeorevent, eventtypes, args.plot, args.data, args.xlimit)


