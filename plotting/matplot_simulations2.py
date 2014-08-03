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
	def __init__(self, simdir, edgeorevent, pvalcutoff=0):
		sys.stderr.write("loading %s\n" % simdir)
		self.datfile=os.path.join(simdir, "%s.dat" % edgeorevent) 
		self.breaks=os.path.join(simdir, "breakpoints.txt")
		self.histstats=os.path.join(simdir, "historystats.txt")
		historyScores=np.loadtxt(self.histstats, dtype=int)
		self.truescore=np.mean(historyScores[0,1:3])
		self.avescore=np.mean(historyScores[1000:,1:3])
		self.minscore=np.min(np.mean(historyScores[1000:,1:3], axis=1))
		self.costdiff=self.avescore-self.truescore
#		sys.stderr.write("truescore: %f, simscore: %f, costdiff: %f\n" % (self.truescore, self.avescore, self.costdiff))
		nodes=np.loadtxt(self.breaks, usecols=(2,3))
		self.truenodes=sum(nodes[:,1]>0)
		self.totalnodes=nodes.shape[0]
		# get TP, etc. for calculating accuracies, etc. 
		eventtypes=['any', 'adj', 'amp', 'del']
		self.TP=range(len(eventtypes))
		self.FP=range(len(eventtypes))
		self.FN=range(len(eventtypes))
		self.TN=range(len(eventtypes))
		mydat=np.loadtxt(self.datfile, skiprows=1, usecols=(1,2,3,4,7), ndmin=2)
		x=np.loadtxt(self.datfile, skiprows=1, usecols=(0,5,6), dtype=np.dtype(str), ndmin=2)
		myetypes=x[:,0]
		(lscore, CNval, t, length, numhists) = range(5)
		Ps=mydat[:,lscore]>pvalcutoff
		Ns=mydat[:,lscore]<=pvalcutoff
		if pvalcutoff>0: 
			Ps=mydat[:,lscore]>pvalcutoff
			Ns=mydat[:,lscore]<=pvalcutoff
			TP = (mydat[:,t] ==1) & Ps
			FP = (mydat[:,t]==0) & Ps
			TN = (Ns & (mydat[:,t] == 0)) | (mydat[:,t]==2)
			FN = (mydat[:,t] ==-1) | ((mydat[:,t] == 1) & Ns)
		else: 
			TP = mydat[:,t] ==1
			FP = mydat[:,t] ==0
			TN = mydat[:,t] == 2 #mydat[:,t] ==0
			FN = mydat[:,t] ==-1
		for i in xrange(len(eventtypes)):
			type =eventtypes[i]
			if type=='any': 
				self.TP[i] = sum(TP)
				self.FP[i] = sum(FP) 
				self.TN[i] = sum(TN)
				self.FN[i] = sum(FN)
			else: 
				self.TP[i] = sum(TP & (myetypes == type))
				self.FP[i] = sum(FP & (myetypes == type))
				self.TN[i] = sum(TN & (myetypes == type))
				self.FN[i] = sum(FN & (myetypes == type))
		
def splitCommaToFloats(s): 
	return (map(float, s.split(',')))

def get_accuracy_values(mysimulations, edgeorevent, eventtypes, numsims, pvalcutoff): 
	datalist=[]
	truenodecounts=[]
	nodecounts=[]
	blocksarray=[]
	simidarray=[]
	costdiffs=[]
	avescores=[]
	minscores=[]
	truescores=[]
	for i in xrange(len(mysimulations)): 
		line=mysimulations[i]
		if len(line.strip().split('\t'))==3:
			(dir, blocks, events)=line.strip().split('\t')
		else: 
			(dir, simid, blocks, events)=line.strip().split('\t')
		mysimdata=SimData(dir, edgeorevent, pvalcutoff)
		datalist.append(mysimdata)
		nodecounts.append(mysimdata.totalnodes)
		truenodecounts.append(mysimdata.truenodes)
		costdiffs.append(mysimdata.costdiff)
		avescores.append(mysimdata.avescore)
		minscores.append(mysimdata.minscore)
		truescores.append(mysimdata.truescore)
		blocksarray.append(blocks)
		simidarray.append(simid)
	myaccuracies={}
	myaccuracies['totalnodes']=nodecounts
	myaccuracies['truenodes']=truenodecounts
	myaccuracies['costdiffs']=costdiffs
	myaccuracies['avescores']=avescores
	myaccuracies['minscores']=minscores
	myaccuracies['truescores']=truescores
	myaccuracies['blocks']=blocksarray
	myaccuracies['simids']=simidarray
	ntypes=len(eventtypes)
	truecounts=range(ntypes)
	totalcounts=range(ntypes)
	f1scores=range(ntypes)
	precisions=range(ntypes)
	recalls=range(ntypes)
	accuracies=range(ntypes)
	adjaccuracies=range(ntypes)
	types=eventtypes
	for type in types: 
		truecount=np.zeros(numsims)
		total=np.zeros(numsims)
		f1score=np.zeros(numsims)
		precision=np.zeros(numsims)
		recall=np.zeros(numsims)
		accuracy=np.zeros(numsims)
		adjaccuracy=np.zeros(numsims)
		for i in xrange(numsims):
			data=datalist[i]	
			truecount[i]=data.TP[type]+data.FN[type]
			total[i]= truecount[i]+data.TN[type]+data.FP[type]
			f1score[i]=float(2*data.TP[type])/float(2*data.TP[type]+data.FP[type] + data.FN[type])
			precision[i]=float(data.TP[type])/float(data.TP[type]+data.FP[type])
			if ((data.TP[type] + data.FN[type]) >0):	
				recall[i]=float(data.TP[type])/float(data.TP[type]+data.FN[type])
			else:
				recall[i]=0
			adjaccuracy[i]=float(data.TP[type])/float(total[i])
			accuracy[i]=float(data.TP[type]+data.TN[type])/float(total[i])
			#else: 
			#	f1score[i]=-1
			#	precision[i]=-1
			#	recall[i]=-1
			#	accuracy[i]=-1
			#	adjaccuracy[i]=-1
		truecounts[type]=truecount
		totalcounts[type]=total
		f1scores[type]=f1score
		precisions[type]=precision
		recalls[type]=recall
		accuracies[type]=accuracy
		adjaccuracies[type]=adjaccuracy
	myaccuracies["truecount"]=truecounts
	myaccuracies["totalcount"]=totalcounts
	myaccuracies["f1score"]=f1scores
	myaccuracies["precision"]=precisions
	myaccuracies["recall"]=recalls
	myaccuracies["accuracy"]=accuracies
	myaccuracies["adjaccuracy"]=adjaccuracies
	return(myaccuracies)

def plot_myaccuracies(ykeys, xkey, types, xlimit, edgeorevent, outname, myaccuracies): 
	numsubplots=len(ykeys)
	fig=plt.figure(figsize=(8,16))
	typcolors=('grey', 'turquoise', 'purple')
	blocks=np.array(myaccuracies['blocks'])
	blockvals=[set(blocks)]
	colorByBlocks=True
	xlabels={'trueconnectivity': 'True connectivity','truecount': 'True %s count' % edgeorevent, 'truenodes': 'true node count', 'connectivity': 'connectivity', 'totalcount': 'Total %s count' % edgeorevent, 'totalnodes': 'Total node count'} 
	xlabel=xlabels[xkey]
	myaxes=[]
	axi=1
	for key in ykeys:
		if axi==1:
			ax=fig.add_subplot(numsubplots, 1, axi)
		else: 
			ax=fig.add_subplot(numsubplots, 1, axi, sharex=myaxes[0], sharey=myaxes[0])
		ax.set_ylabel(key)
		ax.set_xlabel(xlabel)	
		myaxes.append(ax)
		axi+=1
	(xmin, xmax)=myaxes[0].get_xlim()
	if xlimit < xmax: 
		myaxes[0].set_xlim([-1, xlimit])
	myaxes[0].set_ylim([-.1, 1.1])
	if colorByBlocks:
		typcolors=('black', 'red','orange', 'turquoise')
		blocks=np.array(myaccuracies['blocks'], dtype=int)
		blockvals=map(int, list(set(blocks)))
		blockvals.sort()
		color_labels=map(str, blockvals)
		type=0
		if xkey=='trueconnectivity': 
			xdata=np.array(myaccuracies["truecount"])/np.array(myaccuracies["truenodes"])
		elif xkey=='connectivity': 
			xdata=np.array(myaccuracies["totalcount"])/np.array(myaccuracies["totalnodes"])
		else: 
			xdata=np.array(myaccuracies[xkey], ndmin=2)
		xdata=xdata[0]
		for i in xrange(len(blockvals)): 
			sys.stderr.write("Color labels: %s\t%s\n" % (color_labels[i], typcolors[i]))
			blocki=np.where(blocks==blockvals[i])
			sys.stderr.write("blocki: %s\n" % str(blocki))
			axi=0
			for key in ykeys:
				ydata=myaccuracies[key][type]
				myaxes[axi].scatter(xdata[blocki], ydata[blocki], color=typcolors[i], label=color_labels[i])
				myaxes[axi].legend()	
				axi+=1	
	else: 
		for type in types:
			xdata=np.array(myaccuracies["truecount"])/np.array(myaccuracies["truenodecounts"])
		axi=0
		for key in ykeys:
			ydata=myaccuracies[key][type]
			myaxes[axi].scatter(xdata, ydata, color=typcolors[type], label=Global_types[type])
			if len(types)>1: 
				myaxes[axi].legend()	
			axi+=1	
	plt.savefig(outname)


#@interactive 
def main(simulations, edgeorevent, eventtypes, plotfn, datafn, xlimit, xaxis, pvalcutoff):
	mysimulations=open(simulations, 'r').readlines()
	numsims=len(mysimulations)
	ntypes=4
	# myaccuracies will be a dictionary.  Keys: ("truecounts", "f1scores", "precisions", "recalls", "adjaccuracies"), values: arrays of len(mysimulations). 
	myaccuracies=get_accuracy_values(mysimulations, edgeorevent, eventtypes, numsims, pvalcutoff)
	if datafn:
		for type in eventtypes:
			mydata=np.zeros((numsims, len(myaccuracies.keys())-1)) 
			mylabels=[]
			myfmts=[]
			i=0	
			for key in myaccuracies.keys():
				if key != "simids": 
					if key in ["truenodes", "totalnodes", "blocks", "costdiffs", "truescores", "avescores", "minscores"]:
						vals=myaccuracies[key]
					else: 
						vals=myaccuracies[key][type]
					mydata[:,i]=vals
					mylabels.append(key)
					myfmts.append('%f')
					i+=1
			#mydata=np.hstack((np.array(myaccuracies['simids'], ndmin=2).T, mydata))
			#mylabels=['simids'] + mylabels
			#myfmts = ['%s'] + myfmts 
			np.savetxt("%s.%d.txt" % (datafn, type), mydata, header="\t".join(mylabels), fmt=myfmts, delimiter="\t")
	if plotfn:
		ykeys=("recall", "precision", "f1score")
		xkey=xaxis
		plot_myaccuracies(ykeys, xkey, eventtypes, xlimit, edgeorevent, plotfn, myaccuracies)


if __name__ == "__main__": 
	parser=argparse.ArgumentParser(description="makes plots")
	parser.add_argument('simulations', help='A list of simulation directories')
	parser.add_argument('--plot', help='file to print a plot to')
	parser.add_argument('--data', help='File to print data to.')
	parser.add_argument('--xlimit', help='The max value on the x axis', default=500, type=int)
	parser.add_argument('--eventtypes', help='The type of events to include (amp, del, ins)')
	parser.add_argument('--edge', help='whether to analyze data for edges rather than events', action='store_true')
	parser.add_argument('--pvalcutoff', help='the p-value cutoff for what counts as true', default=0, type=float)
	parser.add_argument('--xaxis', help='what you want the x axis to be', choices=['trueconnectivity', 'truecount', 'truenodes', 'connectivity', 'totalcount', 'totalnodes'], default='connectivity')
	args=parser.parse_args()
	if args.edge: 
		edgeorevent="edges"
	else: 
		edgeorevent="events"
	if args.eventtypes: 
		eventtypes=map(int, args.eventtypes.strip().split())
	else:
		 eventtypes=[0]
	main(args.simulations, edgeorevent, eventtypes, args.plot, args.data, args.xlimit, args.xaxis, args.pvalcutoff)


