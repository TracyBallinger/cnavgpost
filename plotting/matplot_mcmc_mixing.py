#!/inside/home/common/bin/python2.7

import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import argparse 
import os, sys
import event_cycles_module as histseg

# from IPython.parallel.util import interactive 

#@interactive 
def plot_costs(historyScores, runs, hids, binwidth):
	hid1=min(hids)
	hid2=max(hids) 
	x=hids
	histseg.Global_BINWIDTH=binwidth
	indices = histseg.historyids_to_indices(hids, historyScores)
	avecosts=np.zeros(len(runs))
	costi=0
	for i in runs:
		indices = histseg.historyids_to_indices(hids+(i*binwidth), historyScores)
		plt.plot(x, np.mean(historyScores[indices,1:3], axis=1), color="grey")
		avecosts[costi]=np.mean(historyScores[indices,1:3])
		costi+=1
		# add plotting function to plot the upper and lower cost bounds with fill inbetween
		plt.fill_between(x, historyScores[indices,1], historyScores[indices,2], facecolor="blue", alpha=0.3)
	plt.xlabel("Iteration")
	plt.ylabel("Mean History Cost")		 
	return np.mean(avecosts)

#@interactive 
def plot_mixing(depdat, inddat, title): 
	plt.ylim(0,1.2)
	x=depdat[:,0]-depdat[:,1]
	for i in xrange(2,inddat.shape[1]):
		plt.plot(x, inddat[:,i], color="pink", alpha=1)
	for i in xrange(2,depdat.shape[1]):
		plt.plot(x, depdat[:,i], color="grey", alpha=1)
	plt.plot(x, np.mean(depdat[:,2:depdat.shape[1]], axis=1), color="black", linewidth=2)
	plt.plot(x, np.mean(inddat[:,2:inddat.shape[1]], axis=1), color="red", linewidth=2)
	plt.xlabel("Number of Iterations")
	plt.ylabel(title + "\nfraction Overlap", multialignment="center")
	#plt.title(title)
	(t, y)=get_time_to_indep(depdat, inddat)
	plt.vlines(t, 0, y, color="green", linewidth=2, zorder=3) 

def plot_ind_mixing(inddat, title, binwidth): 
	plt.ylim(0,1.2)
	x=np.fmod(inddat[:,0], binwidth) 
	for i in xrange(2,inddat.shape[1]):
		plt.plot(x, inddat[:,i], color="pink", alpha=1)
	plt.plot(x, np.mean(inddat[:,2:inddat.shape[1]], axis=1), color="red", linewidth=2)
	aveind=np.mean(inddat[:,2:inddat.shape[1]])
	plt.hlines(aveind, min(x), max(x), color="blue", linewidth=2, zorder=3)
	plt.xlabel("Number of Iterations")
	plt.ylabel(title + "\nIndependent history\n overlaps", multialignment="center")
	return aveind 

def plot_sim_mixing(simdat, title, binwidth): 
	plt.ylim(0,1.2)
	x= np.fmod(simdat[:,1], binwidth) 
	for i in xrange(2,simdat.shape[1]):
		plt.plot(x, simdat[:,i], color="pink", alpha=1)
	plt.plot(x, np.mean(simdat[:,2:simdat.shape[1]], axis=1), color="red", linewidth=2)
	plt.xlabel("Iteration")
	plt.ylabel(title + "\nTrue history\n Overlap", multialignment="center")
	return np.mean(simdat[:,2:simdat.shape[1]])

def plot_dep_mixing(depdat, title, aveind): 
	plt.ylim(0,1.2)
	x=depdat[:,1]-depdat[:,0]
	for i in xrange(2,depdat.shape[1]):
		plt.plot(x, depdat[:,i], color="grey", alpha=1)
	avedep=np.mean(depdat[:,2:depdat.shape[1]], axis=1)
	plt.plot(x, avedep, color="black", linewidth=2)
	plt.ylabel(title + "\nWithin run\n Overlap", multialignment="center")
	plt.xlabel("Iteration")
	i=0
	while i < len(avedep) and avedep[i]>=aveind: 
		i+=1
	t=x[i-1]
	plt.hlines(aveind, min(x), max(x), color="blue", linewidth=2, zorder=3)
	#plt.vlines(t, 0, aveind, color="green", linewidth=2, zorder=3)
	return np.mean(avedep)

def get_time_to_indep(depdat, inddat):
	avedep=np.mean(depdat[:,2:depdat.shape[1]], axis=1)
	#aveind=np.mean(inddat[:,2:inddat.shape[1]], axis=1)
	x=np.fmod(inddat[:,0], 100000) 
	x=inddat[:,0] 
	aveind=np.mean(inddat[inddat.shape[0]-1,2:inddat.shape[1]])
	step=depdat[:,1] - depdat[:,0]
	i=0
	#while i<len(avedep) and avedep[i]>=aveind[i]:
	while i<len(avedep) and avedep[i]>=aveind:
		i+=1 
	i+=-1
	return(step[i], aveind)
	

#sampledir="TCGA-06-0145/"
#@interactive 
def main(sampledir, outfigure, plot, args):
	historyScores=np.loadtxt(os.path.join(sampledir, "historystats.txt"), dtype=int)
	mcmcdatadir="mcmcdata"
	if args.edges: 
		edgecounts=np.loadtxt(os.path.join(sampledir, mcmcdatadir, "edge_counts.dat"))
		indedge=np.loadtxt(os.path.join(sampledir, mcmcdatadir, "edgemixing_ind.dat"))
		depedge=np.loadtxt(os.path.join(sampledir, mcmcdatadir, "edgemixing_dep.dat"))
	evntcounts=np.loadtxt(os.path.join(sampledir, mcmcdatadir, "evnt_counts.dat"))
	indevnt=np.loadtxt(os.path.join(sampledir, mcmcdatadir, "evntmixing_ind.dat"))
	depevnt=np.loadtxt(os.path.join(sampledir, mcmcdatadir, "evntmixing_dep.dat"))
	numruns=depevnt.shape[1]-2
	numiters=depevnt.shape[0]
	hids=indevnt[:,0]
	hid1=int(indevnt[0,1])
	binwidth=10 ** (len(list(str(hid1)))-1)
	sys.stderr.write("binwidth: %d, hid1: %d\n" % (binwidth, hid1))
	hids=np.fmod(hids, binwidth)
	if args.simulation: 
		if args.edges: 
			simedge=np.loadtxt(os.path.join(sampledir, mcmcdatadir, "edgemixing_sim.dat"))
		simevnt=np.loadtxt(os.path.join(sampledir, mcmcdatadir, "evntmixing_sim.dat"))
	(ave_cost, numevnt, avesim_evnt, aveind_evnt, avedep_evnt, numedge, avesim_edge, aveind_edge, avedep_edge) = [0]*9
	if args.edges: 
		numedge=np.mean(edgecounts[0,:])
	numevnt=np.mean(evntcounts[0,:])
	if plot:
		if not args.simulation: 
			if args.edges:
				fig=plt.figure(figsize=(12,15))
				numplots=5
			else: 
				fig=plt.figure(figsize=(12,8))
				numplots=3
			plt.subplot(numplots,1,1)
			ave_cost=plot_costs(historyScores, range(numruns), hids, binwidth) 
			#plt.xlim(0,2500)
			plt.subplot(numplots,1,2)
			aveind_evnt=plot_ind_mixing(indevnt, "Events", binwidth)
			#plt.xlim(0,2500)
			plt.subplot(numplots,1,3)
			avedep_evnt=plot_dep_mixing(depevnt, "Events", aveind_evnt)
			#plt.xlim(-2500,0)
			if args.edges: 
				plt.subplot(numplots,1,4)
				aveind_edge=plot_ind_mixing(indevnt, "Edges", binwidth)
				plt.subplot(numplots,1,5)
				avedep_edge=plot_dep_mixing(depevnt, "Edges", aveind)
			plt.savefig(outfigure)
		if args.simulation:
			if args.edges: 
				fig=plt.figure(figsize=(12,22))
				numplots=7
			else: 
				fig=plt.figure(figsize=(12,15))
				numplots=5
			plt.subplot(numplots,1,1)
			ave_cost=plot_costs(historyScores, range(1,numruns), hids, binwidth) 
			plt.subplot(numplots,1,2)
			avesim_evnt=plot_sim_mixing(simevnt, "Events", binwidth)
			plt.subplot(numplots,1,3)
			aveind_evnt=plot_ind_mixing(indevnt, "Events", binwidth)
			plt.subplot(numplots,1,4)
			avedep_evnt=plot_dep_mixing(depevnt, "Events", aveind_evnt)
			if args.edges: 
				plt.subplot(numplots,1,5)
				avesim_edge=plot_sim_mixing(simedge, "Edges", binwidth)
				plt.subplot(numplots,1,6)
				aveind_edge=plot_ind_mixing(indedge, "Edges", binwidth)
				plt.subplot(numplots,1,7)
				avedep_edge=plot_dep_mixing(depedge, "Edges", aveind_edge)
			plt.savefig(outfigure)
			
	sys.stdout.write("\t".join(map(str,(ave_cost, numevnt, avesim_evnt, aveind_evnt, avedep_evnt, numedge, avesim_edge, aveind_edge, avedep_edge))) + "\n")


if __name__ == "__main__": 
	parser=argparse.ArgumentParser(description="makes plots and prints out the following data: <ave_cost><time_to_independence><overlap at independence><average number of events><time_to_independence><overlap at independence><average number of edges>.  The averages are calculated using the final history across the 10 independent runs.")
	parser.add_argument('sampledir', help='The directory where the data files are')
	parser.add_argument('outfigure', help='The file for the plots that are made')
	parser.add_argument('--plot', help='whether to plot a figure or not.', action='store_true')
	parser.add_argument('--simulation', help='the data is from a simulation.', action='store_true')
	parser.add_argument('--edges', help='analyze the edges also.', action='store_true')
	args=parser.parse_args()
	main(args.sampledir, args.outfigure, args.plot, args)


