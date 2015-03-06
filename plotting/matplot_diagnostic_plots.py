#!/usr/bin/env python 
import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import argparse 
import os, sys
import matplotlib.gridspec as gridspec

#from IPython.parallel.util import interactive 

#@interactive 
def plot_costs(historyScores, numruns, id1=0, id2=0, myax=plt.gca()):
	numiters=historyScores.shape[0]/numruns
	hid1=min(id1, id2)
	hid2=max(id1, id2)
	if hid2==0: 
		hid2=numiters
	x=range(hid1, hid2)
	runs=range(numruns)
	sys.stderr.write("hid1: %d, hid2: %d, numiters: %d\n" % (hid1, hid2, numiters))
	for i in runs:
		s=i*numiters+hid1
		e=i*numiters+hid2
		yvals=np.mean(historyScores[s:e,1:3], axis=1)	
		# add plotting function to plot the upper and lower cost bounds with fill inbetween
		myax.fill_between(x, historyScores[s:e,1], historyScores[s:e,2], facecolor="blue", alpha=0.1)
		myax.plot(x, np.mean(historyScores[s:e,1:3], axis=1), color="black")
	myax.set_xlim([hid1, hid2])
	plt.xlabel("Iteration")
	plt.ylabel("Average History Cost")		 

#@interactive 
def plot_mixing(depdat, inddat, title, myax=plt.gca()): 
	plt.ylim(0,1.2)
	x=depdat[:,1]
	for i in xrange(2,inddat.shape[1]):
		myax.plot(x, inddat[:,i], color="pink", alpha=0.2)
	for i in xrange(2,depdat.shape[1]):
		myax.plot(x, depdat[:,i], color="grey", alpha=0.2)
	depline, =myax.plot(x, np.mean(depdat[:,2:depdat.shape[1]], axis=1), color="black", linewidth=2, label='Similarity to the final history within a single run')
	indline, =myax.plot(x, np.mean(inddat[:,2:inddat.shape[1]], axis=1), color="red", linewidth=2, label='Similarity between independent sampling runs')
	plt.xlabel("Iteration")
	plt.ylabel(title + "\nSimilarity", multialignment="center")
	leg=myax.legend([depline, indline], prop={'size':8})
	leg.draw_frame(False)
	#plt.title(title)
	#(t, y)=get_time_to_indep(depdat, inddat)
	#plt.vlines(t, 0, y, color="green", linewidth=2, zorder=3) 

def plot_sim_mixing(simdat, title): 
	plt.ylim(0,1.2)
	x= np.fmod(simdat[:,1], 10000) 
	for i in xrange(2,simdat.shape[1]):
		plt.plot(x, simdat[:,i], color="pink", alpha=1)
	plt.plot(x, np.mean(simdat[:,2:simdat.shape[1]], axis=1), color="red", linewidth=2)
	plt.xlabel("Iteration")
	plt.ylabel(title + "\nSimilarity to\ntrue history", multialignment="center")

def plot_ind_mixing(inddat, title): 
	plt.ylim(0,1.2)
	x=np.fmod(inddat[:,0], 10000) 
	sys.stderr.write("ind %s x[1:10]: %s\n" % (title, str(x[1:10])))
	for i in xrange(2,inddat.shape[1]):
		plt.plot(x, inddat[:,i], color="pink", alpha=1)
	#plt.plot(x, np.mean(inddat[:,2:inddat.shape[1]], axis=1), color="red", linewidth=2)
	avey=np.mean(inddat[inddat.shape[0]-1,2:inddat.shape[1]])
	plt.hlines(avey, min(x), max(x), color="blue", linewidth=2, zorder=3)
	plt.xlabel("Number of Iterations")
	plt.ylabel("Independent runs\n" + title + " Overlap", multialignment="center")
	return avey 

def plot_dep_mixing(depdat, title, aveind): 
	plt.ylim(0,1.2)
	x=depdat[:,1]-depdat[:,0]
	sys.stderr.write("dep %s x[1:10]: %s\n" % (title, str(x[1:10])))
	for i in xrange(2,depdat.shape[1]):
		plt.plot(x, depdat[:,i], color="grey", alpha=1)
	avedep=np.mean(depdat[:,2:depdat.shape[1]], axis=1)
	#plt.plot(x, avedep, color="black", linewidth=2)
	plt.ylabel("Within run\n" + title + " Overlap", multialignment="center")
	plt.xlabel("Iteration")
	i=0
	while i < len(avedep) and avedep[i]>=aveind: 
		i+=1
	t=x[i-1]
	plt.hlines(aveind, min(x), max(x), color="blue", linewidth=2, zorder=3)
	#plt.vlines(t, 0, aveind, color="green", linewidth=2, zorder=3)
	return t

def plot_lscore_histogram(sampledir, edgeld=None, ax=plt.gca()): 
	zerocounts=[0,0,0]
	myylims=[0,0,0]
	mycolors=['blue', 'green', 'orange']
	eventdatfn=os.path.join(sampledir, "events.dat")
	if os.path.exists(eventdatfn): 
		i=0
		eventscores=np.loadtxt(eventdatfn, usecols=4, skiprows=1, dtype=float)
		n, bins, patches = ax.hist(eventscores, histtype='step', linewidth=2, label='events', color=mycolors[i])
		myylims[i]= 3*max(sorted(n)[1:])
		zerocounts[i]=n[0]
	edgedatfn=os.path.join(sampledir, "edges.dat")
	if os.path.exists(edgedatfn):
		i=1
		edgescores=np.loadtxt(edgedatfn, usecols=4, skiprows=1, dtype=float)
		n, bins, patches = ax.hist(edgescores, histtype='step', linewidth=2, label='edges', color=mycolors[i])
		myylims[i]= 3*max(sorted(n)[1:])  # the second highest value
		zerocounts[i]=n[0] 
	meventdatfn=os.path.join(sampledir, "mrgevents.dat")
	if os.path.exists(meventdatfn):
		i=2
		mevntscores=np.loadtxt(meventdatfn, usecols=4, skiprows=1, dtype=float)
		n, bins, patches = ax.hist(mevntscores, histtype='step', linewidth=2, label='merged events', color=mycolors[i])
		myylims[i]= 3*max(sorted(n)[1:])  # the second highest value
		zerocounts[i] = n[0]
	myylim=max(myylims)
	for i in xrange(len(zerocounts)):
		if zerocounts[i] >0: 
			plt.text((bins[1]-bins[0])/2, myylim*float(zerocounts[i])/float(max(zerocounts)) , int(zerocounts[i]), va='top', ha='center', color=mycolors[i], size=9)
	plt.ylim(0, myylim)
	if len(zerocounts) >0: 
		leg=ax.legend()
		leg.draw_frame(False)
	plt.xlabel("Likelihood Score")
	plt.ylabel("Number of Events")
	
def plot_lota_various_k_or_n(ld, nk, myax=plt.gca()): 
	sys.stderr.write("nk is %s\n" % nk)
	if nk=="n": 
		for i in xrange(ld.countsTPn.shape[1]):
			myax.plot(ld.countsTPn[:,i], linewidth=2, label=ld.nlabels[i])  
		leg=myax.legend(title="Number of runs", prop={'size':8}, ncol=2)
	else: 
		for i in xrange(ld.countsTPk.shape[1]):
			myax.plot(ld.countsTPk[:,i], linewidth=2, label=ld.klabels[i])  
		leg=myax.legend(title="Value of k", prop={'size':8})
		sys.stderr.write("leg: %s\n" % leg) 
	leg.draw_frame(False)
	plt.xlabel("iteration")
	plt.ylabel("Number of events \nwith likelihood > 0.5")
	
def plot_lota_data(lotadata, ldedges=None, gs=""): 
	ld=lotadata
	#histogram of the likelihood scores
	if gs !="":
		gs1=gridspec.GridSpec(1,1, left=gs.left, right=gs.right)
		gs2=gridspec.GridSpec(2,1, top=0.6, hspace=0.1, left=gs.left, right=gs.right)
	else: 
		gs1=gridspec.GridSpec(1,1)
		gs2=gridspec.GridSpec(2,1, top=0.6, hspace=0.1)
	gs1.update(bottom=0.70)
	ax1=plt.subplot(gs1[0,0])
	plot_lscore_histogram(lotadata, ldedges, ax1)	
	# plot number of events with various number of runs	
	ax2=plt.subplot(gs2[0,0])
	plot_lota_various_k_or_n(lotadata, "n", ax2) 
	ax2.get_xaxis().set_visible(False)
	# plot number of events with various values of k 
	ax3=plt.subplot(gs2[1,0], sharex=ax2)
	plot_lota_various_k_or_n(lotadata, "k", ax3) 

def multiplot_mcmc(data, simevnt=None, simedge=None, gs=""):
	numplots=3
	if simevnt is not None and simedge is not None: 
		numplots=5
	gs1=gs
	if gs=="": 
		gs1=gridspec.GridSpec(numplots,1)
	gs1.update(hspace=0.1, wspace=0.6)
	ax1=plt.subplot(gs1[0,0]) 
	plot_costs(data.historyScores, data.numruns, myax=ax1) 
	ax1.get_xaxis().set_visible(False)
	ax2=plt.subplot(gs1[1,0], sharex=ax1)
	plot_mixing(data.depevnt, data.indevnt, "Events", ax2)
	ax2.get_xaxis().set_visible(False)
	ax3=plt.subplot(gs1[2,0], sharex=ax1)
	plot_mixing(data.depedge, data.indedge, "Edges", ax3)
	if simevnt is not None and simedge is not None: 
		ax3.get_xaxis().set_visible(False)
		ax4=plt.subplot(gs1[3,0], sharex=ax1)
		plot_sim_mixing(simevnt, "Events")
		ax4.get_xaxis().set_visible(False)
		ax5=plt.subplot(gs1[4,0], sharex=ax1)
		plot_sim_mixing(simedge, "Edges")

def plot_mcmc_and_lota_data(mcmcdata, lotadata, ldedge=None, simevnt=None, simedge=None): 
	numplots=3
	if (simevnt is not None) and (simedge is not None): 
		numplots=5
	gs1=gridspec.GridSpec(numplots, 1, right=0.45)
	multiplot_mcmc(mcmcdata, gs=gs1)
	gs2=gridspec.GridSpec(1,1, left=0.55)
	plot_lota_data(lotadata, ldedge, gs=gs2) 
	
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
	
class mcmcdata: 
	def __init__(self, sampledir): 
		self.historyScores=np.loadtxt(os.path.join(sampledir, "historystats.txt"))
		self.edgecounts=np.loadtxt(os.path.join(sampledir, "mcmcdata/edge_counts.dat"))
		self.evntcounts=np.loadtxt(os.path.join(sampledir, "mcmcdata/evnt_counts.dat"))
		self.indedge=np.loadtxt(os.path.join(sampledir, "mcmcdata/edgemixing_ind.dat"))
		self.depedge=np.loadtxt(os.path.join(sampledir, "mcmcdata/edgemixing_dep.dat"))
		self.indevnt=np.loadtxt(os.path.join(sampledir, "mcmcdata/evntmixing_ind.dat"))
		self.depevnt=np.loadtxt(os.path.join(sampledir, "mcmcdata/evntmixing_dep.dat"))
		self.numruns=self.depedge.shape[1]-2
	
class lotadata:
	def __init__(self, lotadir):
		self.countsTPk=np.loadtxt(os.path.join(lotadir, "counts_TP_ks.dat"))
		self.countsTotk=np.loadtxt(os.path.join(lotadir, "counts_Tot_ks.dat"))
		self.countsTPn=np.loadtxt(os.path.join(lotadir, "counts_TP_ns.dat"))
		self.countsTotn=np.loadtxt(os.path.join(lotadir, "counts_Tot_ns.dat"))
		headerline=open(os.path.join(lotadir, "counts_Tot_ks.dat")).readline()
		self.klabels=headerline[1:].split()
		headerline=open(os.path.join(lotadir, "counts_Tot_ns.dat")).readline()
		self.nlabels=headerline[1:].split()
		self.lscores=np.loadtxt(os.path.join(lotadir, "lscores_0.000000.txt.gz"))


#sampledir="TCGA-06-0145/"
#@interactive 
def main(sampledir, outfigure, plot):
	md=None
	ld=None
	ld2=None
	simevnt=None
	simedge=None
	if os.path.exists(os.path.join(sampledir, "mcmcdata")): 
		md = mcmcdata(sampledir)
		(edge_t, edge_y)=get_time_to_indep(md.depedge, md.indedge)
		(event_t, event_y)=get_time_to_indep(md.depevnt, md.indevnt)
		ave_cost=0 #get_average_costs(historyScores, numruns, numiters) #np.mean(costs[(costs.shape[0]-1),:])
		numedge=np.mean(md.edgecounts[0,:])
		numevnt=np.mean(md.evntcounts[0,:])
		sys.stdout.write("\t".join(map(str,(ave_cost, event_t, event_y, numevnt, edge_t, edge_y, numedge))) + "\n")
	ldir=os.path.join(sampledir, "lotadir")
	if os.path.exists(ldir):
		ld = lotadata(ldir)
	l2dir=os.path.join(sampledir, "lota_edge")
	if os.path.exists(l2dir):
		ld2 = lotadata(l2dir)
	if os.path.isfile(os.path.join(sampledir, "mcmcdata/edgemixing_sim.dat")):
		simedge=np.loadtxt(os.path.join(sampledir, "mcmcdata/edgemixing_sim.dat"))
	if os.path.isfile(os.path.join(sampledir, "mcmcdata/evntmixing_sim.dat")):
		simevnt=np.loadtxt(os.path.join(sampledir, "mcmcdata/evntmixing_sim.dat"))
		simulation=True
	if plot:
		fig=plt.figure(figsize=(14,13))
		if (ld is not None) and (md is not None): 
			sys.stderr.write("testing plot now\n")
		#	plot_lota_various_k_or_n(ld, "n") 
			plot_mcmc_and_lota_data(md, ld, ld2, simevnt, simedge) 
		elif ld is not None:
			plot_lota_data(ld, ld2)
			#plot_lota_various_k_or_n(ld, "n") 
		elif md is not None: 
			multiplot_mcmc(md, simevnt, simedge)
		plt.savefig(outfigure)
			

if __name__ == "__main__": 
	parser=argparse.ArgumentParser(description="makes plots and prints out the following data: <ave_cost><time_to_independence><overlap at independence><average number of events><time_to_independence><overlap at independence><average number of edges>.  The averages are calculated using the final history across the 10 independent runs.")
	parser.add_argument('sampledir', help='The directory where the data files are')
	parser.add_argument('outfigure', help='The file for the plots that are made')
	parser.add_argument('--plot', help='whether to plot a figure or not.', action='store_true')
	args=parser.parse_args()
	main(args.sampledir, args.outfigure, args.plot)


