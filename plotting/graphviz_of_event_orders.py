#!/usr/bin/env python 
import matplotlib
matplotlib.use('Agg')
import argparse
import sys, os
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec

class EventsData: 
	def __init__(self, filename):
		header="event_id	event_type	avecost	Lscore	CNval	true	length	prevalmean	prevalsd	ordermean	ordersd	numhists"
		(eventid, etype, avecost, lscore, cnval, true, length, prevalmean, prevalsd, ordermean, ordersd, numhist) = range(len(header.split('\t')))
		data=np.loadtxt(filename, skiprows=1, dtype=str)
		self.ids=data[:,eventid]
		self.etype=data[:,etype]
		self.avecost=data[:,avecost].astype(float) 
		self.lscore=data[:,lscore].astype(float) 
		self.lscore=self.lscore/max(self.lscore)
		self.cnval=data[:,cnval].astype(float) 
		self.length=data[:,length].astype(int) 
		self.preval=data[:,prevalmean].astype(float) 
		self.prevalsd=data[:,prevalsd].astype(float)
		self.order=data[:,ordermean].astype(float) 
		self.ordersd=data[:,ordersd].astype(float) 
		self.numhist=data[:,numhist].astype(int) 

class LinksData: 
	def __init__(self, filename): 
		data=np.loadtxt(filename, dtype=str)
		self.ids1=data[:,0]
		self.ids2=data[:,1]
		self.numhists=data[:,2].astype(int)
		self.lscores=data[:,3].astype(float)
		#scale the lscores
		self.lscores=self.lscores/max(self.lscores)
		self.x1=np.ones(len(self.ids1))
		self.y1=np.ones(len(self.ids1))
		self.x2=np.ones(len(self.ids1))
		self.y2=np.ones(len(self.ids1))
		
	def get_coords(self, events, orderOrPreval):
		if orderOrPreval=='preval': 
			for i in xrange(len(self.ids1)): 
				myi=np.where(events.ids==self.ids1[i])
				self.x1[i]=events.preval[myi]
				self.y1[i]=events.cnval[myi]
				myi=np.where(events.ids==self.ids2[i])
				self.x2[i]=events.preval[myi]
				self.y2[i]=events.cnval[myi]
		else: 
			for i in xrange(len(self.ids1)): 
				myi=np.where(events.ids==self.ids1[i])
				self.x1[i]=events.order[myi]
				self.y1[i]=events.cnval[myi]
				myi=np.where(events.ids==self.ids2[i])
				self.x2[i]=events.order[myi]
				self.y2[i]=events.cnval[myi]
		

def scatter_plot_prevalence(events, links):
	etype=events.etype
	myi=np.where(etype=='adj')
	plt.scatter(1-events.preval[myi], events.cnval[myi], s=events.lscore[myi]*100, c='yellow', zorder=2, label='adj')
	myi=np.where(etype=='del')
	plt.scatter(1-events.preval[myi], events.cnval[myi], s=events.lscore[myi]*100, c='blue', zorder=2, label='del')
	myi=np.where(etype=='amp')
	plt.scatter(1-events.preval[myi], events.cnval[myi], s=events.lscore[myi]*100, c='red', zorder=2, label='amp')
	links.get_coords(events, 'preval')
	plt.plot(np.vstack((links.x1, links.x2)), np.vstack((links.y1, links.y2)), c='black', zorder=1)
	plt.xlabel("Timing")
	plt.ylabel("CN change")
	leg=plt.legend()
	leg.draw_frame(False)
	#plot(np.vstack((links.x1, links.x2)), np.vstack((links.y1, links.y2)), c='black', linewidth=links.lscores*10)
	#ecols=etype
	#ecols[etype=='adj']='yellow'
	#ecols[etype=='del']='blue'
	#ecols[etype=='amp']='red'

def scatter_plot_order(events, links):
	etype=events.etype
	myi=np.where(etype=='adj')
	plt.scatter(events.order[myi], events.cnval[myi], s=events.lscore[myi]*100, c='yellow', zorder=2, label='adj')
	myi=np.where(etype=='del')
	plt.scatter(events.order[myi], events.cnval[myi], s=events.lscore[myi]*100, c='blue', zorder=2, label='del')
	myi=np.where(etype=='amp')
	plt.scatter(events.order[myi], events.cnval[myi], s=events.lscore[myi]*100, c='red', zorder=2, label='amp')
	links.get_coords(events, 'order')
	plt.plot(np.vstack((links.x1, links.x2)), np.vstack((links.y1, links.y2)), c='black', zorder=1)
	plt.xlabel("Average Order in histories")
	plt.ylabel("CN change")
	leg=plt.legend()
	leg.draw_frame(False)

def scatter_with_annotation(events, links, annfn): 
	myanns=range(len(events.ids))
	genecounts=np.ones(len(events.ids))
	annfh=open(annfn, 'r')
	for line in annfh: 
		(eid, ann)=line.strip().split()
		genes=ann.replace('/', ',').split(',')
		myi=np.where(events.ids==eid)
		if len(myi[0])>0:
			i=myi[0][0] 
			genecounts[i]=len(genes)
			myanns[i]=ann
	annfh.close()		
	#scatter_plot_prevalence(events, links)	
	etype=events.etype
	myi=np.where(etype=='adj')
	plt.scatter(1-events.preval[myi], events.cnval[myi], s=events.lscore[myi]*100, c='yellow', zorder=2, label='adj')
	myi=np.where(etype=='del')
	plt.scatter(1-events.preval[myi], events.cnval[myi], s=events.lscore[myi]*100, c='blue', zorder=2, label='del')
	myi=np.where(etype=='amp')
	plt.scatter(1-events.preval[myi], events.cnval[myi], s=events.lscore[myi]*100, c='red', zorder=2, label='amp')
	links.get_coords(events, 'preval')
	plt.plot(np.vstack((links.x1, links.x2)), np.vstack((links.y1, links.y2)), c='black', zorder=1)
	plt.xlabel("Timing")
	plt.ylabel("CN change")
	leg=plt.legend()
	leg.draw_frame(False)
	for i in np.where(genecounts<=5)[0]: 
		plt.annotate(myanns[i], xy=(1-events.preval[i], events.cnval[i]), fontsize=8, ha='left', va='top', clip_on=True)	


def scatter_with_given_annotation(events, links, annfn, genelist): 
	myanns=[]
	for i in xrange(len(events.ids)): 
		myanns.append([])
	genecounts=np.ones(len(events.ids))
	annfh=open(annfn, 'r')
	for line in annfh: 
		(eid, ann)=line.strip().split()
		genes=ann.replace('/', ',').split(',')
		myi=np.where(events.ids==eid)
		if len(myi[0])>0:
			i=myi[0][0] 
			genecounts[i]=len(genes)
			for g in genes: 
				if g in genelist: 
					myanns[i].append(g)
	annfh.close()		
	#scatter_plot_prevalence(events, links)	
	ax1=plt.subplot(1,2,1)
	etype=events.etype	
	pointsize= events.lscore*100 #genecounts
	myi=np.where(etype=='adj')
	ax1.scatter(events.cnval[myi], events.preval[myi], s=pointsize[myi], c='yellow', zorder=2, label='adj')
	myi=np.where(etype=='del')
	ax1.scatter(events.cnval[myi], events.preval[myi], s=pointsize[myi], c='blue', zorder=2, label='del')
	myi=np.where(etype=='amp')
	ax1.scatter(events.cnval[myi], events.preval[myi], s=pointsize[myi], c='red', zorder=2, label='amp')
	links.get_coords(events, 'preval')
	ax1.plot(np.vstack((links.y1, links.y2)), np.vstack((links.x1, links.x2)), c='black', zorder=1)
	ax1.set_ylabel("Timing")
	ax1.set_xlabel("CN change")
	ax1.set_ylim(0,1.1)
	handles, labels = ax1.get_legend_handles_labels()
	sys.stderr.write("there are %d handles\n" % len(handles))
	#leg=ax1.legend(handles, labels)
	#leg.draw_frame(False)
	# add annotation to the graph
	for i in xrange(len(myanns)):
		ann=myanns[i]
		if len(ann)>0: 
			plt.annotate(",".join(ann), xy=(events.cnval[i], events.preval[i]), fontsize=8)	
	# plot the number of genes per event
	plt.subplot(1,2,2)
	plt.barh(events.preval, genecounts/1000, height=0.01, alpha=0.5, linewidth=0.1)
	plt.xlabel("Number of Genes hit (thousands)")
	plt.ylim(0,1.1)

#def plot_timing_vs_gene_coordinate(edat, ldat, genename=""): 
	


def make_gv_file(evntsfn, linksfn, outgv): 	
	eventf=open(evntsfn, 'r')
	outf=open(outgv, 'w')
	outf.write("digraph G {\n")
	outf.write("size = \"8.5,11!\";\n")
	maxscore=0
	for line in eventf: 
		event=Event(line)
		label=event.id #data[9] #should be gene names
		cnval=event.cnval
		score=event.lscore 
		maxscore=max(score, maxscore)
		if cnval >1: color="red" # "%f 0 0" % (score) #"red"
		elif cnval <0: color="blue" #"0 0 %f" % (score) #"blue"
		else: color="black" #"%f %f %f" % (score, score, score) #"black"
		outf.write("%s [color=%s, label=\"%s\"];\n" % (id, color, label))
	
	linkf=open(linksfn, 'r')
	style="solid"
	arrowhead="normal"
	for line in linkf:
		link=Link(line)
		node1id = link.id1
		node2id=link.id2
		score=(float(link.lscore)/maxscore)*3
		outf.write("%s -> %s [color=\"black\", style=%s, arrowhead=%s, penwidth=%f]\n" % (node1id, node2id, style, arrowhead, score))
	outf.write("}\n")	
		
	
if __name__ == '__main__': 
	parser=argparse.ArgumentParser(description='Given a .evnts file and a .links file containing the links between the events, it will generate a graph of the history.') 
	parser.add_argument('evnts', help='a .evnts dat file')
	parser.add_argument('links', help='a .links dat file')
	parser.add_argument('pdf', help='The pdf of the scatter plot.')
	parser.add_argument('--gv', help='generate a graphing file for graphviz.')
	parser.add_argument('--ann', help='The annotation file for the events')
	parser.add_argument('--genes', help='This is a text file that has a list of the genes to be annotated on the scatter plot.')
	args=parser.parse_args()

	if args.gv: 
		make_gv_file(args.evnts, args.links, args.gv)
	if args.genes: 
		gf=open(args.genes, 'r')
		genes=[]
		for l in gf.readlines(): 
			genes.append(l.strip())
		myevents=EventsData(args.evnts)
		mylinks=LinksData(args.links)
		fig=plt.figure(figsize=(8,11))
		scatter_with_given_annotation(myevents, mylinks, args.ann, genes)
		plt.savefig(args.pdf) 



