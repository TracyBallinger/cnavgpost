#!/usr/bin/env python 

import argparse 
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import matplotlib.patches as mpatches

def plot_sample_panel(samplist):
	samples=open(samplist, 'r').readlines()
	#samples=samplelist
	numplots=len(samples)
	nc=min(4, numplots)
	nr=max(numplots/nc, 0)
	fig, axs=plt.subplots(numplots/4,4, sharex=True, sharey=True, subplot_kw=dict(xlim=(0,4), ylim=(0,1e9)))
	for (i, s) in enumerate(samples):
		r=i/4
		c=i % 4 
		if nc>0:
			plt.sca(axs[r,c])
			plot_CNhist(s.strip(), 0.1)
		else: 
			plt.sca(axs[r])
			plot_CNhist(s.strip(), 0.1)

def get_tableau10(): 
	tableau10= [(31, 119, 180),(255, 127, 14),(44, 160, 44),(214, 39, 40),(148, 103, 189),(140, 86, 75),(227, 119, 194),(127, 127, 127),(188, 189, 34),(23, 190, 207)]
	for i in range(len(tableau10)):
		r, g, b = tableau10[i]  
		tableau10[i] = (r / 255., g / 255., b / 255.)
	return tableau10

def get_tableau20(): 
	# These are the "Tableau 20" colors as RGB.  
	tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),  
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),  
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),  
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),  
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]  
	# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.  
	for i in range(len(tableau20)):  
		r, g, b = tableau20[i]  
		tableau20[i] = (r / 255., g / 255., b / 255.)
	return tableau20

class CNprofile: 
	def __init__(self, segline): 
		dat=segline.strip().split('\t')
		self.chr=dat[0]
		self.start=int(dat[1])
		self.end=int(dat[2])
		self.lscore=float(dat[3])
		self.numhists=int(dat[4])
		self.cnvals=map(float, dat[5].split(','))
		self.prevals=map(float, dat[6].split(','))
		self.prevalsd=map(float, dat[7].split(','))
		self.label=""
		if len(dat)==9: 
			self.label=dat[8]

def group_seglines_by_loc(seglines):
	mylocs={}
	for line in seglines:
		cnp=CNprofile(line)
		myloc="%s:%d-%d" % (cnp.chr, cnp.start, cnp.end)
		if myloc in mylocs.keys(): 
			mylocs[myloc].append(line)
		else: 
			mylocs[myloc]=[line]
	return mylocs	

def plot_seglines_by_loc(seglines): 
	mylocs=group_seglines_by_loc(seglines)
	mycolors=get_tableau20()
	for (i, loc) in enumerate(mylocs.keys()): 
		plot_seghists(mylocs[loc], color=mycolors[i])

def plot_seglines_by_label(seglines, genes, plotnone=True):
	genecolors={}
	mycolors=get_tableau10()
	for (i,g) in enumerate(genes): 
		genecolors[g]=mycolors[i]
	geneprofs={'None':[]}
	for g in genes: 
		geneprofs[g]=[]
	for line in seglines:
		cnp=CNprofile(line) 
		hasgene=False
		if cnp.label != "":
			genenames=cnp.label.split(",")
			for g in genes: 
				if g in genenames: 
					geneprofs[g].append(cnp)
					hasgene=True
		if not hasgene: 
			geneprofs['None'].append(cnp)
	if plotnone: plot_cnps(geneprofs['None'], 'gray')
	print "gene %s, %d" % ('None', len(geneprofs['None']))
	leghandles=[]
	for g in genes:
		print "gene %s, %d" % (g, len(geneprofs[g]))
		plot_cnps(geneprofs[g], genecolors[g])
		leghandles.append(mpatches.Patch(color=genecolors[g], label="%s (%d)" % (g, len(geneprofs[g]))))
	plt.xlim(1.1,0)
	plt.ylabel("copy number", fontsize=20)
	plt.xlabel("prevalence", fontsize=20)
	# add legend
	plt.legend(handles=leghandles)
	
def plot_cnps(cnps, color): 
	for cnp in cnps: 
		plot_seghist(cnp, color)

def plot_seghists(seglines, color): 
	for line in seglines: 
		plot_seghist(CNprofile(line), color)

def plot_top_seghists(seglines, color): 
	mylines=[]
	cutoff=0.01
	for line in seglines:
		cnp=CNprofile(line)
		if cnp.lscore > cutoff: 
			mylines.append(cnp)
	print "There are %d lines" % len(mylines)
	if len(mylines)>10: 
		mylines=sorted(mylines, key=lambda x: (x.lscore), reverse=True)[:10]
	print "There are %d lines" % len(mylines)
	for cnp in mylines: 
		plot_seghist(cnp, color)

def plot_seghist(cnp, color, plotbp=False):
	cnvals=np.array(cnp.cnvals)
	pvals=np.array(cnp.prevals)
	if cnp.start != cnp.end: 
		pvals=pvals #[cnvals!=-10]
		cnvals=cnvals * -1
		cnvals=np.cumsum(cnvals) #[np.not_equal(cnvals, -10)])
		lwd=cnp.lscore*20
		plt.step(pvals, cnvals, where='post', color=color, linewidth=lwd, alpha=0.5)
	elif plotbp: 
		plt.scatter(pvals, np.zeros(pvals.shape), color=color, marker='*', alpha=0.5)

if __name__=='__main__': 
	parser=argparse.ArgumentParser(description='plot seghists')
	parser.add_argument('seghists', help='a seghists.txt file')
	parser.add_argument('pdf', help='The pdf to save the output to')
	parser.add_argument('--genelist', help='a list of genes to be labeled')
	parser.add_argument('--plotall', action='store_true', default=False, help='Whether to plot all the seghists or just the ones with the genes of interest.')
	args=parser.parse_args()
	seglines=open(args.seghists, 'r').readlines()
	genes=['EGFR', 'CDKN2A', 'MDM2', 'PTEN', 'PDGFRA']
	fig=plt.figure(figsize=(8,8))
	plot_seglines_by_label(seglines, genes, args.plotall)
	plt.savefig(args.pdf)


