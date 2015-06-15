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
		(eventid, lscore, cnval, prevalmean, prevalsd, etype, chrms, numhists) = range(8)
		data=np.loadtxt(filename, dtype=str)
		self.ids=data[:,eventid]
		self.etype=data[:,etype]
		self.lscore=data[:,lscore].astype(float) 
		self.lscore=self.lscore/max(self.lscore)
		self.cnval=abs(data[:,cnval].astype(float))
		self.cnval[self.etype=='del']=self.cnval[self.etype=='del']*-1
		self.numhists=data[:,numhists].astype(int) 
		self.preval=data[:,prevalmean].astype(float) 
		self.prevalsd=data[:,prevalsd].astype(float)
		self.chrms=data[:,chrms]

def get_event_chr_lengths(edat): 
	xvals=[]
	lengths=[]
	mycolors=[]
	bottoms=[]
	chrmcolors=get_chrmcolors()
	for i in xrange(len(edat.chrms)): 
		(chrms, lens) = edat.chrms[i].strip().split(":")
		b=0
		for (c,l) in zip(chrms.split(','), lens.split(',')): 
			cnum=c.replace("chr", '')
			if cnum.lower() =="x" or cnum.lower()=="y":
				cnum=23
			cnum=int(cnum)
			mycolors.append(chrmcolors[cnum-1])
			xvals.append(edat.preval[i])
			l=int(l)
			if l==0: l=20
			lengths.append(l)
			bottoms.append(b)
			b+=l
	return (xvals, lengths, bottoms, mycolors)

def get_chrmlengths(clfile): 
	#xvals=(np.array(range(23), dtype='float')/200)*-1
	clens=range(23)
	for l in open(clfile, 'r'): 
		(c, l) = l.strip().split()
		c=c.replace("chr", "")
		if c.lower() =="x" or c.lower()=="y":
			c=23
		c=int(c)
		clens[c-1]=int(l)
	return(clens)
	#plt.bar(xvals, clens, color=mycols, alpha=0.5, width=0.01, linewidth=0)
		
def plot_event_bars(edat, clfile=""):
	ax=plt.gca()
	(barx, barh, barb, barc) = get_event_chr_lengths(edat)
	plt.bar(barx,barh, bottom=barb, color=barc, alpha=0.5, width=0.01, linewidth=0, zorder=3)
	plt.ylabel("length of event", fontsize=20)
	plt.xlabel("prevalence", fontsize=20)
	plt.xlim(1.1,0)
	plt.yscale('log')
	plt.ylim(1e4, 3e8)
	if clfile!="": 
		mycols=get_chrmcolors()
		chrmls=get_chrmlengths(clfile)
		chrlabels=map(lambda x: "chr"+str(x+1), range(len(chrmls)))
		plt.text(1.1, chrmls[0], chrlabels[0], ha='right') 
		plt.text(1.1, chrmls[-1], chrlabels[-1], ha='right') 
#		ax2=ax #ax.twinx()
#		locs, labels = plt.yticks()
#		ax2.set_yticks(list(locs)+[chrmls[0],chrmls[-1]])
#		ax2.set_yticklabels(list(labels)+[chrlabels[0], chrlabels[-1]])
		for i in xrange(len(chrmls)): 
			plt.axhline(y=chrmls[i], xmin=0, xmax=0.05, linewidth=3, color=mycols[i], alpha=0.5, zorder=1)


mytypecolors={
	'amp':'#e31a1c', 
	'del': '#1f78b4', 
	'amdl': '#6a3d9a', 
	'tran': '#33a02c', 
	'rtran': '#b2df8a', 
	'inv': '#ff7f00',
	'rinv': '#fdbf6f',
	'fuse': '#fb9a99', 
	'fiss': '#a6cee3', 
	'any':'#ffff99', 
	'oth': '#ffff99'
	}

# Add a couple of colors to the "Tableau 20" colors as RGB.  
def get_chrmcolors():
	# These are the "Tableau 20" colors as RGB.  
	tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199), 
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]
	tableau20.append((200, 82,0))
	tableau20.append((95, 158, 209))
	tableau20.append((65, 68, 81))
	#Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts. 
	for i in range(len(tableau20)):
		r, g, b = tableau20[i]
		tableau20[i] = (r / 255., g / 255., b / 255.)
	return tableau20

def plot_event_scatter(edat, annfn): 
	pointsize= edat.lscore*200 #genecounts
	etype=edat.etype
	myax=plt.gca()
	for t in np.unique(etype): 
		myi=np.where(etype==t)
		myax.scatter(edat.preval[myi], edat.cnval[myi], s=pointsize[myi], c=mytypecolors[t], zorder=2, label=t, alpha=0.5)
	plt.xlabel("prevalence", fontsize=20)
	plt.ylabel("copy number change", fontsize=20)
	plt.xlim(1.1, 0)
	plt.ylim(-5,20)
	handles, labels = myax.get_legend_handles_labels()
	sys.stderr.write("there are %d handles\n" % len(handles))
	leg=plt.legend(handles, labels)
	#leg.draw_frame(False)
	# add annotation to the graph
	if annfn != "":
		i=0
		for l in open(annfn, 'r'):
			(eventid, ann) = l.strip().split()
			if ann != "None": 
				plt.annotate(ann, xy=(edat.preval[i], edat.cnval[i]), fontsize=8, rotation=45, ha='left', va='bottom')
			i+=1

def plot_vcf_scatter(vcffile): 
	vcfdat=np.loadtxt(vcffile, usecols=(2,3,4))
	plt.scatter(vcfdat[:,2]*2, vcfdat[:,0]+vcfdat[:,1], alpha=0.5, linewidths=0, s=10)
	plt.ylim(0,200)
	plt.xlim(1.1, 0)
	plt.xlabel("Allelic fraction")
	plt.ylabel("Read Coverage")
	#ax2= ax1.twinx()
	#ax2.hist(vcfdat[:,2]*2, bins=100, histtype='step', color='gray', linewidth=2)
	#plt.xlim(1.1, 0)

def plot_vcf_histogram(vcffile): 
	vcfdat=np.loadtxt(vcffile, usecols=(2,3,4))
	plt.hist(vcfdat[:,2]*2, bins=100, histtype='step', color='gray', linewidth=2)
	plt.xlabel("Allelic fraction")
	plt.ylabel("Frequency")
	
def plot_sample_info_panel(datfile, outpdf, annfn, chrfile, vcffile): 
	dat=EventsData(datfile)
	fig=plt.figure(figsize=(8,11))
	if vcffile:
		gs=gridspec.GridSpec(4, 1, height_ratios=[1,1,1,0.5])
	else:
		gs=gridspec.GridSpec(2, 1, height_ratios=[1,1])
	ax1=plt.subplot(gs[1,0])
	plot_event_scatter(dat, annfn)
	# plot the length of each event
	ax2=plt.subplot(gs[0,0], sharex=ax1)
	plot_event_bars(dat, chrfile)
	plt.xlabel("")
	if vcffile: 
		ax3=plt.subplot(gs[2,0], sharex=ax1)
		plot_vcf_scatter(vcffile)
		ax4=plt.subplot(gs[3,0], sharex=ax1)
		plot_vcf_histogram(vcffile)
		plt.xlim(1.1,0)
	plt.savefig(outpdf) 
	
	
if __name__ == '__main__': 
	parser=argparse.ArgumentParser(description='Given a .evnts file and a .links file containing the links between the events, it will generate a graph of the history.') 
	parser.add_argument('dat', help='a events.dat dat file')
	parser.add_argument('pdf', help='The pdf of the scatter plot.')
	parser.add_argument('--gns', help='an events.gns file that contains annotation to put on the plots')
	parser.add_argument('--vcfdat', help='a vcf dat file for plotting allele variant frequencies')
	parser.add_argument('--chrlengths', help='a text file containing the lengths of the chromosomes.')
	args=parser.parse_args()
	plot_sample_info_panel(args.dat, args.pdf, args.gns, args.chrlengths, args.vcfdat)

