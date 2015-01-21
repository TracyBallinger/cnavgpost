#!/usr/bin/env python 
import numpy as np 
import argparse
import sys, os

class bedEntry: 
	def __init__(self, bedline):
		dat=bedline.strip().split("\t")
		self.chr=dat[0].replace("chr", "")  #get rid of leading chr if it's there
		self.start=int(dat[1])
		self.end=int(dat[2])
		self.name=dat[3]
		self.score=0
	def __str__(self): 
		return "\t".join(map(str, (self.chr, self.start, self.end, self.name, self.score)))

class bambamCNV: 
	def __init__(self, bbline):
		dat=bbline.strip().split("\t")
		self.chr=dat[0]
		self.start=int(dat[1])
		self.end=int(dat[2])
		self.ratio=float(dat[3])
	def __str__(self): 
		return "\t".join(map(str, (self.chr, self.start, self.end, self.ratio)))

def get_bambam_cnratio_for_genes(bedfn, cnvbblist, outfn): 
	cnvfiles=open(cnvbblist, 'r').readlines()
	headerlabels=["genename"]
	myscores=[]
	for line in cnvfiles: 
		cnvbbfn=line.strip()
		basefn=os.path.basename(cnvbbfn)
		headerlabels.append(basefn)
		(names, scores)=get_bambam_cnratio(bedfn, cnvbbfn)
		if myscores == []: 
			myscores=np.array(scores, ndmin=2)
		else: 
			myscores = np.vstack((myscores, np.array(scores)))
		sys.stderr.write("myscores dim: %s\n" % str(myscores.shape))
	rowlabels=np.array(names)[:, np.newaxis]
	sys.stderr.write("names dim %s\n" % str(rowlabels.shape))
	np.savetxt(outfn, np.hstack((rowlabels, myscores.T)), delimiter="\t", fmt='%s', header="\t".join(headerlabels))

def get_bambam_cnratio(bedfn, cnvbbfn): 
	bedentries=[]
	bedfh=open(bedfn, 'r')
	for line in bedfh: 
		bedentries.append(bedEntry(line))
	bedfh.close()
	bedentries=sorted(bedentries, key=lambda x: (x.chr, x.start, x.end))
	cnvbb=[]
	cnvbbfh=open(cnvbbfn, 'r')	
	for line in cnvbbfh: 
		cnvbb.append(bambamCNV(line))	
	cnvbbfh.close()	
	cnventries=sorted(cnvbb, key=lambda x: (x.chr, x.start, x.end))
	bbi=0
	bedi=0
	firstbbi=-1
	scores=[]
	names=[]
	weights=[]
	cnvals=[]
	while bedi < len(bedentries) and bbi < len(cnventries):
		bed=bedentries[bedi]
		overlap=get_overlap(bed, cnventries[bbi])
		if overlap>0: 
#			sys.stderr.write("overlapping\tgene[%d]\t%s\t%d\t%d\tcnv[%d]\t%s\t%d\t%d\toverlap: %d\n" % (bedi, bed.chr, bed.start, bed.end, bbi, cnventries[bbi].chr, cnventries[bbi].start, cnventries[bbi].end, overlap)) 
			weights.append(overlap)
			cnvals.append(cnventries[bbi].ratio)
			if firstbbi==-1: 
				firstbbi=bbi
			bbi=bbi+1
		elif overlap <0: #cnvbb entry comes before bed region
#			sys.stderr.write("comesbefore\tgene[%d]\t%s\t%d\t%d\tcnv[%d]\t%s\t%d\t%d\toverlap: %d\n" % (bedi, bed.chr, bed.start, bed.end, bbi, cnventries[bbi].chr, cnventries[bbi].start, cnventries[bbi].end, overlap)) 
			bbi=bbi+1
		else: 
#			sys.stderr.write("comesafter\tgene[%d]\t%s\t%d\t%d\tcnv[%d]\t%s\t%d\t%d\toverlap: %d\n" % (bedi, bed.chr, bed.start, bed.end, bbi, cnventries[bbi].chr, cnventries[bbi].start, cnventries[bbi].end, overlap)) 
			wts=np.array(weights, dtype=float)/float(bed.end-bed.start+1) 
			wa=np.ma.average(np.array(cnvals), weights=wts)
			bed.score=wa
			scores.append(wa)
			names.append(bed.name)
#			sys.stderr.write("done: %s\n" % str(bed)) 
#			sys.stderr.write("for %s, weights: %s, cnvals: %s, score: %f\n" % (bed.name, str(weights), str(cnvals), wa))
			bedi=bedi+1
			if firstbbi != -1: 
				bbi=firstbbi
			else: 
				bbi=max(bbi-1, 0)
			firstbbi=-1
			weights=[]
			cnvals=[]
	while bedi < len(bedentries): 
		bed=bedentries[bedi]
		scores.append(np.nan)
		names.append(bed.name)
#		sys.stderr.write("done: %s\n" % str(bed)) 
		bedi=bedi+1
	return (names, scores)	

def get_overlap(bed, cnv): 
	if (bed.chr == cnv.chr and
		bed.start < cnv.end and 
		bed.end > cnv.start): 
		overlap = (min(bed.end, cnv.end) - max(bed.start, cnv.start))
	elif (bed.chr < cnv.chr or 
		(bed.chr== cnv.chr and bed.end < cnv.start)):
		overlap = 0 
	else: 
		overlap = -1
	return overlap 



if __name__=='__main__': 
	parser=argparse.ArgumentParser(description='Given a bedfile and a list of bambam cnv files, it will generate a table of the average CN ratio for each gene per cnv file.')
	parser.add_argument('bedf', help="a bed file")
	parser.add_argument('cnvbb', help="a list of bambam cnv files.")
	parser.add_argument('outfile', help="The file to write the output to.")
	args=parser.parse_args()
	get_bambam_cnratio_for_genes(args.bedf, args.cnvbb, args.outfile)
	
