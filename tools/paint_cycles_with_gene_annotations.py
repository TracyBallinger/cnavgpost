#!/inside/home/common/bin/python2.7

import sys, os
import history_segments_module as histseg
import pysam
import re

def merge_segments_into_regions(segstr): 
	mylocs=[]
	segs=segstr.split(",")
	for loc in segs: 
		dat=re.split(":-|--|:|-|", loc)
		if len(dat) == 3: 
			mylocs.append([dat[0], int(dat[1]), int(dat[2])])
		else:
			mylocs.append([dat[0], int(dat[1]), int(dat[1])+1])
			mylocs.append([dat[2], int(dat[3]), int(dat[3])+1])
	mycoords=[]
	while len(mylocs) >0:
		(chr, start, end) = mylocs.pop(0)
		merger=True
		while merger:
			merger=False
			newlocs=[]
			for j in xrange(len(mylocs)): 
				(chr2, start2, end2)= mylocs[j]
				if (chr == chr2 and start <=end2 and end >= start2): 
					start=min(start, start2)	
					end=max(end, end2)	
					merger=True
				else:
					newlocs.append([chr2, start2, end2])
			mylocs=newlocs
		mycoords.append([chr, start, end])
	return mycoords

def consolidate_genes(genelist): 
	mynames=set([])
	for gene in genelist:
		(chr, start, end, name) = gene.split("\t")
		if name not in mynames: 
			mynames.add(name)
	return sorted(list(mynames))

if __name__ == '__main__': 
	import argparse
	parser=argparse.ArgumentParser(description='Given a tabix file with genomic annotations and an .evnt file with rearrangment cycles, it will replace the coordinates in the .evnts file with the annotation names.')
	parser.add_argument('evnt', help='the .evnt file')
	parser.add_argument('tabixfile', help='a tabix file of genes')
	args=parser.parse_args()
	mytabix=pysam.Tabixfile(args.tabixfile, 'r')
	eventf=open(args.evnt, 'r')
	for eline in eventf: 
		myevent=histseg.Event(eline)
		sys.stderr.write(eline)
		coords= merge_segments_into_regions(myevent.segstr)	
		mygenes=[]
		for coord in coords: 
			region="%s:%d-%d" % (coord[0], coord[1], coord[2])
			try: 
				mygenes+=mytabix.fetch(region)  # make sure the chromosome format is the same as what is used in the tabix file 
			except: 
				sys.stderr.write("could not fetch region: %s\n" % (region))	
		genesetstr="None"
		if len(mygenes) >0: 
			geneset=consolidate_genes(mygenes)
			genesetstr=",".join(map(str, geneset))
		sys.stdout.write(str(myevent).strip() + "\t" + genesetstr + "\n")

