#!/inside/home/common/bin/python2.7
import sys, os
import event_cycles_module as histseg
import pysam
import re

def merge_segments_into_regions(event): 
	mylocs=[]
	for seg in event.segs: 
		if seg.adj: 
			mylocs.append((seg.chr, seg.start, seg.start+1))
			mylocs.append((seg.chr2, seg.end, seg.end+1))
		else: 
			mylocs.append((seg.chr, seg.start, seg.end))
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
	parser=argparse.ArgumentParser(description='Given a bedfile with gene annotations and an .evnt file with rearrangment cycles, it will intersect the genes with the cycles to create groups of genes that are linked in rearrangements.')
	parser.add_argument('evnt', help='the .evnt file') 
	parser.add_argument('tabixfile', help='a tabix file of genes')
	args=parser.parse_args()
	mytabix=pysam.Tabixfile(args.tabixfile, 'r')
	eventf=open(args.evnt, 'r')
	for eline in eventf: 
		myevent=histseg.Event(eline)
		myevent.make_segs_from_str()
		coords= merge_segments_into_regions(myevent)	
		mygenes=[]
		for coord in coords: 
			region="%s:%d-%d" % (coord[0], coord[1], coord[2])
			try: 
				mygenes+=mytabix.fetch(region)  # make sure the chromosome format is the same as what is used in the tabix file 
			except ValueError:
				sys.stderr.write("Error in tabix fetch for %s, with %s, %s\n" % (args.tabixfile, region, str(ValueError))) 
		genesetstr="None"
		if len(mygenes) >0: 
			geneset=consolidate_genes(mygenes)
			genesetstr=",".join(map(str, geneset))
		sys.stdout.write(str(myevent).strip() + "\t" + genesetstr + "\n")

