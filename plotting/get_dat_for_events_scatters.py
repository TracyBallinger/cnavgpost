#!/usr/bin/env python 

import argparse
import sys, os 
import cPickle as pickle
import cnavgpost.genehistory.bedFileModule as bed
import cnavgpost.mergehistories.event_cycles_module as ecycles

def get_chrs_lengths(event): 
	if len(event.segs)==0:
		event.make_segs_from_str()
	mychrs={}
	for s in event.segs:
		s.order_ends()
		if s.adj: 
			if s.chr != "None" and s.chr not in mychrs.keys(): 
				mychrs[s.chr]=0
			if s.chr2 != "None" and s.chr2 not in mychrs.keys(): 
				mychrs[s.chr2]=0
		else: 
			if s.chr in mychrs.keys(): 
				mychrs[s.chr]+= s.end-s.start+1
			else: 
				mychrs[s.chr]= s.end-s.start+1
	skeys=sorted(mychrs.keys())
	vals=[]
	for k in skeys: 
		vals.append(mychrs[k])
	return(",".join(skeys)+":"+",".join(map(str, vals)))

def write_event_dat(event): 
	mystr="\t".join(map(str, [
	event.id,
	event.likelihood,
	event.cnval,
	event.prevalmean, event.prevalsd,
	event.CharacterizeEvent()[0],  # this tells us the event type
	get_chrs_lengths(event),
	event.numhists
	])) + "\n"
	return mystr

def main(evntsfn, outbn, cutoff, bedfile): 
	events=pickle.load(open(evntsfn, 'rb'))
	datfh=open(outbn+".dat", 'w')
	myevents=[]
	for e in events: 
		if e.likelihood > cutoff:
			myevents.append(e) 
			datfh.write(write_event_dat(e))
	datfh.close()
	if bedfile: 
		mybeds=[]
		for line in open(bedfile, 'r').readlines(): 
			mybeds.append(bed.BedEntry(line))
		genedat=open(outbn+".gns", 'w')
		for e in myevents:
			geneset=[] 
			for bede in mybeds: 
				if e.check_overlap("chr"+bede.chr, bede.start, bede.end):
					geneset.append(bede.name)
			if len(geneset)==0: 
				genedat.write(e.id + "\tNone\n")
			else: 
				genedat.write(e.id + "\t" + ",".join(sorted(set(geneset)))+"\n")
		genedat.close()


if __name__=='__main__': 
	parser=argparse.ArgumentParser(description='create files needed for making event scatter plots. (The input to ...)')
	parser.add_argument('--cutoff', help='The pvalue cutoff to include an event', type=float, default=0)
	parser.add_argument('--bedfile', help='a bedfile of genes to label.  If this is included, a basename.gns file will be created.')
	parser.add_argument('evnts', help='a .pevnts file')
	parser.add_argument('outname', help='The basename of the output files.')
	args=parser.parse_args()
	main(args.evnts, args.outname, args.cutoff, args.bedfile)
