#!/usr/bin/env python 

import argparse
import cnavgpost.genehistory.bedFileModule as bedmod
import cnavgpost.genehistory.segments_history_module as sgh 


def get_overlapping_seghist(bed, mysegs): 
	seghists=get_overlapping_seghists(bed, mysegs): 
	

def main(segfn, bedfn, pval=0.5): 
	mysegs=sgh.read_in_segs(segfn)
	currseg=mysegs[0]
	for l in open(bedfn, 'r'): 
		bed=BedEntry(l)
		seghist=get_overlapping_seghist(bed, mysegs)
		CNprofs=seghist.CNprofiles
		ltot=get_total_likelihood(CNprofiles)
		if ltot > pval: 
			bestCN=sorted(CNprofiles, key=lambda x: x.likelihood, reverse=True)[0] 
			mypreval=bestCN.pvals[0]
			myprevalsd=bestCN.pvalsd[0]
			mylscore=bestCN.likelihood
		else: 
			(mypreval, myprevalsd, mylscore)=(0,0, 0)
		sys.stdout.write("\t".join(map(str, [
			bed.chr, bed.start, bed.end, bed.name,
			ltot, mylscore, mypreval, myprevalsd])))


if __name__== "__main__": 
	parser=argparse.ArgumentParser(description='give seghists and a bedfile of genes, it will score each gene by the time that it gets first hit.')
	parser.add_argument('seghists', help='a file of seghists')
	parser.add_argument('bedfile', help='a bedfile of genes to score.')
	args=parser.parse_args()
	main(args.seghists, args.bedfile)
