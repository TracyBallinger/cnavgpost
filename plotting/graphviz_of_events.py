#!/usr/bin/env python 
import cnavgpost.mergehistories.event_cycles_module as histseg 
import argparse
import sys

def graphviz_of_events(events, dotfn): 
	outf=open(dotfn, 'w')
	mychromcolors=get_chromcolors()
	outf.write("graph G{\n")
	outf.write("size = \"8.5,11!\";\n")
	arrowhead='normal'
	outf.write("node [shape=box, style=filled, fontsize=8, height=0.1];\n")
	for e in events: 
		e.make_segs_from_str()
		for s in e.segs: 
			if s.start <0: s.start=0
			if s.end <0: s.end=0
			chrom=s.chr
			c=chrom.replace("chr", "")
			if c=="X" or c=="Y": c=23
			if c=="None": c=0
			mycolor=mychromcolors[int(c)]
			label1="%s_%d" % (chrom, s.start)
			outf.write("%s [label=\"%s\", color=%s];\n" % (label1, label1, mycolor))
			if s.seg: 
				chrom=s.chr
				linestyle='bold'
			else:
				chrom=s.chr2 
				linestyle='dashed'
			c=chrom.replace("chr", "")
			if c=="X" or c=="Y": c=23
			if c=="None": c=0
			mycolor=mychromcolors[int(c)]
			label2="%s_%d" % (chrom, s.end)
			outf.write("%s [label=\"%s\", color=%s];\n" % (label2, label2, mycolor))
			if s.cnval<0: lcolor='blue'
			else: lcolor='red'
			outf.write("%s -- %s [style=%s, arrowhead=%s, color=%s, label=%s, fontsize=8]\n" % (label1, label2, linestyle, arrowhead, lcolor, e.id))	
	outf.write("}\n")


def get_chromcolors(): 
	mychromcolors=[
	'aquamarine', 
	'antiquewhite',
	'burlywood',
	'cadetblue',
	'chartreuse',
	'chocolate',
	'cornflowerblue', 
	'cyan', 
	'darkgoldenrod', 
	'darkorange', 
	'darksalmon', 
	'darkseagreen', 
	'deepskyblue',
	'deeppink1', 
	'firebrick', 
	'gold', 
	'forestgreen', 
	'darkorchid', 
	'coral', 
	'darkolivegreen', 
	'dodgerblue4',
	'darkturquoise',
	'crimson', 
	'firebrick1']	
	return mychromcolors

if __name__== '__main__': 
	parser = argparse.ArgumentParser(description = 'turns an events file into a graphviz file so you can see it easily.\nDo: $ neato -Tps file.dot > file.ps\n')
	parser.add_argument('evnts', help="an .evnts file")
	parser.add_argument('dotfn', help="The name of the output file (should end in .dot)")
	args=parser.parse_args()
	events=[]
	for l in open(args.evnts, 'r'): 
		events.append(histseg.Event(l))
	graphviz_of_events(events, args.dotfn)	
