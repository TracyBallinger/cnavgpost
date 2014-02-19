#!/inside/home/common/bin/python2.7
import sys, os
import event_cycles_module as histseg

if __name__ == '__main__': 
	import argparse
	parser=argparse.ArgumentParser(description='Given a .evnts file, will generate a .gv (graphing file for graphviz')
	parser.add_argument('evnts', help='a .evnts file')

	args=parser.parse_args()
	eventf=open(args.evnts, 'r')
	mynodelabels={}
	nodecount=1
	sys.stdout.write("digraph G {\n")
	for eline in eventf: 
		myevent=histseg.Event(eline)
		myevent.make_segs_from_str()
		for seg in myevent.segs: 
			locus1="%s:%d" % (seg.chr, seg.start)
			locus2="%s:%d" % (seg.chr2, seg.end)
			if locus1 not in mynodelabels.keys(): 
				mynodelabels[locus1]=nodecount
				nodecount+=1
			if locus2 not in mynodelabels.keys(): 
				mynodelabels[locus2]=nodecount
				nodecount+=1
			if seg.adj: 
				style="dotted"
				arrowhead="none"
			else: 
				style="solid"
				arrowhead="normal"
			label=str(seg.cnval)
			if seg.cnval>0: color="red"
			elif seg.cnval<0: color="blue"
			else: color="black"
			sys.stdout.write("%s -> %s [style=%s, label=\"%s\", color=%s, arrowhead=%s, penwidth=3.0, shape=point]\n" % (str(mynodelabels[locus1]),  str(mynodelabels[locus2]), style, label, color, arrowhead))
	sys.stdout.write("}\n")	
			
