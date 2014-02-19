#!/inside/home/common/bin/python2.7 
import sys
import event_cycles_module as histseg


if __name__ == '__main__': 
	import argparse
	parser = argparse.ArgumentParser(description='Given a list of edges, will tell you how many of those cancel eachother out') 
	parser.add_argument('edgs')
	args=parser.parse_args()
	
	mysegs=[]
	for line in open(args.edgs): 
		edge=histseg.Event(line)
		edge.make_segs_from_str()
		mysegs += edge.segs
	sortedsegs=sorted(mysegs, key=lambda x: (x.chr, x.start, x.chr2, x.end, x.st1, x.st2))
	preseg=sortedsegs[0]
#	sys.stderr.write(str(preseg)+ "\n")
	numpairs = 0
	for seg in sortedsegs[1:]: 
#		sys.stderr.write(str(seg) + "\n")
		if seg.same_coords(preseg) and (seg.cnval == -1 * preseg.cnval):
			numpairs+=1
		preseg=seg 
	sys.stdout.write("pairs: %d\nsegments: %d\n" % (numpairs, len(sortedsegs)))

