#!/inside/home/common/bin/python2.7
import sys, os
import argparse
import gzip
import history_segments_module as histseg

def braney_to_segments(braneyfile, history_id): 
	braney = gzip.open(braneyfile, 'rb')
	prevevent=""
	numsegs=0
	numadj=0
	seglist=[]
	for braneyline in braney: 
		bseg=histseg.Braney_seg(braneyline)
		if bseg.historyid == history_id:
	#		sys.stdout.write("ptrid: %s, prevevent: %s, numsegs: %d, numadj: %d\n" % (bseg.ptrid, prevevent, numsegs, numadj)) 
		# if there are only adjacencies for an event, and no CN change for a segment, then it's an LOH.  Create a segment to indicate this. 
			if bseg.ptrid == prevevent: 
				seglist.append(bseg)
				if bseg.seg: numsegs+=1
				else: numadj+=1
			else: 
				if numsegs == 0 and numadj==2:
					lohadj=seglist[0]
					if (lohadj.chr == lohadj.chr2) or lohadj.chr =="None" or lohadj.chr2 == "None": 
						s=min(lohadj.start, lohadj.end)
						e=max(lohadj.start, lohadj.end)
						lohadj.start=s
						lohadj.end=e
						lohadj.seg=True
						lohadj.adj=False
						if lohadj.chr == "None": lohadj.chr=lohadj.chr2
						lohadj.chr2="LOH"
				for seg in seglist: 
					if seg.seg: 
						sys.stdout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (seg.chr, seg.start, seg.end, seg.cnval, seg.historyid, seg.preval, seg.ptrid))
				seglist=[bseg]
				if bseg.seg: (numsegs, numadj)=(1,0)
				else: (numsegs, numadj)=(0,1) 
				prevevent=bseg.ptrid
	for seg in seglist: 
		if seg.seg: 
			sys.stdout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (seg.chr, seg.start, seg.end, seg.cnval, seg.historyid, seg.preval, seg.ptrid))
	braney.close()

if __name__ == "__main__": 
	parser = argparse.ArgumentParser(description='given a HISTORY.braney file, will get the order and prevalence information for segments.') 
	parser.add_argument('braney', help='a HISTORIES.braney file')
	parser.add_argument('--history_id', help='which history to output', type=int, default=2500)
	args=parser.parse_args()
	braney_to_segments(args.braney, args.history_id)
