#!/inside/home/common/bin/python2.7

import sys, os
import argparse
from numpy import array 

parser = argparse.ArgumentParser(description='given a file of events and nearest segments to the event (output from find_nearest_segment.py), will return the events with the average order and prevalance across alternate histories.')
parser.add_argument('history_pairs', help='a file of events and their nearest neighbor from alternate histories.')

args=parser.parse_args()
infile=open(args.history_pairs, 'r')

class Segment2: 
    def __init__(self, myline):
		data=myline.strip().split('\t')
		self.chr=data[0]
		self.start=int(data[1])
		self.end=int(data[2])
		self.cnval=float(data[3])
		self.order=int(data[4])
		self.preval=float(data[5])
		self.id=int(data[6])
		self.nnorders=[]
		self.nnpreval=[]
		self.refline="\t".join(data[0:7])

eventhash={}   #key is a ptr or event id, values are segments

for line in infile: 
	seg=Segment2(line)
	data=line.strip().split('\t')
	nnorder=data[10]
	nnpreval=data[11]
	if seg.id not in eventhash.keys(): 
		eventhash[seg.id]=seg
		refseg=seg
	else: 
		refseg=eventhash[seg.id]
	if (nnorder != "NA"):
		refseg.nnorders.append(nnorder)
		refseg.nnpreval.append(nnpreval)
	
for key in eventhash.keys(): 
	myseg=eventhash[key]
	myorders=myseg.nnorders
	myorders.append(myseg.order)
	myprevals=myseg.nnpreval
	myprevals.append(myseg.preval)
	myorders=array(map(float, myorders))
	myprevals=array(map(float, myprevals))
	sys.stdout.write("%s" % (myseg.refline))
	#sys.stdout.write("\t%f\t%f\t%f\t%f\t%d\n" % (myorders.mean(), myorders.std(), myprevals.mean(), myprevals.std(), len(myorders)))
	sys.stdout.write("\t%s\t%s\n" % (",".join(map(str, myseg.nnorders)), ",".join(map(str, myseg.nnpreval))))
	
	
