#!/Library/Frameworks/EPD64.framework/Versions/7.1/bin/python2.7

import sys, os
import argparse

parser = argparse.ArgumentParser(description='given two history data files, will match an entry from the second file with an entry from the first based on location overlap and other features.')
parser.add_argument('historyA', help='a History data file.')
parser.add_argument('historyB', help='a History data file.')
args=parser.parse_args()

# assume the input is sorted by genomic location 
Adata=open(args.historyA, 'r')
Bdata=open(args.historyB, 'r')

class Segment: 
	def __init__(self, myline): 
		data=myline.strip().split('\t')
		self.chr=data[0]
		self.start=int(data[1])
		self.end=int(data[2])
		self.cnval=float(data[3])
		self.order=int(data[4])
		self.preval=float(data[5])
		self.id=int(data[6])
		self.line=myline

lineB= Bdata.readline()
segb=Segment(lineB)
lineA= Adata.readline()
sega=Segment(lineA)
mycoord_distances=[]
myoverlaps=[]
myprevoverlaps=[]

while lineA != '':
	#sys.stderr.write("lineA: %s" % (lineA))
	#sys.stderr.write("lineB: %s" % (lineB))
	# if the two segments overlap, calculate the distance between them and save the overlapping entry
	if (lineB != '') and (segb.chr == sega.chr) and (segb.start < sega.end) and (sega.start < segb.end): 
		d=((sega.start-segb.start)**2 + (sega.end-segb.end)**2) ** (0.5)
		mycoord_distances.append(d)
		myoverlaps.append(segb)
		if (len(myprevoverlaps) > 0): 
			segb=myprevoverlaps.pop(0)
			lineB=segb.line
		else:  
			lineB=Bdata.readline()
			if lineB != '':
				segb=Segment(lineB)	
	# if segment B is past segment A, we've found all the overlaps (if any), so move on to the next segment in A
	elif (segb.chr == sega.chr and segb.start > sega.end) or (segb.chr > sega.chr) or (lineB ==''):
		#sys.stderr.write("distances: %s\n" % (str(mycoord_distances)))
		if myoverlaps:
			mini=[x for x in range(len(mycoord_distances)) if mycoord_distances[x] == min(mycoord_distances)]
			if len(mini)==1: 
				nearestb=myoverlaps[mini[0]]
			else: #if there's more than 1 closest segment, look at cn and prevalence to get the best match
				distcns=[]
				for i in mini:
					seg=myoverlaps[i]
					d2=(sega.cnval - seg.cnval) ** 2
					distcns.append(d2)
				d2i = distcns.index(min(distcns))
				#sys.stderr.write("distcns: %s, d2i: %d\n" % (str(distcns), d2i))
				nearestb=myoverlaps[mini[d2i]]
			sys.stdout.write("%s\t%d\t%d\t%f\t%d\t%f\n" % (lineA.strip(), nearestb.start, nearestb.end, nearestb.cnval, nearestb.order, nearestb.preval)) 
		else: 
			sys.stdout.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (lineA.strip(), "NA", "NA", "NA", "NA", "NA")) 
		lineA= Adata.readline()
		if lineA != '': 
			sega=Segment(lineA)
		myprevoverlaps=myoverlaps
		myprevoverlaps.append(segb)
		# backtrack in the B segments to the first previous overlap
		if (len(myprevoverlaps) > 0): 
			segb=myprevoverlaps.pop(0)
			lineB=segb.line
		mycoord_distances=[]
		myoverlaps=[]
	else: # while B segments are before A segment, 
		while ((segb.chr == sega.chr and segb.end < sega.start) or (segb.chr < sega.chr)) and lineB != '' : 
			if (len(myprevoverlaps) > 0): 
				segb=myprevoverlaps.pop(0)
				lineB=segb.line
			else:  
				lineB=Bdata.readline()
				if lineB != '': 
					segb=Segment(lineB)	


