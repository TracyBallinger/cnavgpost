# Module for merging histories created from the CN-AVG pipeline together
# Events are equivalent if they have identical coordinates across all segments and adjacencies and have the same CN value.  Basically, it has to be the same graph structure with the same Copy number change (The CN will be an integer).

import sys, os
import re, glob, subprocess
import itertools 
import math 
import copy
import cPickle as pickle 
from collections import Counter
import numpy as np
import pysam, gzip
from cnavgpost.mergehistories.braney_lines_module import * 

Global_BINWIDTH=10000
Global_MAXCOST=300 
Global_K=0
Global_EVENTTYPES=['any', 'amp', 'del', 'adj', 'oth']

# an Event is made up of multiple Braneysegs and some extra info
class Event:
	def __init__(self, bseglist):
		if isinstance(bseglist, str): 
			data=bseglist.strip().split('\t')
			self.id=data[0]  # The id is usually the Run id and the first iteration of the run that the event is found in.  It has the form [int].[int]
			self.likelihood=data[1]  # This is a float. 
			(self.numhists, self.numsims) = map(int, data[2].split(','))
			(self.prevalmean, self.prevalsd)=map(float, data[3].split(','))
			(self.ordermean, self.ordersd)=map(float, data[4].split(','))
			self.cnval=float(data[5])
			self.segstr=data[6]
			self.indyRunCounts=data[7]
			(self.numsegs, self.numadjs, self.numdisc) = map(int, [data[8]]+data[9].split(','))	
			self.histories=[] 
			self.histRanges=[]
			self.prevals=[] 
			self.orders=[] 
			self.segs=[]
			self.uppercosts={} 
			self.lowercosts={}
			self.dupsremoved=True
			
		elif isinstance(bseglist, list):
			self.id="%d.%d" % (bseglist[0].historyid, bseglist[0].order)
			self.likelihood=0
			(self.numhists, self.numsims)=(1,1)
			(self.prevalmean, self.prevalsd)=(bseglist[0].preval, "NA")
			(self.ordermean, self.ordersd)=(bseglist[0].order, "NA")
			self.cnval=round(abs(bseglist[0].cnval/bseglist[0].preval), 2)
			self.segstr=""
			(self.numsegs, self.numadjs, self.numdisc) = (0 , 0,0)	
			self.histories=[bseglist[0].historyid]
			self.histRanges=[]
			self.indyRunCounts=""
			self.prevals=[bseglist[0].preval]
			self.orders=[bseglist[0].order]
			self.segs=bseglist
			self.uppercosts=[bseglist[0].upperEventCost]
			self.lowercosts=[bseglist[0].lowerEventCost]
			self.dupsremoved=False
	
	def __str__(self):
		fstr="%s\t%s\t%d,%d\t%s,%s\t%s,%s\t%f\t%s\t%s\t%d\t%d,%d\n" % (self.id, str(self.likelihood), self.numhists, self.numsims, self.prevalmean, str(self.prevalsd), str(self.ordermean), str(self.ordersd), self.cnval, self.segstr, self.indyRunCounts, self.numsegs, self.numadjs, self.numdisc)
		return fstr
	
	def __eq__(self, other): 
		if self.segstr == "" : self.make_segstr()
		if other.segstr == "" : other.make_segstr()
		if self.segstr == other.segstr and self.cnval == other.cnval: 
			return True
		else:
			return False

	
	# update is done after equivalent events across multiple histories have been merged together.  Then the likelihood score and other stats for this event can be calculated. 	
	def update(self, historyScores):
		if self.segstr=="": 
			self.make_segstr()
		if self.numsegs==0: 
			self.count_discontsegs()
		self.sort_values_by_history()
		self.id = "%d.%d" % (self.histories[0], self.orders[0]) 
		self.numhists=len(self.histories)
		(self.numsims, self.indyRunCounts) = getIndRunCounts(self.histories)
		self.histRanges=getRanges(self.histories)
		if historyScores is not None: 
			self.compute_timing_wmeansd(historyScores)
	
	def addseg(self, bseg):
		self.segs.append(bseg)
	
	#Trimming of the event is done to make it more compact before it is written out to a pickled file. 
	def trim(self): 
		if self.segstr=="":
			self.make_segstr()
		self.segs=[]
		if not self.histRanges:
			self.histRanges=getRanges(self.histories)
			(self.numsims, self.indyRunCounts) = getIndRunCounts(self.histories)
		self.histories=[]

	# Unpacking of an event is done after it's been read in from a pickled file.  This basically undoes the trim (above). 
	def unpack(self):
		if not self.segs: 
			self.make_segs_from_str()
		if not self.histories: 	
			self.histories=listout_ranges(self.histRanges)	
	
	# This makes a string from the genomic coordinates of the event (event being a flow in the CN-AVG).  	
	def make_segstr(self):
		self.remove_dup_adj()
		self.merge_adjacent_segs()
		self.ordersegs()
		mystr=""
		mysegs=[]
		for seg in self.segs: 
			if seg.cnval < 0: sign="-"
			else: sign="+"
			if seg.adj: s="%s/%s:%d(%s)-%s:%d(%s)" % (sign, seg.chr, seg.start, seg.st1, seg.chr2, seg.end, seg.st2)	
			else: s="%s/%s:%d-%d" % (sign, seg.chr, seg.start, seg.end)
			mysegs.append(s)	
		mystr=",".join(mysegs)
		self.segstr=mystr
	
	# This orders the segments of the event, putting the lowest genomic coordinate first and then going around the rearrangement cycle. 
	def ordersegs(self): 
		firstseg=sorted(self.segs, key=lambda x: (x.chr, x.start, x.chr2, x.end, x.st1, x.st2))[0]
		sortedsegs=sorted(self.segs, key=lambda x: x.cycleorder)
		ifirst=sortedsegs.index(firstseg)
		second = sortedsegs[ifirst-1]
		if ((firstseg.chr2 == second.chr and firstseg.end == second.start) or (firstseg.chr2 == second.chr2 and firstseg.end == second.end)): 
			self.segs=sortedsegs[:(ifirst+1)][::-1]+sortedsegs[(ifirst+1):][::-1]
		else: 
			self.segs=sortedsegs[ifirst:]+sortedsegs[:ifirst]
	
	def merge_adjacent_segs(self):
		sortedsegs=sorted(self.segs, key=lambda x: x.cycleorder)
		l=len(sortedsegs)
		i=0
		while i < l:  
			seg = sortedsegs[i]
			if (seg.adj and (seg.chr == seg.chr2) and 
				(seg.start== seg.end+1 or seg.start==seg.end-1)):
				pseg=sortedsegs[i-1]
				if (i+1)==l: 
					nseg=sortedsegs[0]
				else: 
					nseg=sortedsegs[i+1]
				if nseg.seg and pseg.seg:
					minloc=min(nseg.start, nseg.end, pseg.start, pseg.end)
					maxloc=max(nseg.start, nseg.end, pseg.start, pseg.end)
					pseg.start=minloc
					pseg.end=maxloc
					sortedsegs.pop(i) #get rid of the adjacency
					if (i+1)==l:
						sortedsegs.pop(0) #get rid of the next segment
					else: 
						sortedsegs.pop(i)
					l=len(sortedsegs)
				else: 
					i+=1
			else: 
				i+=1
		self.segs=sortedsegs
	#	return sortedsegs

	def make_segs_from_str(self):
		mysegs=self.segstr.split(",")
		dummysegline="%s\t%s\t%s\t%f\t%f\t%d\t0\t0\t%d\t%d\t0\t0\t0\t0\t0\n"
		dummyadjline="A\t%s\t%s\t%s\t%s\t%s\t%s\t%f\t%f\t%d\t0\t0\t%d\t%d\t0\t0\t0\t0\t0\n"
		self.segs=[]
		cycleorder=0
		
		minhist=min(listout_ranges(self.histRanges))
		for loc in mysegs: 
			m=re.search('([+|-])/(\w+):(-?\d+)-(-?\d+)', loc)			
			if m: 
				coords=m.group(2,3,4)
				bseg=Braney_seg(dummysegline % (coords[0], coords[1], coords[2], self.cnval, self.prevalmean, minhist, cycleorder, self.ordermean))
				sign=m.group(1)
			else: 
				m=re.search('([+|-])/(\w+):(-?\d+)\(([+|-])\)-(\w+):(-?\d+)\(([+|-])\)', loc)			
				if m: 
					coords=m.group(2,3,4,5,6,7)
					bseg=Braney_seg(dummyadjline % (coords[0], coords[1], coords[2], coords[3], coords[4], coords[5], self.cnval, self.prevalmean, minhist, cycleorder, self.ordermean))
					sign=m.group(1)
			if sign=="+": bseg.cnval = bseg.cnval * -1
			self.segs.append(bseg)
			cycleorder+=1
	
	def add_Event_data(self, other):
		self.histories+= other.histories
		self.uppercosts+= other.uppercosts
		self.lowercosts += other.lowercosts
		self.prevals += other.prevals 
		self.orders += other.orders 
		
	def merge_Event_data(self, other): 
		indicesToAdd=get_index_of_non_intersecting_items(self.histories, other.histories)
		for i in indicesToAdd: 	
			self.histories.append(other.histories[i])
			self.uppercosts.append(other.uppercosts[i])
			self.lowercosts.append(other.lowercosts[i])
			self.prevals.append(other.prevals[i])
			self.orders.append(other.orders[i])
	
	def sort_values_by_history(self):
		histarray=np.array(self.histories) 
		myorderi = list(np.argsort(np.array(self.histories)))
		self.histories=[self.histories[i] for i in myorderi]
		self.uppercosts=[self.uppercosts[i] for i in myorderi]
		self.lowercosts=[self.lowercosts[i] for i in myorderi]
		self.prevals=[self.prevals[i] for i in myorderi]
		self.orders=[self.orders[i] for i in myorderi]

	def remove_dup_adj(self): 
		if not self.dupsremoved:
			nodup_segs=[]
			for seg in self.segs: 
				addseg=True
				if seg.adj: # if the braney_seg is an adjacency, see that it's not a duplicate of one that's been added already (there are two lines for every adjacency in .braney files)
					for adj in nodup_segs: 	
						addseg = (addseg and not seg.is_dup(adj))
					if addseg: 
						nodup_segs.append(seg)
				else: nodup_segs.append(seg)
			self.segs=nodup_segs
		self.dupsremoved=True

#Global_EVENTTYPES=['any', 'amp', 'del', 'adj', 'oth']
	def determineEventType(self): 
		mytype=0
		if len(self.segs) ==0: 
			self.make_segs_from_str()
		for seg in self.segs: 
			if seg.seg:
				# remember for segments, the sign is opposite
				if seg.cnval<0 and mytype != 2:
					mytype=2
				elif seg.cnval>0 and mytype != 1: 
					mytype=1
				else: 
					mytype=4
			else: 
				if mytype ==0:
					mytype=3
		return mytype
	
	def get_Event_length(self):
		seglen=0
		adjlen=0 
		for seg in self.segs: 
			if seg.seg: 
				seglen=max(seglen, seg.end-seg.start+1)
			else:
				if seg.chr==seg.chr2:  
					adjlen=max(adjlen, seg.end-seg.start)
		if seglen>0:
			len=seglen
		else: 
			len=adjlen	
		return len	

	def count_discontsegs(self): 
		if not self.dupsremoved: 
			self.remove_dup_adj() # this will order the segments and adjacencies by location
		myadjs=[]
		discont=0
		if len(self.segs) == 0: 
			self.make_segs_from_str()
		for seg in self.segs: 
			if seg.adj: 
				myadjs.append(seg)
				if seg.chr != seg.chr2: 
					discont+=1
		for i in xrange(len(myadjs)): 
			adj1=myadjs[i]
			for j in xrange(i,len(myadjs)):
				adj2=myadjs[j]
				if adj1.adjacency_cross(adj2): 
					discont+=1
		self.numsegs=len(self.segs)
		self.numadjs=len(myadjs)
		self.numdisc=discont 
	
	def check_overlap(self, chr, start, end):
		for seg in self.segs:
			if seg.seg: 
				if (seg.chr == chr and seg.start < end and seg.end>start): 
					return True
			elif seg.adj: 
				if seg.chr == chr and ((seg.start <= end and seg.start >= start) or (seg.end <= end and seg.end >= start)):
					return True 
		return False
	
	def multiline_str(self):
		mystr=""
		for seg in self.segs: 
			mystr+="%s\t%s\t%d\t%s\t%s\t%f\n" % (str(seg), self.id, len(self.histories), self.hranges, str(self.likelihood))
		return mystr
	
	def compute_timing_wmeansd(self, historyScores):
		hindices = historyids_to_indices(self.histories, historyScores)
		costsarray=historyScores[hindices, 1]
		w=np.exp(-1*Global_K*np.array(costsarray))
		pvals=np.array(self.prevals, dtype=float)
		waverage=np.average(pvals, weights=w)
		var=np.average((pvals - waverage)**2, weights=w)
		self.prevalmean=waverage
		self.prevalsd=math.sqrt(var)
		# do the same calculation for the ordering of the event within the history
		pvals=np.array(self.orders, dtype=float)
		waverage=np.average(pvals, weights=w)
		var=np.average((pvals - waverage)**2, weights=w)
		self.ordermean=waverage
		self.ordersd=math.sqrt(var)

def remove_signs_from_segstr(segstr):
	locs=segstr.split(',')
	newlocs=[]
	for loc in locs:
		m=re.search('([+|-])/(\w+):(-?\d+)\(([+|-])\)-(\w+):(-?\d+)\(([+|-])\)', loc)
		if m:
			(chr1, s, st1, chr2, e, st2)=m.group(2,3,4,5,6,7)
			newlocs.append("%s:%s(%s)-%s:%s(%s)" % m.group(2,3,4,5,6,7))
			sign=m.group(1)
		else:
			m=re.search('([+|-])/(\w+):(-?\d+)-(\-?\d+)', loc)
			if m:
				(chr1, s, e) = m.group(2,3,4)
				newlocs.append("%s:%s-%s" % m.group(2,3,4))
				sign=m.group(1)
		mysign=1
		if sign=="-":
			mysign=-1
	return (",".join(newlocs), mysign)

def sort_segs_in_cycle(seglist):
	segs=sorted(seglist, key=lambda x: x.cycleorder)
	for seg in segs: 
		sys.stderr.write("seg is: %s" % (str(seg)))
	currentseg=segs.pop(0)
	myorderedsegs=[currentseg]
	i=0
	maxiter=len(segs)**2
	tmpsegs=[]
	while (len(segs)>0 and maxiter >=0):
		seg=segs.pop(0)
		maxiter= maxiter -1 
		if (currentseg.chr2 == seg.chr2 and currentseg.end==seg.end): 
			seg.flip_ends()
			myorderedsegs.append(seg)
			currentseg=seg
		elif (currentseg.chr2 == seg.chr and currentseg.end==seg.start): 
			myorderedsegs.append(seg)
			currentseg=seg
		else: 
			segs.append(seg)
	if maxiter <0: 
		sys.stderr.write("len of segs: %d, and segstr: %s" % (len(segs), currentseg))
		for seg in myorderedsegs: 
			sys.stderr.write("orderedseg is: %s" % (str(seg)))
		for seg in segs: 
			sys.stderr.write("seg is: %s" % (str(seg)))
		sys.exit(-1)
	return myorderedsegs

def listout_ranges(ranges): 
	myilist=[]
	for (s,e) in ranges: 
		for i in range(s,e+1):
			myilist.append(i)
	return myilist

def ranges(ilist): 
	G=(list(x) for _,x in itertools.groupby(sorted(ilist), lambda x, c=itertools.count(): next(c)-x))
	return ",".join("-".join(map(str, (g[0], g[-1])[:len(g)])) for g in G)

def get_index_of_non_intersecting_items(list1, list2): 
	myis=[]
	for i in xrange(len(list2)): 
		if list2[i] not in list1: 
			myis.append(i)
	return myis 

def getRanges(vals):
	myranges=[] 
	vals.sort()
	rangestart=vals[0]
	rangeend=vals[0]
	for i in xrange(len(vals)-1): 
		if vals[i+1] > vals[i]+1: 
			rangeend=vals[i]
			myranges.append((rangestart, rangeend))
			rangestart=vals[i+1]
	rangeend=vals[-1]
	myranges.append((rangestart, rangeend))
	return myranges 
	
	
def getIndRunCounts(histories):
		binwidth=Global_BINWIDTH
		numbins=int(max(histories)/binwidth)+1
		hranges=[0]*numbins
		for h in histories: 
			bval=int(h/binwidth)
			hranges[bval]+=1
		vals=[]
		numsims=0
		for i in xrange(numbins): 
			vals.append("%d:%d" % (i, hranges[i]))
			if hranges[i]>0: numsims+=1
		hranges=",".join(vals)
		return (numsims, hranges)	

def get_cnvalue_from_seglist(bseglist):
	cnval=None
	for bseg in bseglist: 
		if bseg.seg: 
			cnval=round(bseg.cnval/bseg.preval, 2)
	if cnval==None: 
		cnval=round(abs(bseg.cnval/bseg.preval), 2)
	return cnval 

 
#### get_overlapping_events #######################
# Input: a chromosome, start and end, and the filename of a .evnts file 
# Output: a list of Events, where either an adjancency or a segment overlaps the region 
def get_overlapping_events(chr, start, end, evntsfn ):
	overlapping_events=[]
	for line in open(evntsfn, 'r'):
		sys.stderr.write("line: %s\n" % (line))
		myevent=Event(line)
		if (myevent.check_overlap(chr, start, end)):
			overlapping_events.append(myevent)

#### get_overlapping_events #######################
# Input: a chromosome, start and end, and the filename of a .braney file in tabix form
# Output: a list of Events, where either an adjancency or a segment overlaps the region 
def get_overlapping_events_tabix(chr, start, end, tabixfn, tabixfn2):
	eventid=""
	bsegs=[]
	eventlines=[]
	myevents=[] # a list of Events that contain a segment overlapping our genomic region
	tabixfile=pysam.Tabixfile(tabixfn)
	bseglines=[]
	try: 
		bseglines=tabixfile.fetch(reference=chr, start=start, end=end)
	except ValueError: 
		sys.stderr.write("Error in tabix fetch for %s, with %s:%d-%d, %s\n" % (tabixfn, chr, start, end, str(ValueError)))
	for bline in bseglines: 
		bsegs.append(Braney_seg(bline))
	badjlines=[]
	if (tabixfn2): 
		tabixfile2=pysam.Tabixfile(tabixfn2)
		try: 
			badjlines=tabixfile2.fetch(reference=chr, start=start, end=end)
		except ValueError: 
			sys.stderr.write("Error in tabix fetch for %s, %s\n" % (tabixfn2, str(ValueError)))
		for bline in badjlines: 
			bsegs.append(Braney_seg(bline))
	if len(bsegs)>0:
		bsegs.sort(key=lambda x: x.ptrid)
		for bseg in bsegs:
			if eventid == "": 
				eventid=bseg.ptrid
				eventlines.append(bseg) 
			elif bseg.ptrid == eventid: 
				eventlines.append(bseg)
			elif bseg.ptrid != eventid: 
				myevents.append(Event(eventlines))
				eventlines=[bseg]
				eventid=bseg.ptrid
		myevents.append(Event(eventlines))
	return myevents

# Input: a file handle to an eventsfile sorted by genomic description
# Output: a list of events where duplicate ones are merged together. 
def unique_c_events_sorted_file(evntsfile):
	unique_events=[]
	line1=evntsfile.readline()
	eventA=Event(line1)
	for line in evntsfile:
		eventB=Event(line)
		if eventA == eventB:
			eventA.add_Event_data(eventB) 
		else:
			eventA.update(None)
			unique_events.append(eventA)
			eventA=eventB
	unique_events.append(eventA)
	return unique_events

def get_split_offs(finalevents):
	splitoffs=[]
	for evnt in finalevents:
		tmpevents=split_off_duplicate(evnt)
		n=len(tmpevents)
		while n>0: 
			xlist=split_off_duplicate(evnt)
			tmpevents+=xlist
			n=len(xlist)
		splitoffs+=tmpevents
	return splitoffs 

def split_off_duplicate(event):
	dups=[h for h, k in Counter(event.histories).items() if k>1]
	if len(dups)>=1:
		newevent=copy.deepcopy(event)
		newevent.histories=[]
		newevent.prevals=[]
		newevent.orders=[]
		newevent.uppercosts=[]
		newevent.lowercosts=[]
		myprevals=[]
		for h in dups:
			histilist=[i for i, j in enumerate(event.histories) if j==h]
			pvals=[(event.prevals[i], i) for i in histilist]
			pvals.sort(key=lambda x: x[0])
			(p, i) = pvals[0]
			newevent.histories.append(event.histories.pop(i))
			newevent.prevals.append(event.prevals.pop(i))
			newevent.orders.append(event.orders.pop(i))
			newevent.uppercosts.append(event.uppercosts.pop(i))
			newevent.lowercosts.append(event.lowercosts.pop(i))
			if len(event.id.split('.'))>1:
				(idhist, order) = map(int, event.id.split('.'))
				newevent.id = "%d.%d" % (idhist, min(newevent.orders))
		return [newevent]
	else: 
		return []

def unique_c_events_sorted_list(eventslist): 
	unique_events=[]
	eventA=eventslist[0]
	for eventB in eventslist[1:]:
		if eventA == eventB:
			eventA.add_Event_data(eventB) 
		else: 
			unique_events.append(eventA)
			eventA=eventB
	unique_events.append(eventA)
	return unique_events
	
def make_events_from_braneyfn(braneyfn): 
	sys.stderr.write("working on %s\n" % braneyfn)
	braneyf=gzip.open(braneyfn, 'rb')
	eventorder=-1
	histid=0
	myevents=[]
	eventlines=[]
	batchbreak=1000 # after x histories, stop and merge them together.  
	breakcounter=0
	eventcount=0
	peventcount=0
	histcount=0
	for braneyline in braneyf: 
		if braneyline.strip() != '':
			braneyseg=Braney_seg(braneyline)
			if eventorder==-1: 
				eventorder=braneyseg.order #the order is unique within but not between histories
				histid=braneyseg.historyid
				eventlines.append(braneyseg)
			elif braneyseg.order == eventorder and braneyseg.historyid == histid:
				eventlines.append(braneyseg)
			elif braneyseg.order != eventorder or braneyseg.historyid != histid: 
				if histid != braneyseg.historyid: 
					histcount +=1
					breakcounter+=1
					peventcount=eventcount
					eventcount=0
				myevent=Event(eventlines)
				myevent.histories=[histcount]  #change this because sometimes when cn-avg.py is run in -c mode histories are added to the end of .braney files, but the history ids start from 0 again. 
				myevent.make_segstr()
				eventcount+=1
				myevents.append(myevent)
				if breakcounter >= batchbreak:
					sortedevents=sorted(myevents, key=lambda x: (x.segstr, x.cnval, x.prevals[0]))
					uniqueevents=unique_c_events_sorted_list(sortedevents)
			#		sys.stderr.write("%d\t%d\t%d\t%d\n" % (histid, peventcount, len(myevents), len(uniqueevents)))
					myevents=uniqueevents
					breakcounter=0
				eventlines=[braneyseg]
				eventorder=braneyseg.order
				histid=braneyseg.historyid
	# append the last event 
	myevent=Event(eventlines)
	myevent.make_segstr()
	myevents.append(myevent)
	sortedevents=sorted(myevents, key=lambda x: (x.segstr, x.cnval, x.prevals[0]))
	uniqueevents=unique_c_events_sorted_list(sortedevents)
	sys.stderr.write("num_events: %d, filtered: %d\n" % (len(myevents), len(uniqueevents)))
	return (uniqueevents)

def get_events_from_cnavgdir(cnavgdir, historyScores, totalp=0):
	allevents=[]
	braneyfiles=glob.glob(cnavgdir+"/"+"HISTORIES_*.braney")
	sys.stderr.write("braneyfiles: %s\n" % (str(braneyfiles)))
	for braneyfn in braneyfiles:
		sim=int(re.match(".*HISTORIES_(\d+)\.braney", braneyfn).group(1))
		events=make_events_from_braneyfn(braneyfn)
		for evnt in events:
			for i in xrange(len(evnt.histories)):
				evnt.histories[i] = evnt.histories[i] + sim*Global_BINWIDTH
			(idhist, order) = map(int, evnt.id.split('.'))
			evnt.id = "%d.%d" % (idhist+(sim*Global_BINWIDTH), order)
		sortedevents=sorted(allevents+events, key=lambda x:(x.segstr, x.cnval, x.prevals[0]))
		allevents=unique_c_events_sorted_list(sortedevents) 
	finalevents=allevents
	# check that the event doesn't need to be split again because it happens twice in some histories
	splitoffs=get_split_offs(finalevents)
	finalevents+=splitoffs
	if totalp==0:
		# get the total likelihood of all events for computing marginal likelihoods
		totalp=compute_likelihood_histories(historyScores[:,0], historyScores)
	for evnt in finalevents:
		evnt.update(historyScores)
		evnt.likelihood=compute_likelihood_histories(evnt.histories, historyScores, totalp)
		#evnt.trim()
	return finalevents

def merge_pevnts_files(pevntsfiles, outputfile, historyScores, totalp): 
	allevents=pickle.load(open(pevntsfiles[0], 'rb'))
	for pevntsfile in pevntsfiles[1:]:
		addinevents=pickle.load(open(pevntsfile, 'rb'))
		allevents=unique_c_events_sorted_list(sorted(allevents+addinevents, key=lambda x: (x.segstr, x.cnval, x.prevals[0])))
		sys.stderr.write("addinevents %d, allevents: %d\n" % (len(addinevents), len(allevents)))
	splitoffs=get_split_offs(allevents)
	sys.stderr.write("There are %d splitoffs\n" % len(splitoffs))
	allevents+=splitoffs
	for evnt in allevents:
		evnt.update(historyScores)
		evnt.likelihood=compute_likelihood_histories(evnt.histories, historyScores, totalp)
		evnt.trim()
	#return(allevents)
	pickle.dump(allevents, open(outputfile, 'wb'), pickle.HIGHEST_PROTOCOL)

def combine_history_statsfiles(cnavgdir): 
	statsfiles=glob.glob(cnavgdir+"/"+"HISTORY_STATS*")
	sys.stderr.write("statsfiles: %s\n" % (str(statsfiles)))
	mystats=np.array([])
	mysims=[]
	runlens=[]
	for statsfile in statsfiles:
		sim=int(re.match(".*HISTORY_STATS_(\d+)", statsfile).group(1))
		mysims.append(sim)
		print "sim is %d" % sim
		historystats=np.loadtxt(statsfile, dtype=int) 
		runlens.append(historystats.shape[0])
	runlen=max(runlens)
	for sim in mysims:
		statsfile=os.path.join(cnavgdir, "HISTORY_STATS_%d" % sim)
		historystats=np.loadtxt(statsfile, dtype=int, ndmin=2) 
		sys.stderr.write("dim of historystats is %s\n" % str(historystats.shape))
		if mystats.size==0: 
			mystats=np.zeros(((max(mysims)+1)*runlen, historystats.shape[1]+1), dtype=int)
		hids=np.array(range(historystats.shape[0]))+ sim*Global_BINWIDTH
		i = sim*runlen
<<<<<<< HEAD
		mystats[i:(i+historystats.shape[0]),:] = np.hstack((np.atleast_2d(hids).T, historystats))
=======
		mystats[i:i+historystats.shape[0],:] = np.hstack((np.atleast_2d(hids).T, historystats))
>>>>>>> 354718568d29786e105b3b1e631a1d0bdfb8a1a7
	return mystats

def get_historyScores(statsfile): 
	historystats=np.fromregex(statsfile, statsFileRegex, dtype=int)
	hids=np.array(range(historystats.shape[0]))
	historyScores=np.hstack((np.atleast_2d(hids).T, historystats))
	return historyScores

def compute_likelihood_histories(historyids, historyScores, denom=1):
	hindices = historyids_to_indices(historyids, historyScores)
	costsarray=np.mean(historyScores[hindices,1:3], axis=1)
	maskedcosts=np.ma.masked_where(costsarray==0, costsarray)
	x=np.sum(np.exp(-1*Global_K*maskedcosts))
	return x/float(denom)

def historyids_to_indices(historyids, historyScores): 
	hids=np.array(historyids, dtype=int)
	iter=np.fmod(hids, Global_BINWIDTH)
	sim=np.round(hids/Global_BINWIDTH)
	runlen=max(np.fmod(historyScores[:,0], Global_BINWIDTH))+1
	newi=iter+sim*runlen
	return newi

def get_breakpoints(events, historyid): 
	breaklocs={}
	for event in events:
		event.unpack() 
		istrue=0
		ispred=1
		if historyid != "": 
			if historyid in event.histories: 
				istrue=1
				if len(event.histories)==1:
					ispred=0
		for seg in event.segs: 
			locs=("%s\t%d" % (seg.chr, seg.start), "%s\t%d" % (seg.chr2, seg.end))
			for loc in locs: 
				if loc not in breaklocs.keys(): 
					breaklocs[loc]=[ispred, istrue]
				else: 
					x=breaklocs[loc]
					x[0] = x[0]+ ispred
					x[1]= x[1]+istrue
	return breaklocs

def get_event_costs_over_time(events, historyScores, run):
	runlen=max(np.fmod(historyScores[:,0], Global_BINWIDTH))+1
	myeventids=["history_cost"]
	hidmin=run*Global_BINWIDTH
	hidmax=hidmin+runlen
	tmpi=np.where((historyScores[:,0]>=hidmin) & (historyScores[:,0]<=hidmax))[0]
	costsarray=np.mean(historyScores[tmpi,1:3], axis=1)
	mydata=costsarray
	for event in events: 
		myeventids.append(event.id)
		event.histories=listout_ranges(event.histRanges)
		hids=np.array(event.histories)
		uppercosts=np.array(event.uppercosts)
		lowercosts=np.array(event.lowercosts)
		event_is=np.where((hids<=hidmax) & (hids>=hidmin))[0]
		eventcosts=np.ones(runlen)* -1
		for i in event_is:
			myi=np.fmod(hids[i], Global_BINWIDTH) 
			eventcosts[myi]=np.mean((event.uppercosts[i], event.lowercosts[i]))
		mydata=np.vstack((mydata, eventcosts))	
	return((myeventids, mydata))


def compute_lscore_over_time(events, historyScores, run, outputfn): 
	runlen=max(np.fmod(historyScores[:,0], Global_BINWIDTH))+1
	myrun=run
	myeventids=["total"]
	myruncounts=[1]
	hidmin=myrun*Global_BINWIDTH
	hidmax=hidmin+runlen
	tmpi=np.where((historyScores[:,0]>=hidmin) & (historyScores[:,0]<=hidmax))[0]
	costsarray=np.mean(historyScores[tmpi,1:3], axis=1)
	totalscores=np.cumsum(np.exp(-1*Global_K*costsarray))
	mydata=totalscores/totalscores[runlen-1]
	for event in events:
		event.histories=listout_ranges(event.histRanges)	
		myeventids.append(event.id)
		sys.stderr.write("working on event: %s\n" % (event.id)) 
		myruncounts.append(event.numsims)
		myhids=hids[np.where((hids<= hidmax) & (hids>=hidmin))]
		myhis = np.fmod(hids, Global_BINWIDTH)
		for i in xrange(myhis.size): 
			totalscore = totalscores[myhis[i]] 
			lscore = compute_likelihood_histories(myhids[:(i+1)], historyScores, totalscore)
			lscores[myhis[i]]=lscore
		for i in xrange(max(myhis),hidmax): 
			totalscore = totalscores[i] 
			lscore = compute_likelihood_histories(myhids, historyScores, totalscore)
			lscores[i]=lscore
		mydata = np.vstack((mydata, lscores))
	np.savetxt(outputfn, mydata.T, delimiter='\t', header="\t".join(myeventids) + "\n" + "\t".join(map(str, myruncounts)))

# This will merge events together how?... 	  	
def merge_events(events): 
	newevent=Event(events[0].segs)
	newevent.histories=events[0].histories
	newevent.uppercosts=events[0].uppercosts
	newevent.lowercosts=events[0].lowercosts
	newevent.orders=events[0].orders
	newevent.prevals=events[0].prevals
	for event in events[1:]:
		for i in xrange(len(event.histories)):
			historyid = event.histories[i]
			if historyid in newevent.histories: 
				i2=newevent.histories.index(historyid)
				if event.orders[i] < newevent.orders[i2]: 
					# ADD SPLITOFF EVENTS HERE?
					newevent.orders[i2]=event.orders[i]
					newevent.prevals[i2]=event.prevals[i]
			else: 
				newevent.histories.append(historyid)
				newevent.uppercosts.append(event.uppercosts[i])
				newevent.lowercosts.append(event.lowercosts[i])
				newevent.orders.append(event.orders[i])
				newevent.prevals.append(event.prevals[i])
	# SORT THE VALUES BY HISTORY ID
	newevent.id = "%d.%d" % (newevent.histories[0], newevent.orders[0]) 
	return newevent

def merge_events_by_type(events, historyScores=None): 
	 #order events by segstr
	sevents=sorted(events, key=lambda x: (x.segstr))
	unique_events=[]
	eventA=sevents[0]
	for eventB in sevents[1:]: 
		if eventB.segstr == eventA.segstr:  #don't need to look at CN because we only care about direction change, and this is in the segstr
			eventA.add_Event_data(eventB)
		else: 
			unique_events.append(eventA)
			eventA=eventB
	unique_events.append(eventA)
	finalevents=unique_events
	splitoffs=get_split_offs(finalevents)
	sys.stderr.write("There are %d splitoffs\n" % len(splitoffs))
	finalevents+=splitoffs
	for e in finalevents: 
		e.update(historyScores)
	return(finalevents)	
