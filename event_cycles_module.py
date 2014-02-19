# Module for history segments output from cnavg pipeline
import sys, os
import re, glob
import itertools 
import math 
import numpy as np
import pysam, gzip
from braney_lines_module import * 

Global_BINWIDTH=10000
Global_EVENTTYPES=['any', 'amp', 'del', 'adj']
# an Event is made up of multiple Braneysegs and some extra info
class Event:
	def __init__(self, bseglist):
		if isinstance(bseglist, str): 
			data=bseglist.strip().split('\t')
			self.segs=[]
			self.segstr=data[6]
			self.hranges="NA" 
			tmpx=map(int, data[2].split(','))
			self.numhists=tmpx[0]
			if len(tmpx)>1: self.numsims=tmpx[1]
			else: self.numsims=0
			self.histories=range(self.numhists)  
			self.cnval=float(data[5])
			(self.prevalmean, self.prevalsd)=map(float, data[3].split(','))
			(self.ordermean, self.ordersd)=map(float, data[4].split(','))
			self.prevals=[self.prevalmean]* self.numhists
			self.orders=[self.ordermean]* self.numhists
			self.id=data[0]
			self.costs=[0]*self.numhists
			self.likelihood=data[1]
			self.dupsremoved=True
		elif isinstance(bseglist, list):
			self.segs=bseglist
			self.segstr=""
			self.histories=[bseglist[0].historyid]
			self.hranges="NA"
			self.numhists=1
			self.numsims=1
			#self.cnval=get_cnvalue_from_seglist(bseglist) 
			self.cnval=round(abs(bseglist[0].cnval/bseglist[0].preval), 2)
			self.prevals=[bseglist[0].preval]
			self.prevalmean=bseglist[0].preval
			self.prevalsd="NA"
			self.orders=[bseglist[0].order]
			self.ordermean=bseglist[0].order
			self.ordersd="NA"
			self.id="%d.%d" % (self.histories[0], self.orders[0]) #bseglist[0].ptrid
			self.costs=[bseglist[0].complexity]
			self.likelihood=0
			self.dupsremoved=False
	
	def __str__(self):
		if self.segstr=="": 
			self.make_segstr()
		if self.hranges=="NA": 
			self.get_hranges()
		#if self.id == "NA":
		#	self.make_id()
		#self.hranges=ranges(self.histories)
		if len(self.segs) == 0: 
			self.make_segs_from_str()
		(numadjs, numdisc) = self.count_discontsegs()
		numsegs=len(self.segs)
		fstr="%s\t%s\t%d,%d\t%s,%s\t%s,%s\t%f\t%s\t%s\t%d\t%d,%d\n" % (self.id, str(self.likelihood), self.numhists, self.numsims, self.prevalmean, str(self.prevalsd), str(self.ordermean), str(self.ordersd), self.cnval, self.segstr, self.hranges, numsegs, numadjs, numdisc)
		return fstr
	
	def __eq__(self, other): 
		if self.segstr == "" : self.make_segstr()
		if other.segstr == "" : other.make_segstr()
		if self.segstr == other.segstr and self.cnval == other.cnval: 
			return True
		else:
			return False
	
	def addseg(self, bseg):
		self.segs.append(bseg)
	
	def trim(self): 
		if self.segstr=="":
			self.make_segstr()
		self.segs=[]
		
#	def make_id(self): 
#		firsthistory=min(self.histories)
#		firstorder=min(x for x in self.orders if self.histories==firsthistory)
#		self.id = "%d.%d" % (firthistory, firstorder)

	def make_segstr(self):
		self.remove_dup_adj()
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
	
	def make_segs_from_str(self):
		mysegs=self.segstr.split(",")
		dummysegline="%s\t%s\t%s\t%f\t%f\t0\t0\t0\t0\t0\t0\t0\n"
		dummyadjline="A\t%s\t%s\t%s\t%s\t%s\t%s\t%f\t%f\t0\t0\t0\t0\t0\t0\t0\n"
		for loc in mysegs: 
			m=re.search('([+|-])/(\w+):(-?\d+)-(-?\d+)', loc)			
			if m: 
				coords=m.group(2,3,4)
				bseg=Braney_seg(dummysegline % (coords[0], coords[1], coords[2], self.cnval, self.prevalmean))
				sign=m.group(1)
			else: 
				m=re.search('([+|-])/(\w+):(-?\d+)\(([+|-])\)-(\w+):(-?\d+)\(([+|-])\)', loc)			
				if m: 
					coords=m.group(2,3,4,5,6,7)
					bseg=Braney_seg(dummyadjline % (coords[0], coords[1], coords[2], coords[3], coords[4], coords[5], self.cnval, self.prevalmean))
					sign=m.group(1)
			if sign=="-": bseg.cnval = bseg.cnval * -1
			self.segs.append(bseg)
	
	def add_Event_data(self, other):
		self.histories += other.histories
		self.costs += other.costs
		self.prevals+= other.prevals
		self.orders+= other.orders
		self.numhists=len(self.histories)
		self.hranges="NA"

	def remove_dup_adj(self): 
		if not self.dupsremoved:
			nodup_segs=[]
			for seg in self.segs: 
				addseg=True
				if seg.adj: # if the braney_seg is an adjacency, see that it's not a duplicate of one that's been added already (there are two lines for every adjacency in .braney files)
					for adj in nodup_segs: 	
						addseg = (addseg and not seg.is_dup(adj))
					if addseg: 
					#	seg.order_ends()
						nodup_segs.append(seg)
				else: nodup_segs.append(seg)
			self.segs=sorted(nodup_segs, key=lambda x: (x.chr, x.start, x.chr2, x.end, x.st1, x.st2))
		self.dupsremoved=True
	
	def get_hranges(self):
		binwidth=Global_BINWIDTH
		numsims=int(max(self.histories)/binwidth)+1
		hranges=[0]*numsims
		for h in self.histories: 
			bval=int(h/binwidth)
			hranges[bval]+=1
		vals=[]
		self.numsims=0
		for i in xrange(numsims): 
			vals.append("%d:%d" % (i, hranges[i]))
			if hranges[i]>0: self.numsims+=1
		self.hranges=",".join(vals)	

#Global_EVENTTYPES=['any', 'amp', 'del', 'adj']
	def determineEventType(self): 
		mytype=0
		if len(self.segs) ==0: 
			self.make_segs_from_str()
		for seg in self.segs: 
			if seg.seg:
				if seg.cnval>0 and mytype != 2:
					mytype=1
				elif seg.cnval<0 and mytype != 1: 
					mytype=2
			else: 
				if mytype ==0:
					mytype=3
		return mytype
	
	def get_Event_length(self):
		len=0 
		for seg in self.segs: 
			if seg.seg: 
				len=max(len, seg.end-seg.start+1)
		return len	

	def count_discontsegs(self): 
		if not self.dupsremoved: 
			self.remove_dup_adj() # this will order the segments and adjacencies by location
		myadjs=[]
		discont=0
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
		return (len(myadjs), discont)		
	
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
	
	def compute_timing_wmeansd(self):
		w=np.exp(-1*np.array(self.costs))
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

def listout_ranges(ranges): 
	myilist=[]
	for x in ranges.split(','): 
		se=x.split('-')
		s=int(se[0])
		e=int(se[-1])
		for i in range(s,e+1):
			myilist.append(i)
	return myilist

def ranges(ilist): 
	G=(list(x) for _,x in itertools.groupby(sorted(ilist), lambda x, c=itertools.count(): next(c)-x))
	return ",".join("-".join(map(str, (g[0], g[-1])[:len(g)])) for g in G)

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

def merge_events(events): 
	newevent=Event(events[0].segs)
	newevent.histories=events[0].histories
	newevent.costs=events[0].costs
	newevent.orders=events[0].orders
	newevent.prevals=events[0].prevals
	newevent.id = events[0].id 
	for event in events[1:]:
		for i in xrange(len(event.histories)):
			historyid = event.histories[i]
			if historyid in newevent.histories: 
				i2=newevent.histories.index(historyid)
				if event.orders[i] < newevent.orders[i2]: 
					newevent.orders[i2]=event.orders[i]
					newevent.prevals[i2]=event.prevals[i]
			else: 
				newevent.histories.append(historyid)
				newevent.costs.append(event.costs[i])
				newevent.orders.append(event.orders[i])
				newevent.prevals.append(event.prevals[i])
		newevent.segs.append(event.segs)
		#for listname in ['histories', 'costs', 'orders', 'prevals']: #(histories, event.histories, costs, event.costs, orders, event.orders, prevals, event.prevals):
		#	eventlist = getattr(event, listname)
		#	addedx = [eventlist[i] for i in addedi]
		#	newlist = getattr(newevent, listname) 
		#	newlist += addedx
	return newevent

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
			unique_events.append(eventA)
			eventA=eventB
	unique_events.append(eventA)
	return unique_events

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
	braneyf=gzip.open(braneyfn, 'rb')
	eventid=""
	histid=0
	myevents=[]
	eventlines=[]
	for braneyline in braneyf: 
		if braneyline.strip() != '':
			braneyseg=Braney_seg(braneyline)
			if eventid=="": 
				eventid=braneyseg.ptrid  #the pointer id is unique within but not between histories
				histid=braneyseg.historyid
				eventlines.append(braneyseg)
			elif braneyseg.ptrid == eventid and braneyseg.historyid == histid:
				eventlines.append(braneyseg)
			elif braneyseg.ptrid != eventid or braneyseg.historyid != histid: 
				myevent=Event(eventlines)
				myevent.likelihood=compute_likelihood(myevent.costs, 1)
				s=str(myevent)
				myevents.append(myevent)
				eventlines=[braneyseg]
				eventid=braneyseg.ptrid
				histid=braneyseg.historyid
	# append the last event 
	myevent=Event(eventlines)
	myevent.likelihood=compute_likelihood(myevent.costs, 1)
	s=str(myevent)
	myevents.append(myevent)
	return (myevents)

def get_events_from_cnavgdir(cnavgdir):
	allevents=[]
	braneyfiles=glob.glob(cnavgdir+"/"+"HISTORIES_*.braney")
	sys.stderr.write("braneyfiles: %s\n" % (str(braneyfiles)))
	for braneyfn in braneyfiles:
		sim=int(re.match(".*HISTORIES_(\d+)\.braney", braneyfn).group(1))
		events=make_events_from_braneyfn(braneyfn)
		sortedevents=sorted(events, key=lambda x: (x.segstr, x.cnval))
		uniqueevents=unique_c_events_sorted_list(sortedevents)
		for evnt in uniqueevents:
			for i in xrange(evnt.numhists):
				evnt.histories[i] = evnt.histories[i] + sim*Global_BINWIDTH
			(idhist, order) = map(int, evnt.id.split('.'))
			evnt.id = "%d.%d" % (idhist+(sim*Global_BINWIDTH), order)
		allevents+=uniqueevents
		sys.stderr.write("num_events: %d, filtered: %d\n" % (len(events), len(allevents)))
    # do a final merge of events across different simulations
	sortedevents=sorted(allevents, key=lambda x: (x.segstr, x.cnval))
	sys.stderr.write("number of allevents: %d, sortedevents: %d\n" % (len(allevents), len(sortedevents)))
	finalevents=unique_c_events_sorted_list(sortedevents)
	totalp=get_total_likelihood(finalevents)
	for evnt in finalevents:
		evnt.compute_timing_wmeansd()
		evnt.likelihood=compute_likelihood(evnt.costs, totalp)
		evnt.trim()
	return finalevents

def get_total_likelihood(events): 
	allhistoryids=[]
	allcosts=[]
	for event in events: 
		for i in xrange(len(event.histories)): 
			hid = event.histories[i]
			if hid not in allhistoryids: 
				allhistoryids.append(hid)
				allcosts.append(event.costs[i])
	totalp = compute_likelihood(allcosts, 1)
	return totalp 
	
