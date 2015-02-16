# Module for history segments output from cnavg pipeline
import sys, os
import re
import numpy as np
import math
import copy 
import cnavgpost.mergehistories.event_cycles_module as ecycles

class SegmentHistory:
	def __init__(self, myline):
		if isinstance(myline, str):
			data=myline.strip().split('\t')
			self.chr=data[0]
			self.start=int(data[1])
			self.end=int(data[2])
			self.chr2=int(data[3])
			self.cnvals=map(round, map(float, data[4].split(',')), 2)
			self.prevals=map(float, data[5].split(','))
			self.orders=map(float, data[6].split(','))
			self.likelihoods=map(float, data[7].split(,))
			self.numhists=map(float, data[8].split(,))
			self.prevalsd=map(float, data[9].split(','))
			self.ordersd=map(float, data[10].split(','))
		elif isinstance(myline, ecycles.Event):
			myevent=myline
			self.chr=myevent.chr
			self.start=myevent.start
			self.end=myevent.end
			self.chr2=myevent.chr2
			self.cnvals=[myevent.cnval]
			self.prevals=[mean(myevent.prevals)]
			self.prevalsd=[mean(myevent.prevalsd)]
			self.orders=[mean(myevent.orders)]
			self.ordersd=[mean(myevent.ordersd)]
			self.likelihoods=[myevent.likelihood]
			self.numhists=[myevent.numhists]

	def __str__(self):
		mystr="%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
			self.chr, self.start, self.end, self.chr2, 
			",".join(map(str, self.cnvals)), 
			",".join(map(str, self.prevals)), 
			",".join(map(str, self.orders)), 
			",".join(map(str, self.likelihoods)), 
			",".join(map(str, self.numhists)), 
			",".join(map(str, self.prevalsd)), 
			",".join(map(str, self.ordersd)))
		return mystr
		
	def __eq__(self, other): 
		return (self.samelocus(other)  and (self.cnvals == other.cnvals)) 
	
	def samelocus(self, other):
		return ((self.chr == other.chr) and (self.start == other.start) and (self.end==other.end))
	
	def overlaps(self, other): 
		return (self.chr == other.chr and self.start <= other.end and self.end >= other.start) 

	def comes_before(self, other): 
		return ((self.chr < other.chr) or (self.chr == other.chr and self.end < other.start))
	
	def comes_after(self, other): 
		return ((self.chr > other.chr) or (self.chr == other.chr and self.start > other.end))
	
	def mergedwith(self, other): 
		merged=False
		if (self.samelocus(other) and self.histories == other.histories): 
			if len(self.cnvals) > len(other.cnvals):
				bigger=self
				smaller=other
			else: 
				bigger=other
				smaller=self
			cnvalsboth=[val for val in smaller.cnvals if val in bigger.cnvals]	
			prevalsboth=[val for val in smaller.prevals if val in bigger.prevals]
			if len(prevalsboth) == len(smaller.prevals) and len(cnvalsboth) == len(smaller.cnvals): 
				merged=True
		return merged	

def merge_in_edge(edge, seghist): 
	new_segments=[]
	if seghist.overlaps(edge): 
		if seghist.start > edge.start: 
			newseg=SegmentHistory(edge)
			newseg.end=seghist.start
			new_segments.append(newseg)
		elif seghist.start < edge.start: 
			newseg=copy.deepcopy(seghist)
			newseg.end=edge.start-1
			new_segments.append(newseg)
			seghist.start=edge.start
		if seghist.end <= edge.end: 
			seghist.appendvals(edge)
			new_segments.append(seghist)
			newseg=SegmentHistory(edge)
			newseg.start=seghist.end+1
			new_segments.append(newseg)
		else: # seghist.end > edge.end: 
			newseg=copy.deepcopy(seghist)
			newseg.end=edge.end
			new_segments.append(newseg)
			seghist.appendvals(edge)
			seghist.start=edge.end+1
			new_segments.append(seghist)
	return(new_segments)
			
			
def split_and_merge_in(self, other):
	new_segments=[] 
	if self.overlaps(other):
		if other.start < self.start: 
			tmp=other
			other=self
			self=tmp
		if other.start > self.start: # as opposed to the starts being equal 
			seg1=copy.deepcopy(self)
			seg1.end=other.start-1
			new_segments.append(seg1)
		newstart=max(self.start, other.start)
		newend=min(self.end, other.end)
		seg2s=copy.deepcopy(self)
		(seg2s.start, seg2s.end)=(newstart, newend)
		seg2o=copy.deepcopy(other)
		(seg2o.start, seg2o.end)=(newstart, newend)
		if other.histories==self.histories: 
			plist = list(seg2s.prevals) + list(other.prevals)
			cnlist=list(seg2s.cnvals)+list(other.cnvals)
			(mergedprevals, mergedcnvals) = sort_unique_twolists(plist, cnlist)
#			sys.stderr.write("mergedprevals: %s, mergedcnvals: %s\n" % (str(mergedprevals), str(mergedcnvals)))
			seg2s.prevals=mergedprevals
			seg2s.prevalsarray=np.array(seg2s.prevals)
			seg2s.cnvals=mergedcnvals
			new_segments.append(seg2s)
		# if they are different histories, but have the same CN changes, then merge them
		elif self.cnvals == other.cnvals: 
			hlist=seg2s.histories + other.histories
			sortedi=sorted(range(len(hlist)), key=hlist.__getitem__)
			seg2s.histories=[hlist[i] for i in sortedi]
			seg2s.numhists=len(seg2s.histories)
			clist=seg2s.costs+seg2o.costs
			seg2s.costs=[clist[i] for i in sortedi]
			new_segments.append(seg2s)
			seg2s.prevalsarray = np.vstack([self.prevalsarray, other.prevalsarray])
		else: # they have different histories and different CN changes
			new_segments.append(seg2s)
			new_segments.append(seg2o)
		if other.end < self.end: 
			seg3=copy.deepcopy(self)
			seg3.start=other.end+1
			new_segments.append(seg3)
		elif other.end > self.end: 
			seg3=copy.deepcopy(other)
			seg3.start=self.end+1
			new_segments.append(seg3)
	else: 
		new_segments=(self, other)
	return new_segments

def sort_unique_twolists(plist, cnlist):
	sortedi = sorted(range(len(plist)), key=plist.__getitem__)
	plist=[plist[i] for i in sortedi]
	cnlist=[cnlist[i] for i in sortedi]
	newplist=[plist[0]]
	newcnlist=[cnlist[0]]
	for i in xrange(len(plist)-1): 
		if plist[i] != plist[i+1]: 
			newplist.append(plist[i+1])
			newcnlist.append(cnlist[i+1])
	return(newplist, newcnlist)

def create_min_distance_array(values, seghistories, histid_crossindex): 
	numhistories=len(histid_crossindex.keys())
	numevents=values.shape[1]
	mydistances=np.ones((numevents, numhistories), dtype=float) * 100 
	# find the minimum difference in CN values between two points for every history combination
	for i in xrange(0, numevents): 
		valA=values[:,i]
		histA=seghistories[i]
		for j in xrange(0, numevents): 
			valB=values[:,j]
			histB=seghistories[j]
			histBi=histid_crossindex[histB]
			dist=np.linalg.norm(valA-valB)
			if dist < mydistances[i][histBi]: 
				mydistances[i][histBi]=dist
	return(mydistances)

def create_similarity_matrix(dist_array, histids, histid_crossindex): 
	numhistories=len(histid_crossindex.keys())
	mymatrix=np.zeros((numhistories, numhistories), dtype=float)
	for histid in histid_crossindex.keys():
		rows=[i for i, x in enumerate(histids) if x==histid]
		newrow=histid_crossindex[histid]
		mymatrix[newrow,:]=dist_array[rows,:].mean(axis=0)
	return(mymatrix)

def score_similarity_matrix(similarity_array):
	numscores=similarity_array.shape[0]
	scores=np.zeros(numscores, dtype=float)
	if numscores ==1:
		scores[0]=np.NAN
	else:	
		for i in xrange(0, similarity_array.shape[0]):
			a=np.ma.array(similarity_array, mask=True)
			a.mask[i,:]=False
			a.mask[:,i]=False
			a.mask[i,i]=True
			scores[i]=np.mean(a)
	return(scores)
