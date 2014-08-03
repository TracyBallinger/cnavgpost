# Module for history segments output from cnavg pipeline
import sys, os
import re
import numpy as np
import math
import copy 
from braney_lines_module import *

class History: 
	for 

class SegmentHistory:
	def __init__(self, myline):
		if isinstance(myline, str):
			data=myline.strip().split('\t')
			self.chr=data[0]
			self.start=int(data[1])
			self.end=int(data[2])
			if re.search(',', data[3]):
				self.cnvals=map(float, data[3].split(','))
				self.prevals=map(float, data[4].split(','))
			else: 
				self.cnvals=[float(data[3])]
				self.prevals=[float(data[4])]
			self.Histories=[]
			self.cnvals=map(round, self.cnvals, 2)
			self.likelihood=float(data[5])
			self.numhists=int(data[6])
			self.costs=[0] * self.numhists
			self.histories=range(int(self.numhists))
		elif isinstance(myline, ):
			self.chr=myline.chr
			self.start=myline.start
			self.end=myline.end
			self.prevals=[myline.preval]
			self.cnvals=[round(myline.cnval/myline.preval, 2)]
			self.likelihood=0
			self.costs=[myline.complexity]
			self.numhists=1
			self.histories=[myline.historyid]
		self.numsims=1
		self.prevalsarray=np.array(self.prevals)

	def samelocus(self, other):
		return ((self.chr == other.chr) and (self.start == other.start) and (self.end==other.end))
	
	def __eq__(self, other): 
		return (self.samelocus(other)  and (self.cnvals == other.cnvals)) 
	
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
	
	def overlaps(self, other): 
		return (self.chr == other.chr and self.start <= other.end and self.end >= other.start) 

	def comes_before(self, other): 
		return ((self.chr < other.chr) or (self.chr == other.chr and self.end < other.start))
	
	def compute_prevalence_stats(self):
		w=np.exp(-1*np.array(self.costs))
#		sys.stderr.write("%s\t%d\t%d\t%s, shape: %s\nweights: %s, shape:%s\n" % (self.chr, self.start, self.end, str(self.prevalsarray), str(self.prevalsarray.shape), str(w), str(w.shape)))
		if self.prevalsarray.ndim==2: 
			waverage=np.average(self.prevalsarray, axis=0, weights=w)
			var=np.average((self.prevalsarray -waverage) ** 2, axis=0, weights=w)
		else: 
			waverage=self.prevalsarray	
			var=np.zeros(self.prevalsarray.shape)
		self.prevals=waverage
#		sys.stderr.write("self.prevals is %s\n" % (str(self.prevals)))
		prevalsd=np.sqrt(var)
		return(prevalsd)

	def get_numsims(self): 
		binwidth=10000
		numsims=10
		hranges=[0]*numsims
		for h in self.histories:
			bval=int(h/binwidth)
			hranges[bval]+=1
		vals=[]
		self.numsims=0
		for i in xrange(numsims):
			vals.append("%d:%d" % (i, hranges[i]))
			if hranges[i]>0: self.numsims+=1
		return ",".join(vals)

	def __str__(self):
		prevalsd= list(self.compute_prevalence_stats())
		simstring=self.get_numsims()
		mystr="%s\t%d\t%d\t%s\t%s\t%s\t%d\t%d\t%s\n" % (self.chr, self.start, self.end, ",".join(map(str, self.cnvals)), ",".join(map(str, self.prevals)), str(self.likelihood), self.numhists, self.numsims, ",".join(map(str, prevalsd)))
		return mystr

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
#				sys.stderr.write("mergedprevals: %s, mergedcnvals: %s\n" % (str(mergedprevals), str(mergedcnvals)))
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
