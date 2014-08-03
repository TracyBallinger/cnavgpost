# Module for history segments output from cnavg pipeline

import sys, os
import re
import subprocess
import math

class Braney_seg:
	def __init__(self, myline):
		data=myline.strip().split('\t')
		if myline[0] != 'A': # check that it is not an adjacency line
			self.adj=False
			self.seg=True
			self.chr=data[0]
			self.start=int(data[1])
			self.st1=None
			self.chr2=self.chr
			self.end=int(data[2])
			self.st2=None
			self.cnval=-1 * float(data[3])
			self.preval=float(data[4])
			self.historyid=int(data[5])
			self.cycleorder=int(data[8])
			self.order=int(data[9])
			self.upperHistCost=int(data[10]) 
			self.lowerHistCost=int(data[11])
			self.upperEventCost=int(data[12]) 
			self.lowerEventCost=int(data[13]) 
			self.ptrid=data[14]
		else: # an adjacency line
			self.adj=True
			self.seg=False
			self.chr=data[1]
			self.start=int(data[2])
			self.st1=data[3]
			self.chr2=data[4]
			self.end=int(data[5])
			self.st2=data[6]
			self.cnval=-1 * float(data[7])
			self.preval=float(data[8])
			self.historyid=int(data[9])
			self.cycleorder=int(data[12])
			self.order=int(data[13])		
			self.upperHistCost=int(data[14]) 
			self.lowerHistCost=int(data[15]) 
			self.upperEventCost=int(data[16]) 
			self.lowerEventCost=int(data[17]) 
			self.ptrid=data[18]
			self.order=float(data[13])		
		self.order_ends()

	def __str__(self):
		if self.seg:
			return "%s\t%d\t%d\t%f\t%f\t%d\t0\t0\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n" % (self.chr, self.start, self.end, self.cnval, self.preval, self.historyid, self.cycleorder, self.order, self.upperHistCost, self.lowerHistCost, self.upperEventCost, self.lowerEventCost, self.ptrid)
		else: 
			return "A\t%s\t%d\t%s\t%s\t%d\t%s\t%f\t%f\t%d\t0\t0\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n" % (self.chr, self.start, self.st1, self.chr2, self.end, self.st2, self.cnval, self.preval, self.historyid, self.cycleorder, self.order, self.upperHistCost, self.lowerHistCost, self.upperEventCost, self.lowerEventCost, self.ptrid)

	def __eq__(self, other):
		vals_to_check=['adj', 'chr','chr2', 'start', 'end', 'st1', 'st2', 'cnval']
		self_dict=self.__dict__
		other_dict=other.__dict__
		equal=True
		for kval in vals_to_check: 
			if self_dict[kval]!= other_dict[kval]:
				equal=False
		return equal

	def same_coords(self, other): 
		vals_to_check=['adj', 'chr','chr2', 'start', 'end', 'st1', 'st2']
		self_dict=self.__dict__
		other_dict=other.__dict__
		equal=True
		for kval in vals_to_check: 
			if self_dict[kval]!= other_dict[kval]:
				equal=False
		return equal

	def is_dup(self, other): 
		if (
			((self.adj and other.adj) and (self.cnval == other.cnval) & (self.preval==other.preval)) 
			and (((self.chr == other.chr2) and (self.start == other.end) & (self.end == other.start) & (self.st1 == other.st2) and (self.st2== other.st1)) 
			or ((self.chr == other.chr) and (self.chr2==other.chr2) and (self.end==other.end) and (self.start == other.start) and (self.st1 == other.st1) and (self.st2 == other.st2)))): 
			return True
		else: 
			return False
	
	def order_ends(self):	
		if self.chr > self.chr2: 
			(x, y, z) = (self.chr, self.start, self.st1)
			(self.chr, self.start, self.st1) =(self.chr2, self.end, self.st2)
			(self.chr2, self.end, self.st2) =(x, y, z)
		elif (self.chr == self.chr2 and self.start > self.end): 
			(y,z)=(self.start, self.st1)
			(self.start, self.st1) =(self.end, self.st2)
			(self.end, self.st2)=(y,z)
	
	def flip_ends(self): 
		(x, y, z) = (self.chr, self.start, self.st1)
		(self.chr, self.start, self.st1) =(self.chr2, self.end, self.st2)
		(self.chr2, self.end, self.st2) =(x, y, z)
		

	def adjacency_cross(self, other): 
		if ((self.chr == self.chr2 == other.chr == other.chr2) and 
			(((self.start > other.start and self.start < other.end and self.end > other.end)) or 
			((self.start < other.start and self.end < other.end and self.end > other.start)))): 
			return True
	
	def connected_to(self, other): 
		if ((self.chr == other.chr and self.start == other.start) or  
			(self.chr == other.chr2 and self.start == other.end) or 
			(self.chr2 == other.chr and self.end == other.start) or 
			(self.chr2 == other.chr2 and self.end == other.end)): 
			return True
		else: 
			return False

def compute_likelihood(costs, totp): 
	mysum=0
	for c in costs: 
		mysum+=math.exp(-c)
	likelihood=mysum/totp
	return likelihood

def braney_to_seghist(braneyline): 
	bseg=Braney_seg(braneyline)
	segline="%s\t%d\t%d\t%f\t%f\t%d\t%d\n" % (bseg.chr, bseg.start, bseg.end, bseg.cnval/bseg.preval, bseg.preval, bseg.complexity, bseg.historyid)
	seghist=Segment_history(segline)
	return seghist

def make_tabix_from_braney(braneyfn, outdir):
	newfn="%s/%s.gz" % (outdir, os.path.basename(os.path.normpath(braneyfn)))
	if not os.path.isfile(newfn):
		subprocess.call(["gunzip -dc %s | grep -v ^A | sort -k1,1 -k2,2n | awk '$3>$2' | bgzip > %s" % (braneyfn, newfn)], shell=True)
	if not os.path.isfile(newfn+".tbi"):
		subprocess.call(["tabix -0 -s 1 -b 2 -e 3 %s" % (newfn)], shell=True)
	newfna="%s/%s-adj.gz" % (outdir, os.path.basename(os.path.normpath(braneyfn)))
	if not os.path.isfile(newfna):
		subprocess.call(["gunzip -dc %s | grep ^A | sort -k2,2 -k3,3n | bgzip > %s" % (braneyfn, newfna)], shell=True)
	if not os.path.isfile(newfna+".tbi"):
		subprocess.call(["tabix -0 -s 2 -b 3 -e 3 %s" % (newfna)], shell=True)
	return (newfn, newfna)

def get_total_history_prob(braneyfn):
	costs = subprocess.check_output(["gunzip -dc %s | grep -v ^A | grep -v '^$' | cut -f6,11 | uniq | cut -f2" % (braneyfn)], shell=True)
	mysum=0
	for cost in costs.strip().split("\n"):
		c=float(cost)
		mysum+= math.exp(-c)
	return mysum

