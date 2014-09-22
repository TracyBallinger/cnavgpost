# Module for history segments output from cnavg pipeline
import sys, os
import re
import numpy as np
import math
import copy 
from cnavgpost.mergehistories.braney_lines_module import *
from cnavgpost.mergehistories.history_segments_module import Segment_history

class Combined_braneyseg(Segment_history):
	def __init__(self, myline):
		if isinstance(myline, str):
			myline=Braney_seg(myline)
		if isinstance(myline, Braney_seg):
			Segment_history.__init__(self, myline)
		self.seg=myline.seg
		self.adj=myline.adj
		self.chr2=myline.chr2
		self.st1=myline.st1
		self.st2=myline.st2
		
	def __str__(self):
		self.compute_prevalence_stats()
		self.get_numsims()
		if self.seg:
			return "%s\t%d\t%d\t%f\t%f\t%f\t%d\t%s\n" % (self.chr, self.start, self.end, self.cnvals[0], self.prevals[0], self.likelihood, self.numhists, self.numsims)
		else: 
			return "A\t%s\t%d\t%s\t%s\t%d\t%s\t%f\t%f\t%f\t%d\t%s\n" % (self.chr, self.start, self.st1, self.chr2, self.end, self.st2, self.cnvals[0], self.prevals[0], self.likelihood, self.numhists, self.numsims )
		
	def samelocus(self, other):
		if self.seg and other.seg: 
			return ((self.chr == other.chr) and (self.start == other.start) and (self.end==other.end) and (self.chr2==other.chr2) and (self.st1==other.st1) and (self.st2==other.st2))
		elif self.adj and other.adj: 
			return (((self.chr == other.chr) and (self.start == other.start) and (self.end==other.end) and (self.chr2==other.chr2) and (self.st1==other.st1) and (self.st2==other.st2)) or 
					((self.chr == other.chr2) and (self.chr2== other.chr) and (self.start == other.end) and (self.end==other.start) and (self.st1 == other.st2) and (self.st2 == other.st1)))

	def comes_before(self, other): 
		if self.seg and other.seg: 
			return ((self.chr < other.chr) or (self.chr == other.chr and self.end < other.start))
		elif self.adj and other.adj: 
			return ((self.chr2<other.chr) or 
				((self.chr==other.chr and self.chr2 == other.chr) and (self.end < other.start))) 

	def addin(self, other): 
		self.prevals = list(self.prevals) + list(other.prevals) 
		self.prevalsarray=np.array(self.prevals)
		self.costs= list(self.costs) + list(other.costs)
		self.histories= list(self.histories) + list(other.histories)
		self.numhists=len(self.histories)

