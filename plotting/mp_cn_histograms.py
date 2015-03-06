#!/usr/bin/env python 

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

class CNhist(Figure): 
	def __init__(self, bbcnvtxt, numbins):
		self.numbins=numbins 
		self.bbfn=bbcnvtxt
		dat=np.loadtxt(self.bbfn, usecols=(1,2,3))
		(x, y) = getHistFromCounts(dat[:,1]-dat[:,0], dat[:,2], self.numbins) 
		Figure.__init__(self,
		

def getHistFromCounts(n, x, numbins):
	mycnts=np.zeros(numbins) 
	mids=np.linspace(min(x), max(x), num=numbins)
	d=(mids[1]-mids[0])/2
	breaks=np.zeros(numbins+1)
	breaks[0:numbins]=mids-d
	breaks[numbins]=mids+d
	for i in xrange(numbins): 
		mycnts[i] = sum(n[x>=breaks[i] & x<breaks[i+1]])
		


