#!/inside/home/common/bin/python
###  #!/inside/home/common/bin/python2.7
import sys, os
import numpy as np
from scipy.stats import f

def calculate_fstat(datafn): 
	data=np.genfromtxt(datafn, skip_header=1, dtype='S5,S5,i8,i8,i8,f8,f8,i8')
	explained=np.mean(data[:,5])
	N=len(data[:,0])
	sample_groups=np.unique(data[:,0])
	k=len(sample_groups)
	unexplained=np.arange(k)
	for i in xrange(k): 
		groupi=sample_groups[i]
		datai=data[np.where(data[:,0]==groupi)]
		unexplained[i]=np.mean(datai)
	fstat=explained/np.mean(unexplained)
	pvalue=f.cdf(fstat, N, k)
	return [fstat, k, N, pvalue]
		

if __name__ == '__main__': 
	import argparse
	parser=argparse.ArgumentParser(description='Calculates an F-statistic from an input datafile (the output of get_within_and_between_history_similarities.py')
	parser.add_argument('data', help='The input data file')
	args=parser.parse_args()
	(fstat, k, N, pvalue)=calculate_fstat(args.data)
	sys.stdout.write("fstat: %f\nk: %d\nN: %d\npvalue: %f" % (fstat, k, N, pvalue))
