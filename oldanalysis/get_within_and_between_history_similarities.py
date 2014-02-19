#!/inside/home/common/bin/python2.7 
import sys, os
import subprocess
import glob, StringIO
import numpy as np
bindir="/inside/home/tballing/cnavg-study/cnavgmbin"
import create_segment_history_file
import compare_histories_by_segment

def calculate_history_distances(bfile1, bfile2):
	eventlines=""
	hi=1
	for braneyfile in (bfile1, bfile2): 
		lines=subprocess.check_output(["%s/get_history_tree_data.py %s --history_id=2500 | awk '{print $0\"\t\"%d}'" % (bindir, braneyfile, hi)], shell=True)
		eventlines+=lines
		hi+=1
	tf=os.tmpfile()
	tf.write(eventlines)
	tf.seek(0)
	out, err = subprocess.Popen(["sort -k1,1 -k2,2n -k3,3n -k8,8n | cut -f1-6,8"], stdin=tf, stdout=subprocess.PIPE, shell=True).communicate()
	segments=create_segment_history_file.get_unisegments(StringIO.StringIO(out))
	segfile=open("tmpsegments" ,'w')
	for segment in segments: 
		segfile.write(str(segment)+"\n")
	segfile.close()
	tf.seek(0)
	out, err = subprocess.Popen(["cut -f8 | sort -u > tmphids"], stdin=tf, shell=True).communicate()
	sys.stderr.write("%s/compare_histories_by_segment.py --prevalence %f tmpsegments tmphids\n" % (bindir, args.prevalence))
	outdat=compare_histories_by_segment.compare_histories_by_segment("tmpsegments", "tmphids", args.prevalence)	
	#out, err = subprocess.Popen(["%s/compare_histories_by_segment.py --prevalence %f tmpsegments tmphids" % (bindir, args.prevalence)], stdout=subprocess.PIPE, shell=True).communicate()
	# get stats on the distances between segment histories across all segments
	lines=StringIO.StringIO("".join(outdat))
	mydata=np.loadtxt(lines, usecols = (1,2,4))  # here we only need one column of distance because the distances should be the same since we are only comparing two histories. 
	lengths=mydata[:,1]-mydata[:,0]+1
	distances=mydata[:,2]
	numsegs=mydata.shape[0]-np.isnan(distances).sum()
	mdistances=np.ma.masked_array(distances, np.isnan(distances))
	stats=(np.median(lengths), lengths.min(), lengths.max(), mdistances.mean(), mdistances.std(), numsegs) 
	return stats

if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description='Given a list of samples, will calculate the difference in histories from simulations between samples and within samples.')
	parser.add_argument('sample_key', help='A file with the format <sampleID><cnavg_output_directory>')
	parser.add_argument('--prevalence', help='A cutoff value for the prevalence of an event.', type=float, default=0)
	parser.add_argument('--draws', help='The number of simulations draws to make from a sample or sample combinations.', type=int, default=10)
	args=parser.parse_args()
	# read in the sample_keys: 
	sampleID_directories={} # key: sampleID, value: cn-avg output directory
	samplefile=open(args.sample_key, 'r')
	sampleID_histories={} # key: sampleID, value: list of .braney files for that sample
	for line in samplefile: 
		(sampleID, outputdir) = line.strip().split('\t')
		sampleID_directories[sampleID]=outputdir
		sampleIDs=sampleID_directories.keys()
		braneyfiles=glob.glob(outputdir+"/HISTORIES_?.braney")
	sys.stdout.write("#sample1\tsample2\tlen_median\tlen_min\tlen_max\tsimilar_mean\tsimilar_sd\tN\n")

	for i in xrange(0, len(sampleIDs)):
		sample1=sampleIDs[i]
		braneyfiles1=glob.glob(sampleID_directories[sample1]+"/HISTORIES_?.braney")
		for j in xrange(i, len(sampleIDs)):
			sample2=sampleIDs[j]
			braneyfiles2=glob.glob(sampleID_directories[sample2]+"/HISTORIES_?.braney")
			for bfile1 in braneyfiles1:
				for bfile2 in braneyfiles2: 
					sys.stderr.write("working on %s, %s\n" % (bfile1, bfile2))
					stats=calculate_history_distances(bfile1, bfile2)
					sys.stdout.write("%s\t%s\t%s\n" % (sample1, sample2, "\t".join(map(str, stats))))


