#!/inside/home/common/bin/python2.7 

import sys, os
import argparse
import subprocess
import random
import glob

parser = argparse.ArgumentParser(description='Given a ')
parser.add_argument('mixes', help='a file with a list of the different mixes of different sequences samples and the mix proportions')
parser.add_argument('sample_key', help='A file with the format <sampleID><cnavg_output_directory>')
parser.add_argument('--prevalence', help='A cutoff value for the prevalence of an event.', type=float, default=0)
args=parser.parse_args()
bindir="/inside/home/tballing/cnavg-study/cnavgmbin"

def create_segment_files(line, mixcnt):
	data=line.strip().split('\t')
	samples=data[0].split(',')
	proportions=map(int, data[1].split(','))
	draws=int(data[2])
	tree_events=""
	history_ids=[] 
	# get the number of simulations per sample
	for i in xrange(0, len(samples)):
		sample=samples[i]
		numsims=proportions[i] 
		history_files=glob.glob(sampleID_directories[sample]+"/"+"HISTORY_TREES_?")
		sys.stderr.write("history_files: %s\n" % (str(history_files)))
		sys.stderr.write("number of history_files: %d, numsims: %d\n" % (len(history_files), numsims))
		# randomly pick numsims history files to add to the mix
		rsamp=random.sample(range(0, len(history_files)), numsims)
		for ri in rsamp:
			treefile="%s/HISTORY_TREES_%d" % (sampleID_directories[sample], ri) 
			braneyfile="%s/HISTORIES_%d.braney" % (sampleID_directories[sample], ri)
			sys.stderr.write("%s/get_history_tree_data.py %s %s --history_id=2500 | awk '{print $0\"\t\"%d.%d}'\n" % (bindir, treefile, braneyfile, i, ri)) 
			eventlines=subprocess.check_output(["%s/get_history_tree_data.py %s %s --history_id=2500 | awk '{print$0\"\t\"%d.%d}'" % (bindir, treefile, braneyfile, i, ri)], shell=True)
			tree_events+=eventlines
			history_ids.append("%d.%d" % (i, ri))
	eventlines=tree_events
	tf=os.tmpfile()
	tf.write(eventlines)
	tf.seek(0)
	out, err = subprocess.Popen(["sort -k1,1 -k2,2n -k3,3n -k8,8n | %s/get_segment_values_dist.py - " % (bindir)], stdin=tf, stdout=subprocess.PIPE, shell=True).communicate()
	segfile=open("mix%d_segments" % (mixcnt), 'w')
	segfile.write(out)
	segfile.close()
	tf.seek(0)
	out, err = subprocess.Popen(["cut -f8 | sort -u > mix%d_hids.txt" % (mixcnt)], stdin=tf, shell=True).communicate()

# read in the sample_keys: 
sampleID_directories={} # key: sampleID, value: cn-avg output directory
samplefile=open(args.sample_key, 'r')
for line in samplefile: 
	(sampleID, outputdir) = line.strip().split('\t')
	sampleID_directories[sampleID]=outputdir

# read in the mixing amounts
# a line in mixfile looks like: A,B	100,0	1
mixfile=open(args.mixes, 'r')
mymixes=[] # a list of the final processed files for each mix
mixcnt=1
for line in mixfile: 
	if line.startswith("#") or not line.strip(): # line is a comment or empty line
		continue
	if not (os.path.exists('mix%d_segments' % (mixcnt)) and os.path.exists('mix%d_hids.txt' % (mixcnt))):
		create_segment_files(line, mixcnt) 
	subprocess.check_call(["%s/compare_histories_by_segment.py --prevalence %f mix%d_segments mix%d_hids.txt > mix%d_hdistances.dat" % (bindir, args.prevalence, mixcnt, mixcnt, mixcnt)], shell=True)
	mixcnt+=1


