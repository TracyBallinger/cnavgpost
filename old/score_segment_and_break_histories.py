#!/inside/home/common/bin/python2.7
import sys, os
from numpy import array 
import copy, glob, subprocess, StringIO
import history_segments_module as histseg
bindir="/inside/home/tballing/cnavg-study/cnavgmbin"

def combine_simulations_to_unisegments(sampledir, history_number):
	eventlines=""
	braneyfiles=glob.glob(sampledir+"/"+"HISTORIES_?.braney")
	for hi in xrange(len(braneyfiles)): 
		braneyfile="%s/HISTORIES_%d.braney" % (sampledir, hi)
		sys.stderr.write("Issuing command: %s/get_history_tree_data.py %s --history_id=%d | awk '{print $0\"\t\"%d}\n'" % (bindir, braneyfile, history_number, hi))
		lines=subprocess.check_output(["%s/get_history_tree_data.py %s --history_id=%d | awk '{print $0\"\t\"%d}'" % (bindir, braneyfile, history_number, hi)], shell=True)
		eventlines+=lines
	tf=os.tmpfile()
	tf.write(eventlines)
	tf.seek(0)
	#tmpf=open("tmpevents", 'w')
	#tmpf.write(eventlines)
	#tmpf.close()
	out, err = subprocess.Popen(["sort -k1,1 -k2,2n -k3,3n -k8,8n"], stdin=tf, stdout=subprocess.PIPE, shell=True).communicate()
	mysegments=get_unisegments(StringIO.StringIO(out))
	for seg in mysegments:
		sys.stdout.write(str(seg) + "\n")

def get_unisegments(infile): 
	unisegments=[]
	line=infile.readline()
	current_seghist=histseg.Segment_history(line)
	working_seghists=[current_seghist]
	#Should put in a check somewhere that the input is properly ordered by chr location. 
	for line in infile: 
		myseg=histseg.Segment_history(line)
		working_tmphists=[]
		for workseg in working_seghists: 
			# Check if this segment overlaps the current one
			if (myseg.chr == workseg.chr) & (myseg.start <= workseg.end) & (myseg.end >= workseg.start): 
				# if this segment completely overlaps the current one, then add it. 	
				if (myseg.start <= workseg.start) & (myseg.end >= workseg.end): 
					workseg.addline(line)
					working_tmphists.append(workseg)
				# if this segment partially overlaps the current one, split the current one into or 2 or 3 non-overlapping ones
				else: 
					if (myseg.start > workseg.start):
						tmpseg=copy.deepcopy(workseg)
						tmpseg.end=myseg.start-1
						working_tmphists.append(tmpseg)
					else: # myseg.start <= workseg.start
						myseg.start=workseg.start
					tmpseg=copy.deepcopy(workseg)
					tmpseg.start=myseg.start
					tmpseg.addline(line)
					tmpseg.end = min(myseg.end, workseg.end)
					working_tmphists.append(tmpseg)
					if (myseg.end < workseg.end):
						tmpseg=copy.deepcopy(workseg)
						tmpseg.start=myseg.end+1
						working_tmphists.append(tmpseg)
			elif (workseg.chr == myseg.chr) & (myseg.end < workseg.start): 
				working_tmphists.append(workseg)
			elif (myseg.chr > workseg.chr) | (myseg.start > workseg.end):
			# print out the current working segment - we're done with it. 
				unisegments.append(workseg) 
			
		# after checking for overlap with all of the current working segments
		workseg=working_seghists[-1]
		if (myseg.chr == workseg.chr) & (myseg.start <= workseg.end) & (myseg.end > workseg.end): 
			myseg.start=workseg.end+1
			working_tmphists.append(myseg)
		elif (myseg.chr > workseg.chr) | (myseg.start > workseg.end):
			working_tmphists.append(myseg)
		working_seghists=working_tmphists

	# print out all of the rest of the segment at the end
	for seg in working_seghists: 
		unisegments.append(seg) 
	return unisegments


if __name__ == "__main__": 
	import argparse
	parser = argparse.ArgumentParser(description='Given a cn-avg output directory, it will merge the last step of multiple simulations and create file of universal segments with lists of CNvalues and prevalences across the different histories')
	parser.add_argument('sample', help='a cn-avg output directory.')
	args=parser.parse_args()
	history_number=2500
	combine_simulations_to_unisegments(args.sample, history_number)
	
