#!/inside/home/common/bin/python2.7 

import sys, os
import re
import numpy as np
import history_segments_module as histseg

def compare_histories_by_segment(segment_histories, history_ids, prevalence): 
	history_ids=open(history_ids, 'r')
	history_outputi={} # keys: a history ID, values: the column that the history will be output
	x=0
	for line in history_ids:
		history_id=line.strip()
		history_outputi[history_id]=x
		x+=1	
	nsimulations=x
	mysegstats=[]
	infile=open(segment_histories, 'r')
	for line in infile: 
		seg_hist=histseg.Segment_history(line)
		# values can be just the cnvals or just the prevals.  Distances between points in different histories will be calculated using the Euclidean distance between the given values. 
		fx=np.array(seg_hist.prevals) > prevalence
		values=np.array([seg_hist.cnvals, seg_hist.prevals])[:,fx]
		histid_crossindex={}  # keys: a history id, values= the column that this history is in 
		histid_i=0
		history_subset=[hist for hist, x in zip(seg_hist.histories, fx) if x]
		for id in history_subset : 
			if id not in histid_crossindex.keys(): 
				histid_crossindex[id] = histid_i
				histid_i += 1
		numhistories=len(histid_crossindex.keys())
		# For each event, calculate the minimum distance between that event and an event from each other history
		# distance array is n x m, n=number of events (points), m=number of histories
		distance_array=histseg.create_min_distance_array(values, history_subset, histid_crossindex)
		# history_similarities is m x m, but not symmetric
		history_similarities=histseg.create_similarity_matrix(distance_array, history_subset, histid_crossindex)
		history_scores=histseg.score_similarity_matrix(history_similarities)
		hist_diffs=np.empty(nsimulations)
		hist_diffs.fill(np.NAN)
		for id in histid_crossindex.keys(): 
			histid_i=histid_crossindex[id]
			outcol=history_outputi[id]
			hist_diffs[outcol]=history_scores[histid_i]
		mysegstats.append("%s\t%d\t%d\t%d\t%s\n" % (seg_hist.chr, seg_hist.start, seg_hist.end, np.sum(fx), "\t".join(map(str, hist_diffs))))
	return mysegstats


if __name__ == "__main__": 
	import argparse
	parser = argparse.ArgumentParser(description='Given a file of universal segments with histories (output from get_segment_values_dist.py), calculate the distance between histories.') 
	parser.add_argument('segment_histories', help='a file with a comma separated list of CN values, timings, prevalences, and history IDs') 
	parser.add_argument('history_ids', help='a list of the history_ids.  The order they are listed is the order that they will be columned in the output.')
	parser.add_argument('--prevalence', help='A cutoff value for the prevalence of an event.', type=float, default=0)
	args=parser.parse_args()
	history_ids=open(args.history_ids, 'r')
	history_outputi={} # keys: a history ID, values: the column that the history will be output
	x=0
	for line in history_ids:
		history_id=line.strip()
		history_outputi[history_id]=x
		x+=1	
	nsimulations=x
	differences = compare_histories_by_segment(args.segment_histories, args.history_ids, args.prevalence) 	
	for dat in differences: 
		sys.stdout.write(dat)
	
