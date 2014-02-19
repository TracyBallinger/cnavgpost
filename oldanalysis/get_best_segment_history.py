#!/inside/home/common/bin/python2.7 

import sys, os
import re
import argparse
import numpy as np
import history_segments_module as histseg

parser = argparse.ArgumentParser(description='Given a file of universal segments with histories (output from get_segment_values_dist.py), this will add an 8th column for the best history per segment') 
parser.add_argument('segment_histories', help='a file with a comma separated list of CN values, timings, prevalences, and history IDs') 

args=parser.parse_args()
infile=open(args.segment_histories, 'r')

def get_best_history(seg_hist):
	# values can be just the cnvals or just the prevals.  Distances between points in different histories will be calculated using the Euclidean distance between the given values. 
    values=np.array([seg_hist.cnvals, seg_hist.prevals])
    histid_crossindex={}  # keys: a history id, values= the column that this history is in 
    histid_i=0
    for id in seg_hist.histories :
        if id not in histid_crossindex.keys():
            histid_crossindex[id] = histid_i
            histid_i += 1
    numhistories=len(histid_crossindex.keys())
    # For each event, calculate the minimum distance between that event and an event from each other history
    # distance array is n x m, n=number of events (points), m=number of histories
    distance_array=histseg.create_min_distance_array(values, seg_hist.histories, histid_crossindex)
	# history_similarities is m x m
    history_similarities=histseg.create_similarity_matrix(distance_array, seg_hist.histories, histid_crossindex)
    history_scores=histseg.score_similarity_matrix(history_similarities)
    best_history=np.argmin(history_scores)
    best_historyid=[id for id, column in histid_crossindex.items() if column==best_history][0]
    return(best_historyid)

for line in infile: 
	myseghist=histseg.Segment_history(line)
	best_history=get_best_history(myseghist)
	sys.stdout.write("%s\t%s\n" % (str(myseghist), str(best_history)))		


