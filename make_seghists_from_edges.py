#!/inside/home/common/bin/python2.7
import sys, os
import pickle
import history_segments_module as histseg


def create_seghists_from_edges(sortededges): 
	myseghists=[]
	working_seghists=[Seghist(sortededges[0])]
	for edge in sortededges[1:]:
		tmpseghists=[]
		for seghist in working_seghists: 
			overlapval=get_overlap_order(edge, seghist)
			if overlapval >0: # the edge overlaps the seghist or shares a breakend with it
				new_seghists=merge_in_edge(edge, seghist)
			 	tmpseghists+=new_seghists
			elif overlapval==-1: # the edge end comes before the seghist
			 	tmpseghists.append(seghist)
			elif overlapval==-2: # the edge start comes after the seghist
				myseghists.append(seghist)
		if overlapval==-1 or overlapval==2: # the edge extends or starts past the last segment	
			saved_edges.append(edge)
		working_seghists=tmpseghists
	
				
			

if __name__ == "__main__": 
	import argparse
	parser = argparse.ArgumentParser(description='Will take the edges from events and split the genome at breakpoints, with histories for each segment.')
	parser.add_argument('pegdes', help='a pickled files of event edges.')
	args=parser.parse_args()
	edges=pickle.load(open(args.pedges, 'rb'))
	sortededges=sort(edges, key=lambda x: (x.segstr, x.cnval))
	seghists=create_seghists_from_edges(sortededges)
	pickle.dump(seghists, stdout, pickle.HIGHEST_PROTOCOL)

	
