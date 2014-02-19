#!/inside/home/common/bin/python2.7
import os, sys
import event_cycles_module as histseg
import pickle



if __name__ == '__main__': 
	import argparse
	parser= argparse.ArgumentParser(description='Will unpickle a pickled events file.') 
	parser.add_argument('pickle', help='The pickled file to view')
	parser.add_argument('--historyid', type=int, help='The history you want to look at')
	args=parser.parse_args()
	events=pickle.load(open(args.pickle, 'rb'))
	if args.historyid: 
		hid=args.historyid
		sys.stderr.write("hid is %d\n" % hid)
		for event in events:
			if hid in event.histories: 
				sys.stdout.write(str(event))
	else: 
		for event in events:
			sys.stdout.write(str(event))
			sys.stdout.write(str(event.prevals) + "\n")
