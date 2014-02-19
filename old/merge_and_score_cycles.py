#!/inside/home/common/bin/python2.7
import sys, os
import event_cycles_module as eventmod 

if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser(description='Given a *.evnt file (contains one event per line), which is sorted by the first location within the event, this will merge identical events and assign them a probability score.')
	parser.add_argument('evnt', help='the *.evnt file, the output from braney_to_event_per_line.py', type=argparse.FileType('r'), default='-')
	parser.add_argument('--totp', help='The total probability, default is 1.', default=1)
	args=parser.parse_args()
	evntsfile=args.evnt
	unique_events=eventmod.unique_c_events_sorted_file(evntsfile)
	for event in unique_events: 
		event.get_hranges()
		event.compute_likelihood(args.totp)
		sys.stdout.write(str(event))

