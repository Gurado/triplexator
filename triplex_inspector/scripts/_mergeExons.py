#!/usr/bin/python

import sys, os
from optparse import OptionParser
import fileinput

global options
global args

def process():
	
	last = []
	for line in sys.stdin.readlines():
			
		cols = line.strip().split("\t")
		# expect a 6 field bed
		if (len(cols)<6 and not cols[0].startswith("chr")):
			if (options.verbose):
				print >> sys.stderr, "[ERROR] %s" % (line)
			continue

		# convert genomic location into integer
		try:
			cols[1] = int(cols[1])
			cols[2] = int(cols[2])
		except:
			if (options.verbose):
				print >> sys.stderr, "[ERROR] %s" % (line)
			continue	
			
		# if this is the first entry ever
		if (len(last) == 0):
			last = cols
		# otherwise check id
		elif (cols[3] != last[3]):
			print "%s\t%d\t%d\t%s\t%s\t%s" % (last[0], last[1], last[2], last[3], last[4], cols[5])
			last = cols
		else:
			# merge entries if overlap
			if (last[5]==cols[5] and ((last[1] <= cols[1] and cols[1] <= last[2]) or (last[1] <= cols[2] and cols[2] <= last[2]))):
				last[1] = min(last[1], cols[1])
				last[2] = max(last[2], cols[2])
			# otherwise output last entry
			else:
				print "%s\t%d\t%d\t%s\t%s\t%s" % (last[0], last[1], last[2], last[3], last[4], cols[5])
				last = cols
	# don't forget the remaining entry
	if (len(last) >0 ):
		print "%s\t%d\t%d\t%s\t%s\t%s" % (last[0], last[1], last[2], last[3], last[4], cols[5])
			
		
# manage option and arguments processing
def main():
	global options
	global args
	usage = '''usage: %prog [options]  
	
Takes a 6-field bed file (sorted by col 4 = id) and merges entries with same name wrt strand
	'''
	
	parser = OptionParser(usage)
	parser.add_option("-q", "--quiet", action="store_false", dest="verbose", default=True,
					help="don't print status messages to stdout")
	parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
					help="print status messages to stdout")


	(options, args) = parser.parse_args()
	if (len(args) != 0):
		parser.print_help()
		parser.error("incorrect number of arguments")
		
	process()

	
if __name__ == "__main__":
	main()
	
