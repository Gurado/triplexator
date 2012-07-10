#!/usr/bin/python

import sys, os
from optparse import OptionParser
import fileinput

global options
global args

def process():
	'''create lookup table for annotation'''
	
	lookup = {}
	for line in fileinput.input([args[1]]):
		cols = line.strip().split("\t")
		if (len(cols) < 12):
			if (options.verbose):
				print >> sys.stderr, 'skipping %s' % (line)
			continue
		try:
			line = int(cols[3])
			anno = cols[9]
		except:
			if (options.verbose):
				print >> sys.stderr, 'skipping %s' % (line)
			continue
			
		if (not lookup.has_key(line)):
			lookup[line] = []
		lookup[line] += [anno]
	
	# process tpx file
	linecount = 0
	for line in fileinput.input([args[0]]):
		linecount+=1
		line = line.strip()
		if (line.startswith("#")):
			print "%s\t%s" % (line, args[2])
			
		elif (lookup.has_key(linecount)):
			print "%s\t%s" % (line, ",".join(lookup[linecount]))
		else:
			print "%s\t-" % (line)
		
def main():
	''' manage option and arguments processing'''
	
	global options
	global args
	usage = '''usage: %prog [options] tpx.bed annotation.intersect column_name
	
Adds a column containing annotation data to the tpx file 
	'''
	
	parser = OptionParser(usage)
	parser.add_option("-q", "--quiet", action="store_false", dest="verbose", default=True,
					help="don't print status messages to stdout")
	parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
					help="print status messages to stdout")


	(options, args) = parser.parse_args()
	if (len(args) != 3):
		parser.print_help()
		parser.error("incorrect number of arguments")
		
	process()

	
if __name__ == "__main__":
	main()
	
