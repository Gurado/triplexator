#!/usr/bin/python

import sys, os
from optparse import OptionParser
import fileinput

global options
global args

def process():
	
	for line in sys.stdin.readlines():
			
		cols = line.strip().split("\t")
		
		if (cols[2]!="exon"):
			continue
		
		exonnr = "."
		try:	
			geneid = cols[8].split('"')[1]
#			exonnr = cols[8].split('"')[5]
		except:
			print >> sys.stderr, cols[8]
			
		print "%s\t%d\t%d\t%s\t%s\t%s" % (cols[0], int(cols[3]), int(cols[4]), geneid, exonnr, cols[6])
		
			
# manage option and arguments processing
def main():
	global options
	global args
	usage = '''usage: %prog [options]  
	
Takes a gtf file and converts it to
a bed format keeping only exon entries
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
	
