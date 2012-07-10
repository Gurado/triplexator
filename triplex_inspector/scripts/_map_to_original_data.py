#!/usr/bin/python

import sys, os
from optparse import OptionParser
import fileinput

global options
global args

def process():

	bed_ht = {}
	
	if (options.reference != ""):
		for line in fileinput.input([options.reference]):
			cols = line.strip().split("\t")
			bed_ht[cols[3]] = [cols[0],int(cols[1]), int(cols[2]), cols[3], int(cols[4]), cols[5]]
		if (options.verbose):
			print >> sys.stderr, "read %d entries from reference %s" % (len(bed_ht), options.reference)
			
	counter = 0
	for line in fileinput.input([args[0]]):
		cols = line.strip().split("\t")
		if (len(cols)<4 or line.startswith("#")):
			continue
		
		if (len(bed_ht) > 0):
			if (not bed_ht.has_key(cols[0])):
				print >> sys.stderr, "sequence not known"
			entry = bed_ht[cols[0]]
		else:
			entry = [cols[0], 0, 0, counter, 1 ,"+"]
			counter += 1 
			
		# discard too long features
		if (options.max_length>0 and int(cols[2])-int(cols[1]) > options.max_length):
			continue
			
		print "%s\t%d\t%d\t%s\t%d\t%s" % (entry[0], entry[1]+int(cols[1]), entry[1]+int(cols[2]), entry[3], entry[4], entry[5])
		



# manage option and arguments processing
def main():
	global options
	global args
	usage = '''usage: %prog [options] original.bed triplexator.result
	
Takes an result file from triplexator and the bed file containing the original coordinates 
and convert the coordinates back to the original sequences

prints stdandard out
	'''
	
	parser = OptionParser(usage)
	parser.add_option("-q", "--quiet", action="store_false", dest="verbose", default=True,
					help="don't print status messages to stdout")
	parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
					help="print status messages to stdout")
	parser.add_option("-m", "--max-length-filter", type="int", dest="max_length", default=0, 
					help="maximal length a triplex feature is allowed to be (long features likely contain repeats or secondary sites)")
	parser.add_option("-r", "--reference-bed", type="string", dest="reference", default="", 
					help="bed file that contains the original sequence regions")

	(options, args) = parser.parse_args()
	if (len(args) != 1):
		parser.print_help()
		parser.error("incorrect number of arguments")
		
	process()

	
if __name__ == "__main__":
	main()
	
