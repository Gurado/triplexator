#!/usr/bin/python

import sys, traceback, os, glob
from optparse import OptionParser
import fileinput
from Bio import SeqIO

global options
global args

def process():

	# read reference bed file
	data = {}
	for line in fileinput.input([args[0]]):
		cols = line.strip().split("\t")
		if (len(cols)<6):
			continue
		try:
			data[cols[3]] = (cols[0],int(cols[1]),int(cols[2]), cols[5])
		except:
			print >> sys.stderr, "error in line:	%s" % (line),
		
	if (options.verbose):
		print >> sys.stderr, "reference file contains %d entries " % (len(data))
		
	# process TTSs from stdin
	for line in sys.stdin.readlines():
		cols = line.rstrip('\n').split("\t")
		if (len(cols)<3):
			continue
			
		key = cols[0]
		if (data.has_key(key)):
			try:
				if (options.bed==3 and data[key][3]=="+"):
					print "%s\t%d\t%d\t%s\t0\t+\t%s\n" % (data[key][0], data[key][1]+int(cols[1]), data[key][1]+int(cols[2]), key, cols[0])  ,
				elif (options.bed==3 and data[key][3]=="-"):
					print "%s\t%d\t%d\t%s\t0\t+\t%s\n" % (data[key][0], data[key][2]-int(cols[2]), data[key][2]-int(cols[1]), key, cols[0])  ,
				elif (options.bed==6 and data[key][3]=="+"):
					print "%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n" % (data[key][0], data[key][1]+int(cols[1]), data[key][1]+int(cols[2]), cols[3], cols[4], cols[5], cols[0], cols[6], cols[7])  ,
				elif (options.bed==6 and data[key][3]=="-" and cols[5]=="+"):
					print "%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n" % (data[key][0], data[key][2]-int(cols[2]), data[key][2]-int(cols[1]), cols[3], cols[4], "-", cols[0], cols[6], cols[7]) ,
				elif (options.bed==6 and data[key][3]=="-" and cols[5]=="-"):
					print "%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n" % (data[key][0], data[key][2]-int(cols[2]), data[key][2]-int(cols[1]), cols[3], cols[4], "+", cols[0], cols[6], cols[7]) ,
			except:
				print  >> sys.stderr, "something went wrong\n%s\n%s\n%s" % (key, cols, data[key])
				sys.exit(1)
				
# manage option and arguments processing
def main():
	global options
	global args
	usage = '''usage: %prog [options] LOI.bed
	
Maps a TTS result file to the genomic coordinates.
Reads std-in, writes std-out.
'''
	parser = OptionParser(usage)
	parser.add_option("-q", "--quiet", action="store_false", dest="verbose", default=True,
					help="don't print status messages to stdout")
	parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
					help="print status messages to stdout")
	parser.add_option("-b", "--bed", type="int", dest="bed", default=3, 
					help="bed-format used")
	
	(options, args) = parser.parse_args()
	if (len(args) != 1):
		parser.print_help()
		parser.error("incorrect number of arguments")
		
	datafolder = os.path.abspath(args[0])

	if (options.verbose):
		print >> sys.stderr, "convert source (LOI.bed) : %s" % (args[0])

	process()
	
if __name__ == "__main__":
	main()
