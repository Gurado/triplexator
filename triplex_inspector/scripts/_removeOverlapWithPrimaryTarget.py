#!/usr/bin/python

import sys, os, glob
from optparse import OptionParser
import fileinput

global options
global args

def process():

	tts_ht = {}
	for line in fileinput.input([args[0]]):
		cols = line.strip().split("\t")
		if (len(cols)<4):
			continue
		try:
			tts_ht[cols[3]] = {'chrom':cols[0], 'start':int(cols[1]), 'end':int(cols[2])}
		except:
			if (option.verbose):
				print >> sys.stderr, 'skipping line TTS.file %s' % (line)
			
	removed = 0
	for line in fileinput.input([args[1]]):	
		if (line.startswith("#")):
			print line,
			continue
		cols = line.strip().split("\t")
		try:
			tts = cols[0]
			chrom = cols[3]
			start = int(cols[4])
			end = int(cols[5])
			# check if Sequence-ID of tts is known
			if (tts_ht.has_key(tts)):
				primary = tts_ht[tts]
				# skip lines where the target overlaps with primary one
				if (primary['chrom'] == chrom and (primary['start'] <= start and start <= primary['end']) or (primary['start'] <= end and end <= primary['end'])):
					removed += 1
					if (options.verbose):
						print >> sys.stderr, 'Remove off-target %s:%d-%d due to overlap with primary target region %s:%d-%d' % (chrom, start, end, primary['chrom'], primary['start'] , primary['end'] )

					continue
				else:
					print line,
			else:
				print >> sys.stderr, 'Sequence-ID not found in TTS-file %s' % (line)
		except:
			print >> sys.stderr, 'skipping line TPX.file %s' % (line)
		
	if (options.verbose):
		print >> sys.stderr, '%d lines removed' % (removed)
			
# manage option and arguments processing
def main():
	global options
	global args
	usage = '''usage: %prog [options] TTS.bed TPX.file

removes entries from the tpx file that overlap a primary target region
(provided in bed format). Writes to std-out.
	'''
	parser = OptionParser(usage)
	parser.add_option("-q", "--quiet", action="store_false", dest="verbose", default=True,
					help="don't print status messages to stdout")
	parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
					help="print status messages to stdout")
	
	(options, args) = parser.parse_args()
	if (len(args) != 2):
		parser.print_help()
		parser.error("incorrect number of arguments")
		
	if (options.verbose):
		print >> sys.stderr, "remove Overlap in %s (TPX file_ with %s (TTS file)" % (args[1], args[0])
		
	process()

	
if __name__ == "__main__":
	main()
