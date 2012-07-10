#!/usr/bin/python

import sys, os, glob
from optparse import OptionParser
import fileinput
import numpy as np

global options
global args

def process():
	
	c = 0
	ht_seq = {}
	ht_len = {}
	for line in fileinput.input([args[0]]):
		c += 1
		if (line[0]=="#"):
			continue
			
		cols = line.strip().split("\t")
		seq = cols[8]
		if (options.strict):
			seq = seq.replace("T","N")
			seq = seq.replace("t","N")
			seq = seq.replace("C","N")
			seq = seq.replace("c","N")
			
		if (not ht_seq.has_key(seq)):
			ht_seq[seq] = []
		ht_seq[seq] += [c]
		
		
		if (not ht_len.has_key(len(seq))):
			ht_len[len(seq)] = set()
		ht_len[len(seq)].add(seq)

	if (options.verbose):
		print >> sys.stderr, "read %d lines from file %s" % (c, args[0])


	minl = min(ht_len.keys())
	if (options.lbound > 0):
		minl = min(minl, options.lbound)
		
	if (options.verbose):
		print >> sys.stderr, "minimum length considered %d" % (minl)

	c=0
	copynumber = []
	tts_map = open(args[1],"w")	
	tts_fa = open(args[2],"w")
	for l in range(minl, max(ht_len.keys())+1):
		if (ht_len.has_key(l)):
			for seq in ht_len[l]:
				tts_fa.write(">tts_%d\n%s\n" % (c, seq))
				tts_map.write("tts_%d\t%s\t%d\t%s\n" % (c, seq, len(ht_seq[seq]), ",".join(["%s" % el for el in ht_seq[seq]])) )
				copynumber += [len(ht_seq[seq])]
				c+=1 
	tts_map.close()
	tts_fa.close()
	
	if (options.verbose):
		print >> sys.stderr, "mean copynumber %f" % (np.mean(copynumber))
		print >> sys.stderr, "median copynumber %f" % (np.median(copynumber))
		print >> sys.stderr, "stddev copynumber %f" % (np.std(copynumber))


# manage option and arguments processing
def main():
	global options
	global args
	usage = '''usage: %prog [options] TTS TTS.map TTS.fasta

takes the TTS file and collapses all duplicates writing a map (TTS_map) and 
a fasta file
	'''
	parser = OptionParser(usage)
	parser.add_option("-q", "--quiet", action="store_false", dest="verbose", default=True,
					help="don't print status messages to stdout")
	parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
					help="print status messages to stdout")
	parser.add_option("-l", "--lower-bound", type="int", dest="lbound", default=0, 
					help="minimum length of tts")
	parser.add_option("-s", "--strict", action="store_true", dest="strict", default=False, 
					help="convert T and C to N")
	
	(options, args) = parser.parse_args()
	if (len(args) != 3):
		parser.print_help()
		parser.error("incorrect number of arguments")
		
	datafolder = os.path.abspath(args[0])

	if (options.verbose):
		print >> sys.stderr, "TTS file  : %s" % (args[0])
		print >> sys.stderr, "TTS.map file  : %s" % (args[1])
		print >> sys.stderr, "TTS.fasta file  : %s" % (args[2])
		
	process()

	
if __name__ == "__main__":
	main()
