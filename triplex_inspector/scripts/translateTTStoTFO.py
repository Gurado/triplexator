#!/usr/bin/python

# takes all triplex result files in a folder and translate the TTS to TFOs of the pyrimidine motif

import sys, os, glob, fileinput
from optparse import OptionParser
from operator import itemgetter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

global options, args

# one of the three motifs is enough to define a counterpart, 
# however the GT motif is most resptrictive 
# since it is allowed to form a triplex in parallel as well as
# in antiparallel orientation
def translateTTS(tts):
	tfo =""
	for i in tts:
		if (i=='A' or i=='a' ):
			tfo += 'T'
		elif (i=='G' or i=='g'):
			tfo += 'C'
		else:
			tfo += 'N'
	return tfo

def process():

	records = {}
	if (options.fasta != ""):
		handle = open(options.fasta, "rU")
		records = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
		handle.close()	


	if (options.is_fasta):
		chrs = {}
		sid = ""
		tts = ""
		for line in sys.stdin.readlines():
			if (line[0]==">"):
				if (tts!=""):
					tfo = translateTTS(tts)
					print ">%s\n%s" % (sid, tfo)
					tts = ""
				sid = line.strip(">").strip()
				continue
			tts += line.strip()
		if (tts!=""):
			tfo = translateTTS(tts)
			print ">%s\n%s" % (sid, tfo)
			tts = ""
			if (options.verbose):
				print >> sys.stderr, "Entry %s translated" % (sid)	
				
	else:
		# read all bed entries and store with regard to chromosomes
		chrs = {}
		c = 0
		sid = ""
		tts = ""
		for line in sys.stdin.readlines():
			c += 1
			if (line[0]=="#"):
				continue
			
			# the TTS is encoded in the last column
			cols = line.strip().split("\t")
			strand = cols[5]
			sid = cols[3]
			if (len(records)==0):
				tts = cols[options.ttscol]
			else:
				try:
					tts = records[cols[0]][int(cols[1]):int(cols[2])]
					if (strand == "-"):
						tts = tts.reverse_complement()
				except:
					print >> sys.stderr, "[Error] could not find %s entry in the reference sequence file " % (cols[0])
				
			tfo = translateTTS(tts)
			print ">%s\n%s" % (sid, tfo)
			if (options.verbose):
				print >> sys.stderr, "Entry %s translated" % (sid)	
			
# manage option and arguments processing
def main():
	global options
	global args
	usage = '''usage: %prog [options]  

Translates a TTS result file (read from std-in)  
into a matching TFO file that can 
be used to screen a genome for putative off-targets.
	'''
	parser = OptionParser(usage)
	parser.add_option("-q", "--quiet", action="store_false", dest="verbose", default=True,
					help="don't print status messages to stdout")
	parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
					help="print status messages to stdout")
	parser.add_option("-c", "--column", type="int", dest="ttscol", default=-2, 
					help="column containing the TTS sequence")
	parser.add_option("-i", "--is-fasta", action="store_true", dest="is_fasta", default=False, 
					help="indicates that the input is a fasta rahter than a tts result file")
	parser.add_option("-f", "--fasta", type="string", dest="fasta", default="", 
					help="location with fasta sequences (overwrites any specified column)")
	
	(options, args) = parser.parse_args()
	if (len(args) != 0):
		parser.print_help()
		parser.error("incorrect number of arguments")

	process()
	
if __name__ == "__main__":
	main()
	
