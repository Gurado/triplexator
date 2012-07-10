#!/bin/python

import fileinput, sys
from optparse import OptionParser

def print_entry(chrom, start, end, value):
	print "%s\t%d\t%d\t%.4f" % (chrom, start, end, value)
			
def process():
	chrom = ""
	step = 1
	start = 0
	counter = 0
	last = -1

	for line in sys.stdin.readlines():
		if (line.startswith("fixedStep")):
			if (chrom != ""):
				print_entry(chrom, start, counter, last)
				
			chrom = line.strip().split("chrom=")[1].split(" ")[0]
			step = int(line.strip().split("step=")[1].split(" ")[0])
			if (len(line.split("start="))>1):
				start = int(line.strip().split("start=")[1].split(" ")[0])
			else:
				start = 0
			counter = start
			last = -1
			
		else:
			v = float(line.strip())
			if (v != last):
				if (last != -1):
					print_entry(chrom, start, counter, last)
				start = counter
				last = v
			counter += step
	if (last != -1):
		print_entry(chrom, start, counter, last)

# manage option and arguments processing
def main():
	global options
	global args
	usage = '''usage: %prog [options] 
	
	converts a fixed step wig file from stdin to a bedgraph (stdout)
	'''
	parser = OptionParser(usage)
	parser.add_option("-q", "--quiet", action="store_false", dest="verbose", default=True,
		help="don't print status messages to stdout")
	parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
		help="print status messages to stdout")
	
	(options, args) = parser.parse_args()

	process()
	
if __name__ == "__main__":
	main()
