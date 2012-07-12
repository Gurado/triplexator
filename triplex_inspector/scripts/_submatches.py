#!/usr/bin/python

import sys, os, traceback
from optparse import OptionParser
import fileinput
import math, re
from operator import itemgetter
from quicksect import IntervalNode

global options, args
global submatches, primaries, tfotargets, annotations, parameters, GGGtracksPattern

def find(start, end, tree):
    ''' Returns a list with the overlapping intervals '''
    out = []
    tree.intersect( start, end, lambda x: out.append(x) )
    return [ (x.start, x.end, x.linenum) for x in out ]

def processSubTarget(q, oId, intersect_tree, offtarget_data):
	# get statistics value for this submatch
	start = int(q[1])
	end = int(q[2])
	display_start = start
	display_end = end
							
	# in case the purine tract of the primary region is on the negative strand reverse the span
	if (primaries[oId]["strand"] == "-"):
		display_start = primaries[oId]["length"] - end
		display_end = primaries[oId]["length"] - start
	
	onLength = end-start
	onErrorRate = float(q[5])
	onError = int(round(onErrorRate*onLength))
	onGuanineRate = q[9].count('G')*1./onLength	
	onMaxMatchRun = 0
	onMatchRun = 0
	for i in range(len(q[9])):
		if (q[9][i] != 'G' and q[9][i] != 'A'): # error string
			if (onMatchRun > onMaxMatchRun):
				onMaxMatchRun = onMatchRun
			onMatchRun = 0
		else:
			onMatchRun += 1
	if (onMatchRun > onMaxMatchRun):
		onMaxMatchRun = onMatchRun

	# prepare annotation categories
	categories = [ {} for annoIndex in range(len(annotations)+1) ] # catalog for all annotation and irrespective of overlapping any annotation as well
	for annoIndex in range(len(annotations)+1): 
		categories[annoIndex]['shorterOffs'] = [ 0 for s in range(options.shorterUpTo ) ]
		categories[annoIndex]['erroneousOffs'] = [ 0 for s in range(options.errorsUpTo ) ]
		categories[annoIndex]['asGood'] = 0
    
    # get all the off-targets for this submatch
	overlap = find(start, end , intersect_tree)
	for offtarget in overlap:
		tfo_region = tuple(offtarget[0:2])
		# check minimum length (same for all off-targets)
		offLength = min(end, tfo_region[1])-max(start, tfo_region[0])
		if (offLength < parameters["minLength"]):
#					print >> sys.stderr, "%s %d-%d %d-%d < %d (minLength)" % (oId, start, end, offtarget[0], offtarget[1], offLength)
			continue # violates constraint - not an off-target

		# check individual off-targets
		for tfoErrors in offtarget_data[tfo_region].keys():
			
			# check guanine constraint, error rate and match-block
			g = 0.
			max_mb = 0
			mb = 0
			errors = 0
			error_plus = 0
			for i in range(max(start, tfo_region[0]), min(end, tfo_region[1])):
				if (tfoErrors[i] != ' '): # error string
					if (mb > max_mb):
						max_mb = mb
					mb = 0
					errors += 1
					if (tfoErrors[i] == 'd'  or tfoErrors[i] == 't' ):
						error_plus += 1
				else:
					mb += 1
					# no error, check for G in target
					if (q[9][i-start].upper()=="G"):
						g+=1.

			if (mb > max_mb):
				max_mb = mb

			if (errors*100. / offLength > parameters["errorRate"]):
	#					print >> sys.stderr, "%s %d-%d %d-%d < %.2f% (errorRate)" % (oId, start, end, offtarget[0], offtarget[1], errors*100. / offLength)
				continue # violates constraint - not an off-target
			
			offGuanineRate = g*100. / offLength 
			if (offGuanineRate < parameters["minGuanineRate"] or offGuanineRate > parameters["maxGuanineRate"]):
	#					print >> sys.stderr, "%s %d-%d %d-%d < %.2f (minGuanineRate)" % (oId, start, end, offtarget[0], offtarget[1], offGuanineRate)
				continue # violates constraint - not an off-target

			if (max_mb < parameters["matchBlock"]):
	#					print >> sys.stderr, "%s %d-%d %d-%d < %d (matchBlock)" % (oId, start, end, offtarget[0], offtarget[1], max_mb)
				continue # violates constraint - not an off-target
			
	    	# ok this is a valid off-target
	    	# catalog it!
	    						
			# check intersection with annotations
			for annoIndex in range(len(annotations)):
				if (offtarget_data[tfo_region][tfoErrors][1+annoIndex] > 0):
					if (error_plus == 0 and onLength-offLength == 0):
						categories[annoIndex]['asGood'] += offtarget_data[tfo_region][tfoErrors][1+annoIndex]
					elif (error_plus == 0 and onLength-offLength > 0 and onLength-offLength <= options.shorterUpTo):
						categories[annoIndex]['shorterOffs'][onLength-offLength-1] += offtarget_data[tfo_region][tfoErrors][1+annoIndex]
					elif (onLength-offLength == 0 and error_plus > 0 and error_plus <= options.errorsUpTo):
						categories[annoIndex]['erroneousOffs'][error_plus-1] += offtarget_data[tfo_region][tfoErrors][1+annoIndex]
	
			# and the default (irrespective of overlapping any annotations)
			if (error_plus == 0 and onLength-offLength == 0):
				categories[-1]['asGood'] += offtarget_data[tfo_region][tfoErrors][0]
			elif (error_plus == 0 and onLength-offLength > 0 and onLength-offLength <= options.shorterUpTo):
				categories[-1]['shorterOffs'][onLength-offLength-1] += offtarget_data[tfo_region][tfoErrors][0]
			elif (onLength-offLength == 0 and error_plus > 0 and error_plus <= options.errorsUpTo):
				categories[-1]['erroneousOffs'][error_plus-1] += offtarget_data[tfo_region][tfoErrors][0]

	# count number of GpGpG tracks
	GGGtracks = len(re.findall(GGGtracksPattern, q[9]));
	
	# output statistics
	output = ""
	for annoIndex in range(len(annotations)):				
		output += '\n\t["%s","%d-%d",%d,%.2f,%d,%d,%d,%.2f,"%s",%d,%s,%s],' % (oId, display_start, display_end, onLength, onGuanineRate, GGGtracks, onMaxMatchRun, onError, onErrorRate, annotations[annoIndex], categories[annoIndex]['asGood'], ','.join(["%s" % el for el in categories[annoIndex]['shorterOffs']]), ','.join(["%s" % el for el in categories[annoIndex]['erroneousOffs']]))
	output += '\n\t["%s","%d-%d",%d,%.2f,%d,%d,%d,%.2f,"%s",%d,%s,%s],' % (oId, display_start, display_end, onLength, onGuanineRate, GGGtracks,onMaxMatchRun, onError, onErrorRate, '-', categories[-1]['asGood'], ','.join(["%s" % el for el in categories[-1]['shorterOffs']]), ','.join(["%s" % el for el in categories[-1]['erroneousOffs']]))
	
	return output

def processTarget(oId, submatches, primaries, tfotargets, annotations, parameters):
	''' process one target region for all targets and their offtargets'''
	
	output = ""
	query  = submatches[oId] # list of submatches of the same id
	# quit if this primary target cannot be related to a primary target region
	if (not primaries.has_key(oId)):
		print >> sys.stderr, "[ERROR] submatch does not have a corresponding primary target : %s" % (oId)
		exit(1)
		
	# create an empty root
	intersect_tree = IntervalNode( -1, -1, -1 )
	# accumulate off-target data
	offtarget_data = []
	if (tfotargets.has_key(oId)):
		# put all off-targets into an interval tree:
		offtarget_data = tfotargets[oId]
	
	# build an interval tree from the rest of the data
	for tfo_region in offtarget_data.keys():
		start, end = tfo_region
		intersect_tree = intersect_tree.insert(start, end)
		
	# query all submatches
	output = ""
	for q in query:
		output += processSubTarget(q, oId, intersect_tree, offtarget_data)

	return output

def readData():
	global submatches, primaries, tfotargets, annotations, parameters, GGGtracksPattern
	
	GGGtracksPattern = re.compile('G{3}')
	
	# the first file provides information about the primary target sites
	primaries = {}
	for line in fileinput.input([args[0]]):
		cols = line.strip().split("\t")
		if (len(cols)<6):
			continue
		try:
			primaries[cols[3]] = {'chrom':cols[0], 'start':int(cols[1]), 'end':int(cols[2]), 'strand':cols[5], 'length':(int(cols[2])-int(cols[1])) }
		except:
			if (options.verbose):
				print >> sys.stderr, 'skipping line TTS.file %s' % (line)
	if (options.verbose):
		print >> sys.stderr, 'TTS.file read (%d entries)' % (len(primaries))

	# the second file provides all eligible submatches for a primary target site
	submatches = {}
	try:
		for line in fileinput.input([args[1]]):
			if (len(line.strip())==0 or line[0]=="#"):
				continue
			cols = line.rstrip('\n').split("\t")
			if (not submatches.has_key(cols[0])):
				submatches[cols[0]] = []
			submatches[cols[0]] += [cols]
		if (options.verbose):
			print >> sys.stderr, "TTS.submatches read (%d entries)"  % (len(submatches))
	except:
		print >> sys.stderr, "[ERROR] submatches file could not be read : %s" % (args[1])
		exit(1)
		
		
	# the third file contains the off-target informations for the primary site
	tfotargets = {}
	signature = [] # the signature provides information about the annotations used
	annotations = []
	try:
		for line in fileinput.input([args[2]]):
			# keep signature
			if (line.startswith("#")):
				signature = line.strip("#").strip().split("\t")
				continue
	
			elif (len(line.strip())==0):
				continue
				
			cols = line.strip().split("\t")
			tfoId = cols[0]
	
			tfoErrors = [" " for i in range(primaries[tfoId]['length'])]
			offLength = int(cols[5])-int(cols[4])
			
			# requires the error to be with respect to the TFO
			errors = cols[8]
			if (errors!=""):
				i = 0
				while (i<len(errors)):
					if (errors[i] in "dobt"):
						# get number
						j=i+1
						while(j<len(errors) and errors[j] in "0123456789"):
							j += 1 
						tfoErrors[int(cols[1])+int(errors[i+1:j])] = errors[i]  
						i=j
			tfoErrors = ''.join(tfoErrors)

			# add new TFO entry if required
			if (not tfotargets.has_key(tfoId)):
				tfotargets[tfoId] = {}
				
			# add new tfo_region if required
			tfo_region = tuple([int(cols[1]),int(cols[2])])
			if (not tfotargets[tfoId].has_key(tfo_region)):
				tfotargets[tfoId][tfo_region] = {}
				
			# add new error category (string) if required
			if (not tfotargets[tfoId][tfo_region].has_key(tfoErrors)):
				#                                         all entries & overlap with each annotation 
				tfotargets[tfoId][tfo_region][tfoErrors] = [0]+[0 for i in cols[12:]]

			# count this off-target				
			tfotargets[tfoId][tfo_region][tfoErrors][0] += 1
			
			# add annotation data
			for annoIndex in range(len(cols[12:])):
				if (cols[12+annoIndex]!='-'):
					tfotargets[tfoId][tfo_region][tfoErrors][1+annoIndex] += len(cols[12+annoIndex].split(","))
	
		if (options.verbose):
			print >> sys.stderr, "TRIPLEX.file read (%d entries)"  % (len(tfotargets))
	except:
		print >> sys.stderr, "[ERROR] off-target file could not be read : %s" % (args[2])
		traceback.print_exc(file=sys.stderr)
		exit(1)
	
	# get annotation naming
	annotations = signature[12:]
	
	# extract parameters from log file
	parameters = {}
	logfile = args[3]
	try:
		log = open(logfile).read()
		parameters["errorRate"] = int(log.split("- maximum error-rate :")[1].split("%")[0])
		parameters["minLength"] = int(log.split("- minimum length :")[1].split("nucleotides")[0])
		parameters["minGuanineRate"] = int(log.split("- minimum guanine content with respect to the target : ")[1].split("%")[0])
		parameters["maxGuanineRate"] = int(log.split("- maximum guanine content with respect to the target : ")[1].split("%")[0])
		parameters["matchBlock"] = int(log.split("- number of consecutive matches required in a feature :")[1].split("\n")[0])
		
	except:
		print >> sys.stderr, 'log file not found or in unexpected format: %s' % (logfile)
		sys.exit(1)
	
def process():
	''' processing entry point ''' 
	# get all data first
	readData();

	# get workload
	oIds = submatches.keys()
	oIds.sort()

	# prepare output variable
	output = ""

	for oId in oIds:
		output += processTarget(oId, submatches, primaries, tfotargets, annotations, parameters)
	
	# concatenate output
	output = '{ "aaData": [' + output[0:len(output)-1] + "\n],"
	# add columns
	output +='''
	"aoColumns": [
			{ "sTitle": "region Id", "sToolTip":"internal id used to identify regions of overlapping primary targets" },
			{ "sTitle": "target region", "sToolTip":"target-region spanned by the primary target" },
			{ "sTitle": "length", "sType": "numeric", "sToolTip":"number of nucleotides spanned by the primary target" },
			{ "sTitle": "G-ratio", "sType": "numeric", "sToolTip":"proportion of guanines in the purine tract of the target" }, 
			{ "sTitle": "GGG tracks", "sType": "numeric", "sToolTip":"number of GpGpG tracks in the target, which can be indicative for putative self-association" }, 
			{ "sTitle": "match run", "sType": "numeric", "sToolTip":"the longest continuous polypurine/polypurimidine run in the target" },
			{ "sTitle": "errors", "sType": "numeric", "sToolTip":"number of interruption in the target" },
			{ "sTitle": "error-rate", "sType": "numeric", "sToolTip":"number of interruption in the target normalized over length" },
			{ "sTitle": "annotation", "sToolTip":"type of annotation supplied by the user" },
			{ "sTitle": "as good", "sType": "numeric", "sToolTip":"number of off-targets that have the same affinity as the primary target" }'''
	for i in range(options.shorterUpTo):
		output += ',\n			{ "sTitle": "-%d nt(s)", "sType": "numeric", "sToolTip":"number of off-targets that are %d nucleotides shorter" }' % (i+1,i+1)
	for i in range(options.errorsUpTo):
		output += ',\n			{ "sTitle": "+%d error(s)", "sType": "numeric", "sToolTip":"number of off-targets that have %d more errors" }' % (i+1,i+1)
	# add table properties
	output += '''],
	"bPaginate": true,
	"bAutoWidth": false,
	"bInfo": true,
	"bLengthChange": true,
    "bFilter": true,
	"iDisplayLength": 10,
	"bJQueryUI": true
	'''
	output += '} \n'
	print output
		
			
# manage option and arguments processing
def main():
	global options
	global args
	usage = '''usage: %prog [options] TTS.bed submatches.TTS TRIPLEX.file TRIPLEX.log
	
Cross references all primary targets in a region with any off-target
subject to the constraints (also used in triplexator) and writes
a JSON object to standard out pointing out several statistics
of interest for designing molecules that target duplexes in a triplex manner.
	'''
	parser = OptionParser(usage)
	parser.add_option("-q", "--quiet", action="store_false", dest="verbose", default=True,
					help="don't print status messages to stdout")
	parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
					help="print status messages to stdout")
	parser.add_option("-S", "--shorter-up-to", type="int", dest="shorterUpTo", default=2, 
					help="how much shorter off-targets should appear in the table")
	parser.add_option("-E", "--errors-up-to", type="int", dest="errorsUpTo", default=1, 
					help="how much additional errors in off-targets should appear in the table")

	(options, args) = parser.parse_args()
	if (len(args) != 4):
		parser.print_help()
		parser.error("incorrect number of arguments")
		
	if (options.verbose):
		print >> sys.stderr, "proces primary targets from \n(TTS.bed) : %s\n(submatches.TTS) : %s\n(TRIPLEX.file) : %s\n(TRIPLEX.log) : %s" % (args[0],args[1],args[2],args[3])

	process()
	
if __name__ == "__main__":
	main()
	
