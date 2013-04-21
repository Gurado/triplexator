#!/usr/bin/python

import sys, os, traceback
from optparse import OptionParser
import fileinput
from operator import itemgetter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import json
import math
from quicksect import IntervalNode
import numpy as np

global options
global args
global annotation_palettes

#---------------
def find(start, end, tree):
    "Returns a list with the overlapping intervals"
    out = []
    tree.intersect( start, end, lambda x: out.append(x) )
    return [ (x.start, x.end, x.linenum) for x in out ]


#---------------
def readdata():
	""" read all supplied data files """
	# read fasta if supplied
	records = {}
	if (options.fastaFile != ""):
		handle = open(options.fastaFile, "rU")
		records = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
		handle.close()	
		if (options.verbose):
			print >> sys.stderr, "FASTA: read %d entries"  % (len(records))
			
	# read region of interest BED file
	loiBed = {}
	try:
		for line in fileinput.input([args[0]]):
			if (len(line.strip())==0 or line[0]=="#"):
				continue
			cols = line.strip().split("\t")
			loiBed[cols[3]] = cols
		if (options.verbose):
			print >> sys.stderr, "LOI BED: read %d entries"  % (len(loiBed))
	except:
		print >> sys.stderr, "[ERROR] LOI BED could not be read : %s" % (args[0])
		exit(1)
		
	# read regions to consider for off-targets (BED file)
	offBed = {}
	try:
		for line in fileinput.input([args[1]]):
			if (len(line.strip())==0 or line[0]=="#"):
				continue
			cols = line.strip().split("\t")
			offBed[cols[3]] = [cols[0], int(cols[1]), int(cols[2])]+cols[3:]
		if (options.verbose):
			print >> sys.stderr, "Off BED: read %d entries"  % (len(offBed))
	except:
		print >> sys.stderr, "[ERROR] Off BED could not be read : %s" % (args[1])
		exit(1)

	# read putative on-target region file (triplexator tts file)
	targets = {}
	try:
		for line in fileinput.input([args[2]]):
			if (len(line.strip())==0 or line[0]=="#"):
				continue
			cols = line.rstrip('\n').split("\t")
			if (not targets.has_key(cols[6])):
				targets[cols[6]] = []
			targets[cols[6]] += [cols]
		if (options.verbose):
			print >> sys.stderr, "TTS: read %d loci of interest"  % (len(targets))
	except:
		print >> sys.stderr, "[ERROR] primary target file could not be read : %s" % (args[2])
		exit(1)
		
	# read off-targets file (triplexator tpx file augmented with annotation data)
	tfotargets = {}
	tfoerrors = {} # position with errors
	signature = [] # signature provides information about the annotation used
	annotations = []
	linecount = 0
	try:
		for line in fileinput.input([args[3]]):
			cols = line.strip().split("\t")
			
			# keep signature
			if (line.startswith("#")):
				signature = line.strip("#").strip().split("\t")
				continue
		
			# skip non-qualifying lines
			elif (len(cols) <= 8):
				continue

			linecount += 1
	
			# convert line columns into proper datatypes
			if (not tfotargets.has_key(cols[0])):
				tfotargets[cols[0]] = []
				
			tfotargets[cols[0]] += [tuple([cols[0]]+[int(cols[1])]+[int(cols[2])]+[cols[3]]+[int(cols[4])]+[int(cols[5])]+cols[6:])]
			
			if (options.verbose and linecount % 100000 == 0):
				print >> sys.stderr, "... reading (%d entries done)"  % (linecount)

			
		if (options.verbose):
			print >> sys.stderr, "TRIPLEX: read %d offtarget entries for %d primary targets"  % (linecount, len(tfotargets))
			
		# sort off-targets according to location
		for target in tfotargets.keys():
			tfotargets[target] = sorted(tfotargets[target], key=itemgetter(3,4))
		
		if (options.verbose):
			print >> sys.stderr, "TRIPLEX: offtargets sorted"
			
	except:
		print >> sys.stderr, "[ERROR] off-target file could not be read : %s" % (args[3])
		exit(1)

	# annotations are located after position 12
	annotations = signature[12:]

	# read primary target qualifying subregions and its number of equivalent off-targets
	submatches = {}
	try:
		jd = json.loads(open(args[4],'r').read())["aaData"]
		for entry in jd:
			if (not submatches.has_key(entry[0])):
				submatches[entry[0]] = []
			span = entry[1].split("-")
			# get only entries where there is no annotation limit
			if (entry[8]=="-"):
				submatches[entry[0]] += [[int(span[0]),int(span[1]), int(entry[9])]] # start, end, equivalent off-targets
			
		if (options.verbose):
			print >> sys.stderr, "Submatches: read %d submatches"  % (len(jd))
	except:
		print >> sys.stderr, "[ERROR] submatches file could not be read : %s" % (args[4])
		exit(1)
		
	# read chromosome sizes
	chromsizes = []
	try:
		for line in fileinput.input([args[5]]):
			cols = line.strip().split()
			if (len(cols)>=2):
				chromsizes += [(cols[0], int(cols[1]))]
	except:
		print >> sys.stderr, "[ERROR] could not read chromosome sizes : %s" % (args[5])
		exit(1)
		
	return (records, loiBed, offBed, targets, tfotargets, annotations, submatches, chromsizes);
	
#---------------
def _make_position_array(features, oId, tlength, tstrand, annotation):
	"""
	translates the error information that are with respect to the off-target sites
	to the primary site
	"""
	
	positions = {}
	for t in "-dobt":
		positions[t] = [0 for r in range(tlength)] 

	for feature in features:
		errors = feature['Errors']
		ostrand = feature['Strand']
		# skip entry if it is not linked to requested annotation
		if (annotation != "-" and feature[annotation]=="-"):
			continue
		tmp_positions = [ "." for r in range(tlength)] # place holder 
		for i in range(feature['TFO start'], feature['TFO end']):
			tmp_positions[i] = "-" # indicates a match

		if (errors!=""):
			i = 0
			while (i<len(errors)):
				if (errors[i] in "dobt"):
					# get position (number)
					j=i+1
					while(j<len(errors) and errors[j] in "0123456789"):
						j += 1 
					tmp_positions[feature['TFO start']+int(errors[i+1:j])] = errors[i] # requires the error to be with respect to the TFO
					i=j
		
		# reverse list if primary target is located on the negative strand
		if (tstrand=="-"):
			tmp_positions = tmp_positions[::-1]
		# add the hits to the counter 
		for pos in range(tlength):
			if (tmp_positions[pos] != "."):
				positions[tmp_positions[pos]][pos] += 1	
	return positions
	
#---------------
def processChromatin(offBed, f_configuration):
	""" process chromatin data into blocks of specific resolution """

	if (options.verbose):
		print "Start chromatin processing"
	
	if (not f_configuration["conf"].has_key("genomesize")):
		calculateChromosomeSpecifics(chromsizes, f_configuration)
	step_size = int(f_configuration["conf"]["genomesize"] * options.chromatinResolution / (360. - f_configuration["conf"]["target_radius"] )) 

	if (options.verbose):
		print "chromatin stepsize: %d nts" % (step_size)
			
	samfile = ""
	bigwig = ""
	
	# check if chromatin file has been supplied
	if (options.chromatinFormat == "bam"):
		import pysam # lazy import of required module
		# try open the bam file 
		try:
			samfile = pysam.Samfile(options.chromatin, "rb" )
		except:
			print >> sys.stderr, "[ERROR] chromatin file (bam) could not be read : %s" % (options.chromatin)
			traceback.print_exc()
			exit(1)
	elif (options.chromatinFormat == "bigwig"):
		from bx.intervals.io import GenomicIntervalReader # lazy import of required module
		from bx.bbi.bigwig_file import BigWigFile  # lazy import of required module
		# try open the bigwig file
		try:
			bigwig = BigWigFile(open(options.chromatin))
		except:
			print >> sys.stderr, "[ERROR] chromatin file (bigwig) could not be read : %s" % (options.chromatin)	
    			traceback.print_exc()
			exit(1)
	
	if (options.verbose):
		print "... chromatin file opened"
	
	max_value = 0.
	
	# open output file	
	if (os.path.exists(options.output_dir+"chromatin.txt")):
		if (options.verbose):
			print "using existing chromatin file"
		for line in fileinput.input([options.output_dir+"chromatin.txt"]):
			c_value = float(line.strip().split("\t")[3])
			
			if (c_value > max_value):
				max_value = c_value
				
	else:
		f_chromatin = open(options.output_dir+"chromatin.txt","w")
		# process calculate mean chromatin density value for region
		for region in offBed.values():
			chrom = region[0]
			if (options.verbose):
				print "Collecting chromatin data for chromosom %s" % (chrom)
			
			
			for span_start in range(region[1], region[2], step_size):
				span_end = min(span_start + step_size, region[2])
				c_value = 0.
				
				if (span_start == span_end):
					continue
				
				if (options.chromatinFormat=="bam"):
					try:
						on_chromatin_features = 0.
						for pileupcolumn in samfile.pileup(chrom,span_start,span_end):
							on_chromatin_features += pileupcolumn.n
						c_value = on_chromatin_features/(span_end-span_start)
					except:
						if (options.verbose):
							print >> sys.stderr, "[WARN] bam file exception : %s:%d-%d" % (chrom,span_start,span_end)
				elif (options.chromatinFormat=="bigwig"):
					try:
						bwsummary = bigwig.query(chrom,span_start,span_end, 1 ) 
						c_value = bwsummary[0]["mean"]
					except:
						if (options.verbose):
							print >> sys.stderr, "[WARN] bigwig file exception : %s:%d-%d" % (chrom,span_start,span_end)

				if (math.isnan(c_value)):
					c_value = 0.

				# add pseudocount to circumvent log issues
				c_value += 0.001
				
				if (c_value > max_value):
					max_value = c_value

				f_chromatin.write("%s\t%d\t%d\t%.3f\n" % (chrom, span_start,span_end, c_value))
		
		# close filehandle
		f_chromatin.close()
	
	if (options.verbose):
		print "Maximal chromatin value: %.3f" % (max_value)
	
	# add track to configuration file
	if (not f_configuration.has_key("plots")):
		f_configuration["plots"] = {}

	f_configuration["plots"]["type"] = "heatmap"
	f_configuration["plots"]["color"] = "spectral-9-div"
	f_configuration["plots"]["stroke_thickness"] = "1"
	f_configuration["plots"]["stroke_color"] = "black"
		
	if (not f_configuration["plots"].has_key("plot")):
		f_configuration["plots"]["plot"] = []
	
	stroke_thickness = 0
	if (options.chromatinResolution >= 0.5):
		stroke_thickness = 1


	f_configuration["plots"]["plot"] += ['''
<plot>
	show             = conf(show_histogram)
	type             = histogram
	fill_under       = yes
	fill_color       = lgrey
	file             = %s
	r0               = 0.95r
	r1               = 0.997r
	scale_log_base   = 5
	stroke_thickness = %d
	color			 = black
	min              = 0.0
	max              = %.3f
	
	<axes>
		show = data
		thickness = 1
		color = lgrey
		<axis>
			spacing     = 0.1r
		</axis>
	</axes>
</plot>
''' % (options.output_dir+"chromatin.txt", stroke_thickness, max_value)]
	
	if (options.verbose):
		print "Finished chromatin processing"

#---------------
def outputErrors(loi, targetid, chrom, tstart, tstop, tsize, cstrand, errorstring, f_configuration):
	""" output hightlights where errors occur in the primary target """
	
	# mark errors in primary target
	f_ontarget_highlight = open(options.output_dir+targetid+os.sep+"primarytarget_highlight.txt","w")
	if (errorstring!=""): # errors in the primary target
		for s in errorstring.strip('d').split('d'):
			if (cstrand == "+"):
				f_ontarget_highlight.write("target\t%d\t%d\n" % (int(s), int(s)+1))
			else:
				f_ontarget_highlight.write("target\t%d\t%d\n" % (tsize-int(s)-1, tsize-int(s)))

	f_ontarget_highlight.close()
	
	f_configuration[targetid]["highlight"] += ['''
<highlight>
show       = conf(show_highlights)
file       = %s
ideogram   = yes
fill_color = red_nontriad
</highlight>

<highlight>
show       = conf(show_highlights)
file       = %s
r0         = 0.999r
r1         = 0.21r
fill_color = red_nontriad_faint
ideogram   = no
</highlight>
''' % (options.output_dir+targetid+os.sep+"primarytarget_highlight.txt", options.output_dir+targetid+os.sep+"primarytarget_highlight.txt")]


#---------------
def writeTileData(submatches, targetid, trackfile, f_configuration):
	""" adds the tile plots to the configuration file  """
	
	if (not f_configuration.has_key("tiles")):
		f_configuration["tiles"] = {}
		f_configuration["tiles"]["total_space"] 	=  0.7
		f_configuration["tiles"]["start_radius"]	=  0.91
		f_configuration["tiles"]["lane_padding"]	=  0.006
		f_configuration["tiles"]["tile_padding"]	=  0.001
	
	if (not f_configuration["tiles"].has_key(targetid)):
		f_configuration["tiles"][targetid] = {}
		f_configuration["tiles"][targetid]["region"] = {}
		

	# sample subtarget range
	spans = {} 
	lanes = {}
	span_lookup = {} # lookup table for quick access of submatches of a given length
	for i in range(len(submatches)):
		submatch = submatches[i]
		span = submatch[1]-submatch[0]
		if (not spans.has_key(span)):
			spans[span] = 0
			lanes[span] = 0
			span_lookup[span] = []
			
		lanes[span] += 1
		spans[span] += submatch[2]
		span_lookup[span] += [i]

	# remove all spans without offtargets execpt the shortest one
	spankeys = sorted(list(spans.keys()))
	try:
		while (spans[spankeys[-1]]==0 and spans[spankeys[-2]]==0):
			spankeys.pop()
	except:
		pass

	# pick at most 10 lanes 
	spankeys = spankeys[-10::]
	
	if (options.verbose):
		print "tile tracks for length: %s" % (" ".join(["%s" % el for el in spankeys]))

	# number of total tile lanes
	total_lanes = 0.
	for span in spankeys:
		total_lanes += min(span, lanes[span])
	
	total_space 	= f_configuration["tiles"]["total_space"]
	current_radius 	= f_configuration["tiles"]["start_radius"]
	lane_padding	= f_configuration["tiles"]["lane_padding"]
	tile_padding    = f_configuration["tiles"]["tile_padding"]
	total_radius    = f_configuration["image"]["radius"]
	tilespace       = total_radius*total_space
	
	#           total-space - (number-of-tracks * space-between-tracks) / total-lanes
	space_per_lane = (total_space - (len(spankeys)-1)*lane_padding - tile_padding*(total_lanes-len(spankeys))) / total_lanes
	

	f_configuration["tiles"][targetid]["thickness"] = math.floor((total_radius*0.94-30)*space_per_lane ) # 0.94 & -30 correspond to ideogram radius and height
	
	if (options.verbose):
		print "radius:%d tile-space: %.0f lanes: %d lane-padding: %.1f tile-padding:%.1f total-lanes:%d tile-thickness:%.1f space-per-lane:%.4f" % (total_radius, tilespace, len(spankeys), ((len(spankeys)-1)*lane_padding)*total_radius, tile_padding*total_lanes*total_radius, total_lanes, f_configuration["tiles"][targetid]["thickness"], space_per_lane)
	
	for span in spankeys:
		layers = min(span, lanes[span])
		r1 = current_radius
		r0 = current_radius-(layers*space_per_lane)-tile_padding*layers

		if (options.verbose):
			print "tile %d layers: %d" % (span, layers)

		f_configuration[targetid]["plot"] += ['''# tile plot
<plot>
show   = conf(show_tiles)
type   = tile
file   = %s

layers = %d
layers_overflow = collapse
layers_overflow_color = blue

thickness = %.fp
padding = %dp
orientation = in
stroke_thickness = 1
stroke_color     = black
r0    = %.3fr
r1    = %.3fr
color = rdylgn-11-div-1

<rules>
<rule>
importance = 100
condition = _SIZE_ != %d
show = no
</rule>

<rule>
importance = 99
condition = _VALUE_ < 1
color = rdylgn-11-div-11
</rule>

<rule>
importance = 98
condition = _VALUE_ <= 5
color = rdylgn-11-div-10
</rule>

<rule>
importance = 97
condition = _VALUE_ <= 10
color = rdylgn-11-div-9
</rule>

<rule>
importance = 96
condition = _VALUE_ <= 50
color = rdylgn-11-div-8
</rule>

<rule>
importance = 95
condition = _VALUE_ <= 100
color = rdylgn-11-div-7
</rule>

<rule>
importance = 94
condition = _VALUE_ <= 500
color = rdylgn-11-div-6
</rule>

<rule>
importance = 93
condition = _VALUE_ <= 1000
color = rdylgn-11-div-5
</rule>

<rule>
importance = 92
condition = _VALUE_ <= 5000
color = rdylgn-11-div-4
</rule>

<rule>
importance = 91
condition = _VALUE_ <= 10000
color = rdylgn-11-div-3
</rule>

<rule>
importance = 90
condition = _VALUE_ <= 50000
color = rdylgn-11-div-2
</rule>

</rules>

<<include %sbackground.conf>>
</plot>
''' % (trackfile, layers, f_configuration["tiles"][targetid]["thickness"] , tile_padding*total_radius, r0, r1, span+1, options.output_dir)]

		# save specifics
		f_configuration["tiles"][targetid]["region"][span] = {"radius":(r0, r1), "layers":layers}
		
		# update radius
		current_radius = r0-lane_padding
		
	# add ticks to plot
	writeTicks(targetid, f_configuration)


#---------------
def writeNuclotidePositionData(loi, targetid, chrom, tstart, tstop, tsize, cstrand, positions, f_configuration):
	""" output nucleotide position absolute & relative """
	
	f_ontarget_positions_abs = open(options.output_dir+targetid+os.sep+"nucleotide_abs.txt","w")
	f_ontarget_positions_rel = open(options.output_dir+targetid+os.sep+"nucleotide_rel.txt","w")
	max_v = 0.
	for pos in range(tsize):
		f_ontarget_positions_rel.write("target\t%d\t%d\t" % (pos, pos+1))
		dist = []
		for t in "dobt-":
			dist += [positions[pos][t]]
			
		sumdist = float(sum(dist))
		
		if (sumdist > max_v):
			max_v = sumdist
		
		for dpos in range(len(dist)):
			dist[dpos] = math.ceil(100.*dist[dpos]/sumdist)/100.
		
		f_ontarget_positions_rel.write("%s\n" % (",".join(["%.2f" % el for el in dist])))
	
	# normalize absolute counts with max_v (to %)
	for pos in range(tsize):
		f_ontarget_positions_abs.write("target\t%d\t%d\t" % (pos, pos+1))
		dist = []
		for t in "dobt-":
			dist += [positions[pos][t]]
			
		for dpos in range(len(dist)):
			dist[dpos] = math.ceil(100.*dist[dpos] /max_v )/100.
		
		f_ontarget_positions_abs.write("%s\n" % (",".join(["%.2f" % el for el in dist])))

	# workaround to avoid rotating plot to start 		
	f_ontarget_positions_abs.write("target\t%d\t%d\t" % (tsize, tsize))
	f_ontarget_positions_rel.write("target\t%d\t%d\t" % (tsize, tsize))
		
	f_ontarget_positions_abs.write("%s\n" % (",".join(["0.00" for t in "dobt-"])))
	f_ontarget_positions_rel.write("%s\n" % (",".join(["0.00" for t in "dobt-"])))
	
	f_ontarget_positions_rel.close()
	f_ontarget_positions_abs.close()

	if (not f_configuration[targetid].has_key("plot")):
		f_configuration[targetid]["plot"] = []

	f_configuration[targetid]["plot"] += ['''
<plot>
file = %s
show = conf(show_histogram)
type = histogram
r0 = 0.95r
r1 = 0.997r
color      = black
fill_color = red,orange,yellow,green,blue
fill_under = yes
thickness  = 1
sort_bin_values = no
extend_bin = no

<rules>
<rule>
importance = 100
condition  = eval(_CHR_ ne "target")
show       = no
</rule>
</rules>

</plot>

#<plot>
#file = %s
#show = conf(show_histogram)
#type = histogram
#r0 = 0.9r
#r1 = 0.945r
#color      = black
#fill_color = red,orange,yellow,green,blue
#fill_under = yes
#thickness  = 1
#sort_bin_values = no
#extend_bin = no
#
#<rules>
#<rule>
#importance = 100
#condition  = eval(_CHR_ ne "target")
#show       = no
#</rule>
#</rules>
#
#</plot>
''' % (options.output_dir+targetid+os.sep+"nucleotide_abs.txt", options.output_dir+targetid+os.sep+"nucleotide_rel.txt")]
		
#---------------
def calculateChromosomeSpecifics(chromsizes, f_configuration):
	""" calculates the overall chromosome size """

	counter_bases = 0
	for chromsize in chromsizes:
		counter_bases += chromsize[1]

	f_configuration["conf"]["genomesize"] = counter_bases

#---------------
def OnTargetWedge(targetid, chrom, start, end, chromsizes, f_configuration):
	""" calculates the part of the chromosome (loi) that constitutes the on-target """
	
	counter_bases = 0.
	target_center = 0.
	target_pos = 0
	counter = 0

	# target region should span 20 degrees
	target_size = end - start  # in nucleotides

	if (not f_configuration["conf"].has_key("genomesize")):
		calculateChromosomeSpecifics(chromsizes, f_configuration)
	zoom = f_configuration["conf"]["genomesize"] * f_configuration["conf"]["target_radius"] / 360. * 1 / target_size 
	f_configuration["conf"]["zoom"] = zoom
	
	if (options.verbose):
		print "zoom factor for target size: %.2f" % (f_configuration["conf"]["zoom"])
	
	chromosomes = [] # holds all chromosome names
	chromosomes_order = ["^"] # hold chromosome order
	for chromsize in chromsizes:

		if (chromsize[0] != chrom):
			chromosomes += [chromsize[0]]
			counter_bases += chromsize[1]
			chromosomes_order += [chromsize[0]]
		else:
			chrom_break1 = math.floor(start/f_configuration["image"]["chromosomes_units"])
			chrom_break2 = math.ceil(end/f_configuration["image"]["chromosomes_units"])
						
			 # add chromosomal upstream region of target in large enough
			if (chrom_break1 > 0.):
				chromosomes += ["%s[a]:0-%.3f" % (chromsize[0], chrom_break1)]
				chromosomes_order += ["a"]

			chromosomes += ["target"]
			chromosomes_order += ["target"]		
							
			# add chromosomal downstream region of target in large enough
			if (chrom_break2 < chromsize[1]/f_configuration["image"]["chromosomes_units"]): 
				chromosomes += ["%s[b]:%.3f-)" % (chromsize[0], chrom_break2)]
				chromosomes_order += ["b"]		
				
			target_center = counter_bases + start + 0.5 * (end - start) * zoom
			counter_bases += start + (end - start) * zoom + (chromsize[1] - end) 
			target_pos = counter
		counter += 1
		
	f_configuration["targets"][targetid]["chromosomes"]       = ";".join(chromosomes)
	f_configuration["targets"][targetid]["chromosomes_scale"] = "target:%d" % (zoom) 
	f_configuration["targets"][targetid]["chromosomes_order"] = "%s" % (",".join(chromosomes_order)) 
	
	f_configuration["conf"]["zoomedgenomesize"]               = counter_bases

	r_fraction = ( -90 - target_center/counter_bases * 360.) 
	f_configuration["image"]["angle_offset"]                  = r_fraction % 360.

	if (options.verbose):
		print "rotate image by %.2f degree to align target at the top" % (r_fraction % 360)

#---------------	
def makeOffTargetPlots(targetid, maxbin, tracks_per_lane, f_configuration):
	""" create heatmaps for off-target tracks """
	
	stroke_thickness = 0
	if (options.offtargetResolution >= 1.):
		stroke_thickness = 1
	
	spans = sorted(list(f_configuration["tiles"][targetid]["region"].keys()))
	for span in spans:
	
		radius = f_configuration["tiles"][targetid]["region"][span]["radius"]

		track_spacing = 0.003
		#                avail size           spacing between tracks    number of tracks  
		track_height = (radius[1]-radius[0] - (tracks_per_lane+1) * track_spacing )/tracks_per_lane
		track_shift = track_spacing # no shift
		# draw inwards
		r1 = radius[1] - track_shift
		r0 = radius[1] - track_shift - track_height
		
		filename = f_configuration["tiles"][targetid]["region"][span]["offtargetfile"]
		
		if (options.verbose):
			print "Preparing off target tracks for on-target size %d" % (span)
	
		f_configuration[targetid]["plot"] += ['''
<plot>
	show             = conf(show_heatmaps)
	type             = heatmap
	file             = %s
	r0               = %.3fr
	r1               = %.3fr
	stroke_thickness = %d
	color            = orrd-9-seq
	max              = %.3f
	scale_log_base   = 5
</plot>
''' % (filename, r0, r1, stroke_thickness, maxbin)]
	
	
#---------------
def writeOfftargetData(targetid, off_targets, chromsizes, tracks_per_lane, f_configuration):
	''' writes the off-target data tracks ''' 
	
	if (options.verbose):
		print "Start off-target tracks processing"

	# get step size	
	if (not f_configuration["conf"].has_key("genomesize")):
		calculateChromosomeSpecifics(chromsizes, f_configuration)
	step_size = int(f_configuration["conf"]["genomesize"] * options.offtargetResolution / (360. - f_configuration["conf"]["target_radius"] )) 

	if (options.verbose):
		print "off-target stepsize: %d nts" % (step_size)

	# convert chromsizes into dict
	cs = {}
	for chrom in chromsizes:
		cs[chrom[0]] = chrom[1]

	# create file handle	
	file_handles = {}
	spans = sorted(list(f_configuration["tiles"][targetid]["region"].keys()))
	for span in spans:
		f_configuration["tiles"][targetid]["region"][span]["offtargetfile"] = options.output_dir+targetid+os.sep+"offtarget_"+(str(span))+".txt"
		file_handles[span] = open(f_configuration["tiles"][targetid]["region"][span]["offtargetfile"],"w")
	
	try:
		maxbin = 0.
		for chrom in off_targets.keys():
			for start in range(0, cs[chrom], step_size):
				bins = {} # count the number of off-targets for each span size
				for span in spans:
					bins[span] = 0.001; # pseudocount
					
				end = min(start + step_size, cs[chrom])
				
				for offtarget in find(start, end, off_targets[chrom]):
					offlength = offtarget[1] - offtarget[0]
					for l in range(spans[0], offlength+1):
						if (bins.has_key(l)):
							bins[l] += 1
				
				for span in spans:
					file_handles[span].write("%s\t%d\t%d\t%d\n" % (chrom, start, end, bins[span]))
				
					maxbin = max(bins[span], maxbin)
	except:
		traceback.print_exc()

	# close file handles
	for span in spans:
		file_handles[span].close()
		
	return maxbin
			
#---------------	
def makeAnnotationPlots(targetid, maxbin, annotations, annotation_index, tracks_per_lane, f_configuration):
	""" create heatmaps for off-target tracks intersecting annotations """
	
	annotation = annotations[annotation_index]
	stroke_thickness = 0
	if (options.offtargetResolution >= 1.):
		stroke_thickness = 1
	
	spans = sorted(list(f_configuration["tiles"][targetid]["region"].keys()))
	for span in spans:
	
		radius = f_configuration["tiles"][targetid]["region"][span]["radius"]

		track_spacing = 0.003
		#                avail size           spacing between tracks    number of tracks  
		track_height = (radius[1]-radius[0] - (tracks_per_lane+1) * track_spacing )/tracks_per_lane
		track_shift = (annotation_index+2) * track_spacing + track_height * (annotation_index+1) # +1 due to non-annotated track
		# draw inwards
		r1 = radius[1] - track_shift
		r0 = radius[1] - track_shift - track_height

		filename = f_configuration["tiles"][targetid]["region"][span][annotation+"file"]
		
		if (options.verbose):
			print "Preparing off target annotation tracks for on-target size %d" % (span)
	
		f_configuration[targetid]["plot"] += ['''
<plot>
	show   = conf(show_heatmaps)
	type   = heatmap
	file             = %s
	r0               = %.3fr
	r1               = %.3fr
	stroke_thickness = %d
	color            = %s
	max              = %.3f
	scale_log_base   = 5
</plot>
''' % (filename, r0, r1, stroke_thickness, annotation_palettes[annotation_index % len(annotation_palettes)], maxbin)]

#---------------
def writeAnnotationOverlapData(targetid, off_targets, chromsizes, annotations, annotation_index, tracks_per_lane, f_configuration):

	annotation = annotations[annotation_index]

	if (options.verbose):
		print "Start annotation processing: %s " % (annotation)

	# get step size	
	if (not f_configuration["conf"].has_key("genomesize")):
		calculateChromosomeSpecifics(chromsizes, f_configuration)
	step_size = int(f_configuration["conf"]["genomesize"] * options.offtargetResolution / (360. - f_configuration["conf"]["target_radius"] )) 

	if (options.verbose):
		print "annotation stepsize: %d nts" % (step_size)

	# convert chromsizes into dict
	cs = {}
	for chrom in chromsizes:
		cs[chrom[0]] = chrom[1]

	# create file handle	
	file_handles = {}
	spans = sorted(list(f_configuration["tiles"][targetid]["region"].keys()))
	for span in spans:
		f_configuration["tiles"][targetid]["region"][span][annotation+"file"] = options.output_dir+targetid+os.sep+"offtarget_"+(str(span))+"_"+annotation+".txt"
		file_handles[span] = open(f_configuration["tiles"][targetid]["region"][span][annotation+"file"],"w")
	
	try:
		for chrom in off_targets.keys():
			for start in range(0, cs[chrom], step_size):
				bins = {} # count the number of off-targets for each span size
				for span in spans:
					bins[span] = 0.001; # pseudocount
					
				end = min(start + step_size, cs[chrom])
				
				for offtarget in find(start, end, off_targets[chrom]):
					offlength = offtarget[1] - offtarget[0]
					for l in range(spans[0], offlength+1):
						if (bins.has_key(l)):
							bins[l] += offtarget[2][annotation_index]
				
				for span in spans:
					file_handles[span].write("%s\t%d\t%d\t%d\n" % (chrom, start, end, bins[span]))
				
	except:
		traceback.print_exc()
		print >> sys.stderr, "Ops"

	# close file handles
	for span in spans:
		file_handles[span].close()
	
#---------------
def makeSequenceTrack(targetid, f_configuration, filename):
	
	f_configuration[targetid]["plot"] += ['''
<plot>
show   = conf(show_text)
type  = text
file  = %s
color = black
label_font  = bold
r0    = dims(ideogram,radius_outer) + 0.03r
r1    = dims(ideogram,radius_outer) + 0.06r
label_size = 13p
padding    = 0r
rpadding   = 0.25r
label_rotate   = no
label_parallel = yes
</plot>

''' % (filename)]
	

#---------------
def writeSequenceData(targetid, output_sequence_w, output_sequence_c, f_configuration):
	""" output the sequence data for the target region """
	
	f_sequence = open(options.output_dir+targetid+os.sep+"primarytarget_sequence.txt","w")
	for i in range(len(output_sequence_w)):
		f_sequence.write("target\t%d\t%d\t%s\t\n" % (i,i+1,output_sequence_c[i])) # crick # order important!
		f_sequence.write("target\t%d\t%d\t%s\t\n" % (i,i+1,output_sequence_w[i])) # watson

	f_sequence.close()
	
	makeSequenceTrack(targetid, f_configuration, options.output_dir+targetid+os.sep+"primarytarget_sequence.txt")

#---------------
def onTargets(records, loiBed, offBed, targets, tfotargets, submatches, annotations, chromsizes, f_configuration):
	""" process the data tracks """
	
	maxbin = 0. # maximal entries per bin (to scale all plots to the same level)
	
	# process one on-target region/locus of interest (LOI) at a time
	for loi in targets.keys():
		# info for whole region
		bed = loiBed[loi]
		chrom = bed[0]
		cstart = int(bed[1])
		cstop = int(bed[2])
		clabel = bed[3]
		cstrand = bed[5]
		loi_fasta = records[loi].seq	

		# process individual eligible target sites in this region (TTSpool)
		for target in targets[loi]:
			tstart = int(target[1])
			tstop = int(target[2])
			targetid = target[3]
			tstrand = target[5]
			tsize = tstop-tstart
			
			# add target specifics to f_configuration
			f_configuration[targetid] = {}
			f_configuration[targetid]["highlight"] = []
			f_configuration[targetid]["plot"] = []

			# calcuate ontarget wedge
			OnTargetWedge(targetid, chrom, tstart, tstop, chromsizes, f_configuration)
			
			# get target data filehandle		
			f_ontarget = open(options.output_dir+targetid+os.sep+"primarytarget.txt","w")
			
			targetSeq = ""
			if (cstrand == "+"):
				targetSeq = loi_fasta[tstart-cstart: tstop-cstart]
			else:
				targetSeq = loi_fasta.reverse_complement()[tstart-cstart: tstop-cstart]
			
			# output error marks
			outputErrors(loi, targetid, chrom, tstart, tstop, tsize, cstrand, target[8], f_configuration)
	
			# process off-targets w/r/t to on-target
			positions = [{} for r in range(tsize)] # indicates which nucleotide positions contribute to certrain error types
			for pos in range(tsize):
				for t in "-dobt":
					positions[pos][t] = 0

			counter = 1
			off_targets = {} # keeping offtargets 
			for chromsize in chromsizes:
				off_targets[chromsize[0]] = IntervalNode( -1, -1, -1 )
				
			for loc in tfotargets[target[3]]:
				locOffset1 = loc[1]
				locOffset2 = loc[2]
				jsonOffset = 0
				locStr = ["-" for i in range(locOffset1,locOffset2)]
				# add error encoding
				errors = loc[8]
				if (errors!=""):
					i = 0
					while (i<len(errors)):
						if (errors[i] in "dobt"):
							# get number
							j=i+1
							while(j<len(errors) and errors[j] in "0123456789"):
								j += 1 
							locStr[int(errors[i+1:j])]=errors[i] # requires the error to be with respect to the TFO
							i=j
				locFullStr = [" " for i in range(tsize)]
				if (tstrand=='+'):
					jsonOffset = locOffset1
					for i in range(len(locStr)):
						locFullStr[locOffset1+i] = locStr[i];
				else:
					jsonOffset = tsize-locOffset1-len(locStr)
					for i in range(len(locStr)):
						locFullStr[tsize-locOffset1-i-1] = locStr[i];
				
				# add nucleotide interactions for stacked barchart
				for pos in range(tsize):
					if (locFullStr[pos] != " "):
						positions[pos][locFullStr[pos]] += 1	
				
				locChr = offBed[loc[3]][0]
				locStart = int(offBed[loc[3]][1]) + int(loc[4])
				locStop = int(offBed[loc[3]][1]) + int(loc[5])

				annoarray = [ 0 for i in range(len(annotations)) ]
				for i in range(len(annotations)):
					value = 0
					if (loc[12+i] != "-"):
						value = len(loc[12+i].split(","))
					annoarray[i] = value

				# add offtarget to interval tree
				off_targets[locChr] = off_targets[locChr].insert(locStart, locStop, annoarray)		
				counter += 1 
			
			# output positions
			writeNuclotidePositionData(loi, targetid, chrom, tstart, tstop, tsize, cstrand, positions, f_configuration)
			
			output_sequence_w = targetSeq.upper()
			output_sequence_c = targetSeq.complement().upper()
			writeSequenceData(targetid, output_sequence_w, output_sequence_c, f_configuration)

			# output all submatches of the primary target
			for submatch in submatches[targetid]:
				f_ontarget.write("target\t%d\t%d\tvalue=%d,url=inspector_detail.html?rId=%s&region=%d-%d\n" % (submatch[0], submatch[1], submatch[2], targetid, submatch[0], submatch[1]))
		
			# add tile configuration
			writeTileData(submatches[targetid], targetid, options.output_dir+targetid+os.sep+"primarytarget.txt", f_configuration)
				
			tracks_per_lane = 1 + len(annotations)
			# output off-target data
			maxbin = max(maxbin, writeOfftargetData(targetid, off_targets, chromsizes, tracks_per_lane, f_configuration))
			
			
			for annotation_index in range(len(annotations)):
				writeAnnotationOverlapData(targetid, off_targets, chromsizes, annotations, annotation_index, tracks_per_lane, f_configuration)
			
		f_ontarget.close()
	
	# output circos tracks
	# process one on-target region/locus of interest (LOI) at a time
	for loi in targets.keys():
		# process individual eligible target sites in this region (TTSpool)
		for target in targets[loi]:
			tstart = int(target[1])
			tstop = int(target[2])
			targetid = target[3]
			
			# output off-target tracks
			makeOffTargetPlots(targetid, maxbin, tracks_per_lane, f_configuration)

			# output annotations
			for annotation_index in range(len(annotations)):
				makeAnnotationPlots(targetid, maxbin, annotations, annotation_index, tracks_per_lane, f_configuration)


#---------------
def writeConfiguration(f_configuration):
	""" write circos configuration file """
	
	# open filehandle
	for targetid in f_configuration["targets"].keys():

		foutput = open(options.output_dir+targetid+os.sep+"circos.conf","w")

		foutput.write('''
show_links        = no
show_highlights   = yes
show_text         = yes
show_heatmaps     = yes
show_scatter      = no
show_histogram    = yes
show_tiles        = yes
show_ticks        = yes
show_tick_labels  = yes
show_grid         = yes
		

<colors>
<<include etc/colors.conf>>
red_faint   = 255,0,0,0.8
green_faint = 0,255,0,0.8
blue_faint  = 0,0,255,0.8
grey_faint  = 200,200,200,0.5
green_triad = 200,233,192
red_nontriad = 253,212,158
red_nontriad_faint = 253,212,158,0.8

</colors>
<fonts>
<<include etc/fonts.conf>>
</fonts>

<patterns>
<<include etc/patterns.conf>>
</patterns>

<<include etc/housekeeping.conf>>

<<include %sideogram.conf>>
<<include %s>>

karyotype = %s

<image>
dir = %s
file  = %s.png
# radius of inscribed circle in image
png   = yes
svg   = no # set to yes to generate vector graphics
# radius of inscribed circle in image
radius         = %dp

# by default angle=0 is at 3 o'clock position
angle_offset   = %.1f

auto_alpha_colors = yes
auto_alpha_steps  = 5

background     = white
24bit = yes

image_map_use      = yes
image_map_name     = %s_map
</image>

chromosomes_units = %.0f
chromosomes_display_default = no
''' % (options.output_dir, f_configuration["ticks"][targetid], f_configuration["targets"][targetid]["karyotype"], options.output_dir, targetid, f_configuration["image"]["radius"], f_configuration["image"]["angle_offset"], targetid, f_configuration["image"]["chromosomes_units"]))

		foutput.write("chromosomes = %s\n" % f_configuration["targets"][targetid]["chromosomes"] )
		foutput.write("chromosomes_order = %s\n" % f_configuration["targets"][targetid]["chromosomes_order"] )
		foutput.write("chromosomes_scale = %s\n" % f_configuration["targets"][targetid]["chromosomes_scale"] )
		
		# add plots
		if (f_configuration.has_key("plots") or f_configuration[targetid].has_key("plots")):
			foutput.write("<plots>\n")
			for key in f_configuration["plots"]:
				if (key!="plot"):
					foutput.write("%s = %s\n" % (key, f_configuration["plots"][key]))				
					
			if (f_configuration["plots"].has_key("plot")):
				for plot in f_configuration["plots"]["plot"]:
					foutput.write(plot)

			if (f_configuration[targetid].has_key("plot")):
				for plot in f_configuration[targetid]["plot"]:
					foutput.write(plot)
					
			
			foutput.write("</plots>\n")
	
		# add hightlights
		if (f_configuration.has_key("highlights") or f_configuration[targetid].has_key("highlights")):
			foutput.write("<highlights>\n")
			for key in f_configuration["highlights"]:
				if (key!="highlight"):
					foutput.write("%s = %s\n" % (key, f_configuration["highlights"][key]))				
					
			if (f_configuration["highlights"].has_key("highlight")):
				for highlight in f_configuration["highlights"]["highlight"]:
					foutput.write(highlight)

			if (f_configuration[targetid].has_key("highlight")):
				for highlight in f_configuration[targetid]["highlight"]:
					foutput.write(highlight)

			foutput.write("</highlights>\n")

		
		# close filehandle
		foutput.close()

#---------------
def makeIdeogram(targets, chromsizes, f_configuration):
	
	# process one on-target region/locus of interest (LOI) at a time
	for loi in targets.keys():

		# process individual eligible target sites in this region (TTSpool)
		for target in targets[loi]:
			tstart = int(target[1])
			tstop = int(target[2])
			targetid = target[3]
			tsize = tstop-tstart
			
			f_karyotyp = open(options.output_dir+targetid+os.sep+"karyotyp.txt","w")
			for chrom in chromsizes:
				# check color handling
				colorId = chrom[0].split('_')[0].lower()
				if (not colorId.startswith("chr") or len(colorId)>5):
					colorId = "chrun"
				f_karyotyp.write("chr - %s %s 0 %d %s\n" % (chrom[0], chrom[0], chrom[1], colorId))

			# add target ideogram
			f_karyotyp.write("chr - target target 0 %d green_triad url=report.html\n" % (tsize))
			f_configuration["targets"][targetid]["karyotype"] = options.output_dir+targetid+os.sep+"karyotyp.txt"
	
			f_karyotyp.close()
	
	f_ideogram_conf = open(options.output_dir+"ideogram.conf","w")
	f_ideogram_conf.write(
'''<ideogram>

<spacing>
default = 0.002r
break   = 0.05r
axis_break_at_edge = yes
axis_break         = yes
axis_break_style   = 2

<break_style 2>
stroke_color     = black
stroke_thickness = 5p
thickness        = 2r
</break>

</spacing>

# ideogram positions
radius           = 0.94r
thickness        = 30p
fill             = yes
fill_color       = black
stroke_thickness = 2
stroke_color     = black

# ideogram labels
show_label       = yes
label_font       = default
# labels outside circle
label_radius     = dims(ideogram,radius) + 0.005r
label_with_tag   = no
label_size       = 20p
label_parallel   = yes
label_case       = upper

#ideogram bands
show_bands            = yes
fill_bands            = yes
band_stroke_thickness = 2
band_stroke_color     = white
band_transparency     = 4
</ideogram>''')
	f_ideogram_conf.close()
	
	f_configuration["chromosomes_order_by_karyotype"] = "yes"
	f_configuration["chromosomes_units"] = 1000000
	f_configuration["chromosomes_display_default"] = "yes"

#---------------
def createSubdirectories(targets, f_configuration):

	# process one on-target region/locus of interest (LOI) at a time
	for loi in targets.keys():

		# process individual eligible target sites in this region (TTSpool)
		for target in targets[loi]:
			targetid = target[3]

			subdir = options.output_dir+targetid+os.sep
			
			if not os.path.exists(subdir):
				os.makedirs(subdir)

			if (not f_configuration.has_key("targets")):
				f_configuration["targets"] = {}
			
			if (not f_configuration["targets"].has_key(targetid)):
				f_configuration["targets"][targetid] = {}



#---------------
def writeTicks(targetid, f_configuration):
	''' create chromosomal ticks '''

	if (not f_configuration.has_key("ticks")):
		f_configuration["ticks"] = {}
		
	f_configuration["ticks"][targetid] = options.output_dir+targetid+os.sep+"ticks.conf"

	f_ticks = open(f_configuration["ticks"][targetid],"w")
	f_ticks.write('''
<ticks>
show             = conf(show_ticks)
tick_label_font  = light
radius           = dims(ideogram,radius_outer) + 0.025r
label_offset     = 5p
label_size       = 16p
multiplier       = 1e-6
color            = black
thickness        = 1p

# 25 Mb ticks, all chromosomes, with labels
<tick>
spacing        = 25u
size           = 12p
show_label     = yes
format         = %d
min_label_distance_to_edge = 10p
</tick>

# 5 Mb ticks, all chromosomes, with labels
# labels must be separated by at least 1px, which
#   avoids overlap on human chrs
<tick>
label_separation = 1p
spacing          = 5u
size             = 7p
show_label       = yes
label_size       = 10p
format           = %d
min_label_distance_to_edge = 10p
</tick>

<tick>
use            = no
spacing_type   = relative
rspacing       = 0.999
size           = 6p
grid_start     = 0.21r
grid_end       = dims(ideogram,radius_outer) + 0.025r
grid_color     = grey
grid_thickness = 1p
grid           = yes
</tick>

# 2% & 10% relative ticks on human chromosomes
<tick>
radius         = 0.915r
spacing_type   = relative
rspacing       = 0.02
size           = 2p
show_label     = no
color          = grey
</tick>

<tick>
label_separation = 1p
label_size     = 10p
radius         = 0.915r
spacing_type   = relative
rspacing       = 0.1
size           = 4p
show_label     = yes
label_relative = yes
rmultiplier    = 100
format         = %d
suffix         = %
</tick>

<tick>
radius         = 0.2r
spacing_type   = relative
rspacing       = 0.02
size           = 2p
show_label     = no
color          = grey
</tick>

''')

	f_ticks.write('''
</ticks>
''')
	f_ticks.close()

#---------------
def initImage(f_configuration):

	global annotation_palettes

	f_configuration["image"] = {}
	f_configuration["image"]["radius"] = 1500
	f_configuration["image"]["angle_offset"] = -90
	f_configuration["image"]["chromosomes_units"] = 1000000.

	f_configuration["conf"] = {}
	f_configuration["conf"]["target_radius"] = 20 # target region in degree
	
	annotation_palettes = ["bupu-9-seq","greens-9-seq","pubu-9-seq"]
	
	# add to configuration file
	if (not f_configuration.has_key("highlights")):
		f_configuration["highlights"] = {}

	f_configuration["highlights"]["z"] = "0"
	f_configuration["highlights"]["fill_color"] = "green"

	if (not f_configuration["highlights"].has_key("highlight")):
		f_configuration["highlights"]["highlight"] = []

	# add common plots
	if (not f_configuration.has_key("plots")):
		f_configuration["plots"] = {}

	# write background
	f_background_conf = open(options.output_dir+"background.conf","w")
	f_background_conf.write('''
<backgrounds>
show = yes
<background>
color = grey_faint
</background>
</backgrounds>
''')
	f_background_conf.close()
	
#---------------
def process():
	# configuration data container
	f_configuration = {}
	
	#define the image parameters
	initImage(f_configuration)
	
	# get the data
	(records, loiBed, offBed, targets, tfotargets, annotations, submatches, chromsizes) = readdata()

	# get genome size
	calculateChromosomeSpecifics(chromsizes, f_configuration)
	
	# -----------------------------
	# create data shared between targets
	
	# process chromatin data
	if (options.chromatin != ""):
		processChromatin(offBed, f_configuration)

	# -----------------------------
	# create data individual to targets
	createSubdirectories(targets, f_configuration)
		
	# make ideogram
	makeIdeogram(targets, chromsizes, f_configuration)

	# add ontarget tracks
	onTargets(records, loiBed, offBed, targets, tfotargets, submatches, annotations, chromsizes, f_configuration)

	# write configuration file
	writeConfiguration(f_configuration)
	
			
# manage option and arguments processing
def main():
	global options
	global args
	usage = '''usage: %prog [options] LOI.bed LOC.bed TTSpool.file TRIPLEX.file primary.targets.json chromosome.defs

Generates data files for circos plot generation
	'''
	
	parser = OptionParser(usage)
	parser.add_option("-q", "--quiet", action="store_false", dest="verbose", default=True,
					help="don't print status messages to stdout")
	parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
					help="print status messages to stdout")
	parser.add_option("-s", "--fasta-file", type="string", dest="fastaFile", default="", 
					help="location with fasta sequences (overwrites any specified column)")
	parser.add_option("-c", "--chromatin", type="string", dest="chromatin", default="", 
					help="location of a chromatin data")
	parser.add_option("-F", "--chromatin-format", type="string", dest="chromatinFormat", default="", 
					help="either bam (requires pysam) or bigwig (requires bx-python)")
	parser.add_option("-C", "--chromatin-resolution", type="float", dest="chromatinResolution", default=0.5, 
					help="chromatin resolution in degree, minimum = 0.01, maximum = 10")
	parser.add_option("-O", "--offtarget-resolution", type="float", dest="offtargetResolution", default=0.5, 
					help="offtarget resolution in degree, minimum = 0.01, maximum = 10")
#	parser.add_option("-a", "--annotation", type="string", dest="annotation", default="", 
#					help="location of the annotation data")
	parser.add_option("-o", "--output-dir", type="string", dest="output_dir", default="", 
					help="directory where the output files will be saved into")

	(options, args) = parser.parse_args()
	if (len(args) != 6):
		parser.print_help()
		parser.error("incorrect number of arguments")

	if (options.chromatin != "" and options.chromatinResolution < 0.01 or  options.chromatinResolution > 10):
		parser.error("chromatin resolution outside valid interval : 0.01 <= x <= 10; was %.3f" % options.chromatinResolution)
		
	if (options.chromatin != "" and not (options.chromatinFormat == "bam" or options.chromatinFormat == "bigwig")):
		parser.error("chromatin filetype not supported: %s" % options.chromatinFormat)		

	if (options.offtargetResolution < 0.01 or  options.offtargetResolution > 10):
		parser.error("offtarget resolution outside valid interval : 0.01 <= x <= 10; was %.3f" % options.offtargetResolution)
	
	# convert output dir to absolute path to avoid circos conflicts
	if (options.output_dir != ""):
		options.output_dir = os.path.abspath(options.output_dir)+os.sep
	else:
		options.output_dir = os.path.abspath(".")+os.sep

	print options.output_dir
	
	if (options.verbose):
		print >> sys.stderr, "collect Output from \n(LOI.bed) : %s\n(LOC.bed) : %s\n(TTSpool.file) : %s\n(TRIPLEX.file) : %s" % (args[0],args[1],args[2],args[3])
		if (options.chromatin != ""):
			print >> sys.stderr, "chromatin data : %s \nformat : %s\n" % (options.chromatin, options.chromatinFormat)			

	process()

	
if __name__ == "__main__":
	main()
	
