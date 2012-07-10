#!/usr/bin/python

import sys, os, glob
from optparse import OptionParser
import fileinput

global options
global args

# read all the required data
def readdata():

	# the primary target site bed file provides length information about the primary sites
	tts = {}
	for line in fileinput.input([args[0]]):
		cols = line.strip().split("\t")
		if (len(cols)<6):
			continue
		try:
			tts[cols[3]] = {'chrom':cols[0], 'start':int(cols[1]), 'end':int(cols[2]), 'strand':cols[5]}
		except:
			if (options.verbose):
				print >> sys.stderr, 'skipping line TTS.file %s' % (line)
	
	# the tpx file provides all information about the off targets
	tpx = {}
	signature = [] # the signature provides information about the annotation used
	for line in fileinput.input([args[1]]):
		cols = []
		# keep signature
		if (line.startswith("#")):
			signature = line.strip("#").strip().split("\t")
			continue
		
		cols = line.strip().split("\t")
		if (len(cols)<6):
			if (options.verbose):
				print >> sys.stderr, 'skipping line %s' % (line)
			continue
		
		ht = {}	
		for i in range(len(cols)):
			key = signature[i]
			value = ""
			try:
				value = int(cols[i])
			except:
				try:
					value = float(cols[i])
				except:
					value = cols[i]
			ht[key]=value
			
		oId = ht["Sequence-ID"]
		score = ht["Score"]

		# check if Sequence-ID is known
		if (not tts.has_key(oId)):
			print >> sys.stderr, 'Sequence-ID not found in TTS-file %s' % (line)
			continue
			
		# register entry
		if (not tpx.has_key(oId)):
			tpx[oId] = {}
			
		# categorize with score
		if (not tpx[oId].has_key(score)):
			tpx[oId][score] = []
		tpx[oId][score] += [ht]
		
	if (options.verbose):
		print >> sys.stderr, "...finished"

	return (signature, tpx, tts)
	

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

def generate_data_sets(signature, tpx, tts):
	""" 
	generate tab-separated data files for R processing
	"""
	
	annotations = signature[12:]

	# process each primary target on its own
	for oId in tpx.keys():
		scores_ht = tpx[oId]
		scores = scores_ht.keys()
		scores.sort()
		tlength = tts[oId]['end']-tts[oId]['start']
		
		# data format for triads plot
		outp = open(options.output+oId+"_data.tsv","w")
		outp.write("%s\t%s\t%s\n" % ("Triads", "Value", "Type"))
		for score in scores:
			features = scores_ht[score]
			outp.write("%d\t%d\t%s\n" % (score, len(features), "-"))
			
			for annotation in annotations:
				value = 0
				for feature in features:
					if (feature[annotation] !="-"):
						value += len(feature[annotation].split(","))
				outp.write("%d\t%d\t%s\n" % (score, value, annotation))
	
		outp.close()
		
	annotations += ["-"]
	# process each primary target on its own
	for oId in tpx.keys():
		scores_ht = tpx[oId]
		scores = scores_ht.keys()
		scores.sort()
		tlength = tts[oId]['end']-tts[oId]['start']
		tstrand = tts[oId]['strand']
		
		# data format for error plots
		outp = open(options.output+oId+"_errors.tsv","w")
		outp.write("%s\t%s\t%s\t%s\t%s\n" % ("Triads", "Position", "Value", "Error", "Annotation"))
		# iterate over all triads categories (score)
		for score in scores:
			features = scores_ht[score]
			for annotation in annotations:
				# create array to register errors
				positions = _make_position_array(features, oId, tlength, tstrand, annotation)
				for t in positions.keys():
					for pos in range(tlength):
						if (positions[t][pos]>0):
							# shift position to start at 1 rather than zero for R plotting (add 1)
							if (t == '-'):
								outp.write("%d\t%d\t%d\tnone\t%s\n" % (score, pos+1, positions[t][pos], annotation))
							elif (t =='o'):
								outp.write("%d\t%d\t%d\toligo\t%s\n" % (score, pos+1, positions[t][pos], annotation))
							elif (t =='b'):
								outp.write("%d\t%d\t%d\tboth\t%s\n" % (score, pos+1, positions[t][pos], annotation))
							elif (t =='t'):
								outp.write("%d\t%d\t%d\ttriplex\t%s\n" % (score, pos+1, positions[t][pos], annotation))
							elif (t =='d'):
								outp.write("%d\t%d\t%d\tduplex\t%s\n" % (score, pos+1, positions[t][pos], annotation))

		outp.close()

def plot_targets_per_triad(signature, tpx, tts):
	""" 
	create histogram with the number of off-targets stabilized by a given number of triads
	"""

	# process each primary target separately
	for oId in tpx.keys():
		tlength = tts[oId]['end']-tts[oId]['start']
		
		#write Rscript
		outp = open(options.output+oId+"_triads.R","w")
		rscript = '''
			library(ggplot2)
			odf = read.delim("%s")
			for (i in levels(odf$Type)){
				if (sum(odf[which(odf$Type==i),]$Value)==0){
					odf[which(odf$Type==i),]$Value=0.0001
				}
			}
			labels = seq(min(odf$Triads),%d)
		''' % (options.output+oId+"_data.tsv", tlength)
		# write pdf if requested
		if (options.format=="both" or options.format=="pdf"):
			rscript += '''
				pdf("%s", width=10, height=5, bg="white")
				ggplot(odf)+aes(x=Triads, weight=Value, fill = ..count..)+geom_histogram(binwidth=1)+facet_grid(.~Type) + scale_y_continuous('off-targets') + scale_fill_gradient("Count", low = "green", high = "red") + xlim(min(odf$Triads), %d) + ylim(0,max(1,odf$Value)) + scale_x_continuous(breaks=labels-0.5, labels=labels)
				dev.off()
			''' % (options.output+oId+"_triads.pdf", tlength)
			
		# write png if requested
		if (options.format=="both" or options.format=="png"):
			rscript += '''
				png(filename = "%s", width=1000, height=550, units="px", pointsize=12, bg="white")
				ggplot(odf)+aes(x=Triads, weight=Value, fill = ..count..)+geom_histogram(binwidth=1)+facet_grid(.~Type) + scale_y_continuous('off-targets') + scale_fill_gradient("Count", low = "green", high = "red") + xlim(min(odf$Triads), %d) + ylim(0,max(1,odf$Value)) + scale_x_continuous(breaks=labels-0.5, labels=labels) + opts(title="Off-targets expected to form triplexes stabilized by a given number of nucleotide triads.")
				dev.off()
			''' % (options.output+oId+"_triads.png", tlength)
	
		outp.write(rscript)
		outp.close()
				
def plot_tfo_contribution_per_annotation(signature, tpx, tts):
	""" 
	create a plot that shows the position of the tfo contributing to this off-target
	subject to the annotation classes
	"""
	
	annotations = ["-"]+signature[12:]

	# process each primary target on its own
	for oId in tpx.keys():
		tlength = tts[oId]['end']-tts[oId]['start']
	
		#write Rscript
		outp = open(options.output+oId+"_errors_annotation.R","w")
		rscript = '''
			library(ggplot2)
			odf = read.delim("%s")
			cbbFillPalette = scale_fill_manual(values=c("#0072B2", "#E69F00", "#56B4E9", "#009E73", "#999999", "#D55E00", "#CC79A7"))
			odf$Error <- factor(odf$Error, levels =c("duplex","oligo","both","triplex","none"))
			max_y = 1
			a = odf[which(odf$Annotation=="-"),]
			for (i in names(table(a$Position))){
				max_y = max(max_y,sum(a[which(a$Position==i),]$Value))
			}

			for (i in levels(odf$Annotation)){
				if (sum(odf[which(odf$Annotation==i),]$Value)==0){
					odf[which(odf$Annotation==i),]$Value=0.0001
				}
			}
		'''	% (options.output+oId+"_errors.tsv")
	
		# write pdf if requested
		if (options.format=="both" or options.format=="pdf"):
			rscript += '''
				pdf("%s", width=10, height=5, bg="white")
				ggplot(odf)+aes(x=Position, fill=Error, weight=Value)+geom_bar(position="stack", binwidth=1)+facet_grid(.~Annotation) + cbbFillPalette + xlim(0, %d) + ylim(0,max_y)
				dev.off()
			''' % (options.output+oId+"_errors_annotation.pdf", tlength)
			
		# write png if requested
		if (options.format=="both" or options.format=="png"):
			rscript += '''
				png(filename = "%s", width=1000, height=550, units="px", pointsize=12, bg="white")
				ggplot(odf)+aes(x=Position, fill=Error, weight=Value)+geom_bar(position="stack", binwidth=1)+facet_grid(.~Annotation) + cbbFillPalette + xlim(0, %d) + ylim(0,max_y) + opts(title = "Contribution of nucleotide position in the primary target (x-axis) to off-target effects (y-axis), categorised by annotation")
				dev.off()
			''' % (options.output+oId+"_errors_annotation.png", tlength)

		outp.write(rscript)
		outp.close()


def plot_tfo_contribution_per_triad(signature, tpx, tts):
	"""
	create a plot that shows the position of the tfo contributing to this off-target
	subject to the number of triads stabilizing the triplex
	"""
	# process each primary target on its own
	for oId in tpx.keys():
		scores_ht = tpx[oId]
		scores = scores_ht.keys()
		scores.sort()
		tlength = tts[oId]['end']-tts[oId]['start']

		#write Rscript
		outp = open(options.output+oId+"_errors_triads.R","w")
		rscript = '''
			library(ggplot2)
			odf = read.delim("%s")
			cbbFillPalette = scale_fill_manual(values=c("#0072B2", "#E69F00", "#56B4E9", "#009E73", "#999999", "#D55E00", "#CC79A7"))
			odf = odf[which(odf$Annotation=="-"),]
			odf$Error <- factor(odf$Error, levels =c("duplex","oligo","both","triplex","none"))
			max_y = 1
			for (j in names(table(odf$Triads))){
				a = odf[which(odf$Triads==j),]
				for (i in names(table(a$Position))){
					max_y = max(max_y,sum(a[which(a$Position==i),]$Value))
				}
			}
		'''	% (options.output+oId+"_errors.tsv")
	
		# write pdf if requested
		if (options.format=="both" or options.format=="pdf"):
			rscript += '''
				pdf("%s", width=10, height=5, bg="white")
				ggplot(odf)+aes(x=Position, fill=Error, weight=Value)+geom_bar(position="stack", binwidth=1)+facet_grid(.~Triads) + cbbFillPalette + xlim(0, %d) + opts(axis.text.x=theme_text(size=10, angle=45, hjust=1)) + ylim(0,max_y)
				dev.off()
			''' % (options.output+oId+"_errors_triads.pdf", tlength)
			
		# write png if requested
		if (options.format=="both" or options.format=="png"):
			rscript += '''
				png(filename = "%s", width=1000, height=550, units="px", pointsize=12, bg="white")
				ggplot(odf)+aes(x=Position, fill=Error, weight=Value)+geom_bar(position="stack", binwidth=1)+facet_grid(.~Triads) + cbbFillPalette + xlim(0, %d) + opts(axis.text.x=theme_text(size=10, angle=45, hjust=0)) + ylim(0,max_y) + opts(title = "Contribution of nucleotide position in the primary target (x-axis) to off-target effects (y-axis), binned by triplex length.")
				dev.off()
			''' % (options.output+oId+"_errors_triads.png", tlength)

		outp.write(rscript)
		outp.close()


def process():
	"""
	workflow
	"""
	
	(signature, tpx, tts) = readdata()
	generate_data_sets(signature, tpx, tts)
	plot_targets_per_triad(signature, tpx, tts)
	plot_tfo_contribution_per_annotation(signature, tpx, tts)
	plot_tfo_contribution_per_triad(signature, tpx, tts)
	
# manage option and arguments processing
def main():
	global options
	global args
	usage = '''usage: %prog [options] TTS.bed TPX.file

takes a bed file containing the location of the primary targets 
and a triplexator output file containing information about the 
off-targets and generates files that can be used with R
to generate plots
	'''
	parser = OptionParser(usage)
	parser.add_option("-q", "--quiet", action="store_false", dest="verbose", default=True,
					help="don't print status messages to stdout")
	parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
					help="print status messages to stdout")
	parser.add_option("-o", "--output", type="string", dest="output", default="", 
					help="where to write the output files")
	parser.add_option("-f", "--format", type="string", dest="format", default="both", 
					help="image format to produce: pdf, png or both")
	
	(options, args) = parser.parse_args()
	if (len(args) != 2):
		parser.print_help()
		parser.error("incorrect number of arguments")
		
	if (options.output != ""): 
		options.output += "/"
		
	if (options.verbose):
		print >> sys.stderr, "Generate plots for %s (TTS.bed) with results %s (TPX.file)" % (args[0], args[1])
			
	process()

	
if __name__ == "__main__":
	main()
